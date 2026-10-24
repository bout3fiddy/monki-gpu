"""Shared GPU preparation, dispatch grouping, ownership and float64 round reductions."""

import ctypes
import hashlib
import platform
from pathlib import Path

import numpy as np

from ..input import ReferenceCase
from ..output import SimulationResult
from ..reference import file_hash
from .backend import Backend
from .cpu import CPU
from .gpu_source import specialized_source
from .gpu_tables import (
    LAYOUT_VERSION,
    TALLY_SHARDS,
    pack_tables,
)
from .tallies import BatchControls, RawTallies

ROOT = Path(__file__).resolve().parents[3]
FLOAT_POINTER = ctypes.POINTER(ctypes.c_float)
UINT_POINTER = ctypes.POINTER(ctypes.c_uint32)


class GPUBackend(Backend):
    """Retain GPU buffers across rounds; all photon histories execute on the GPU.

    Only setup uses the Fortran CPU worker, without running CPU photons. GPU
    transport is float32; host round accumulation is float64. Use one caller per
    instance. A failed dispatch discards the session and its partial work.
    """

    backend_name = "GPU"
    native_prefix = "monki_gpu"
    required_platform = ""
    default_library = ROOT / "build/gpu/library"
    default_kernel = ROOT / "kernels/gpu/transport"

    def __init__(
        self,
        device: int = 0,
        *,
        save_path: bool = True,
        profile: bool = False,
        seed: int = 23274,
        photons_per_dispatch: int = 65536,
        tally_shards: int = 32,
        max_steps: int = 100_000,
        library: Path | None = None,
        kernel: Path | None = None,
        upstream: Path = ROOT / "fortran/monki/upstream",
    ):
        for name, value, minimum, maximum in (
            ("device", device, 0, 2**32 - 1),
            ("seed", seed, 0, 2**32 - 1),
            ("tally_shards", tally_shards, 1, 64),
            ("photons_per_dispatch", photons_per_dispatch, 1, 65536),
            ("max_steps", max_steps, 1, 100_000),
        ):
            if type(value) is not int or not minimum <= value <= maximum:
                raise ValueError(f"{name} must be an integer in [{minimum}, {maximum}]")

        if type(save_path) is not bool:
            raise ValueError("save_path must be a bool")

        if type(profile) is not bool:
            raise ValueError("profile must be a bool")

        self.profile = profile
        self.device = device
        self.save_path = save_path
        self.seed = seed
        self.photons_per_dispatch = photons_per_dispatch
        self.tally_shards = tally_shards
        self.max_steps = max_steps
        self.library = library or self.default_library
        self.kernel = kernel or self.default_kernel
        self.upstream = upstream
        self.metadata: SimulationResult | None = None
        self.provenance: dict = {}
        self._native = None
        self._handle = None
        self._shape = (0, 0, 0, 0, 0)
        self._configuration = None
        self._round = 0
        self._tally_count = 0
        self._round_capacity = 1

    def prepare(self, scene: object) -> None:
        """Reuse Victor's optical setup, compile the kernel, and upload resident tables."""

        self.close()

        if not isinstance(scene, ReferenceCase):
            raise TypeError(f"{self.backend_name}.prepare expects a ReferenceCase")

        if scene.backward:
            raise NotImplementedError(f"{self.backend_name} backward transport is not implemented")

        if platform.system() != self.required_platform:
            raise RuntimeError(
                f"{self.backend_name} requires {self.required_platform} and a compatible GPU"
            )

        if not self.library.is_file():
            raise FileNotFoundError(
                f"Build {self.backend_name} with uv run python "
                f"scripts/build_{self.backend_name.lower()}.py"
            )

        # The preparer shares the original input, wavelength interpolation, phase
        # tables and metadata. No run_batch call is made on this CPU instance.
        with CPU(save_path=self.save_path, upstream=self.upstream) as preparer:
            preparer.prepare(scene)
            packed, self._shape = pack_tables(preparer.transport_tables())
            metadata = preparer.metadata
            preparation = preparer.provenance.copy()

        # Independent tally banks reduce contention without changing photon histories.
        packed[TALLY_SHARDS] = self.tally_shards
        nmat, nx, ny, nz, angles = self._shape
        self._tally_count = nmat * nx * ny * angles

        if self.save_path:
            self._tally_count += nx * ny * (nz + 1) * angles

        # Small rounds otherwise underfill the GPU. Reserve bounded, independent
        # output slots; no round is computed until explicitly requested.
        bytes_per_round = self._tally_count * self.tally_shards * 4
        self._round_capacity = max(
            1,
            min(
                32,
                self.photons_per_dispatch // (scene.photons_per_column * nx * ny),
                (64 * 1024 * 1024) // bytes_per_round,
            ),
        )
        native = ctypes.CDLL(str(self.library))
        abi_version = getattr(native, f"{self.native_prefix}_abi_version", None)

        if abi_version is None or abi_version() != 2:
            raise RuntimeError(
                f"Rebuild the {self.backend_name} bridge with uv run python "
                f"scripts/build_{self.backend_name.lower()}.py"
            )

        for operation in ("create", "run", "gpu_seconds", "close", "device_name"):
            setattr(
                native,
                f"monki_gpu_{operation}",
                getattr(native, f"{self.native_prefix}_{operation}"),
            )

        native.monki_gpu_create.argtypes = [
            ctypes.c_char_p,
            ctypes.c_uint32,
            FLOAT_POINTER,
            ctypes.c_size_t,
            ctypes.c_size_t,
            ctypes.c_uint32,
            ctypes.c_char_p,
            ctypes.c_size_t,
        ]
        native.monki_gpu_create.restype = ctypes.c_void_p
        native.monki_gpu_run.argtypes = [
            ctypes.c_void_p,
            UINT_POINTER,
            FLOAT_POINTER,
            ctypes.c_char_p,
            ctypes.c_size_t,
        ]
        native.monki_gpu_run.restype = ctypes.c_int
        native.monki_gpu_gpu_seconds.argtypes = [ctypes.c_void_p]
        native.monki_gpu_gpu_seconds.restype = ctypes.c_double
        native.monki_gpu_close.argtypes = [ctypes.c_void_p]
        native.monki_gpu_close.restype = None
        native.monki_gpu_device_name.argtypes = [ctypes.c_void_p]
        native.monki_gpu_device_name.restype = ctypes.c_char_p
        error = ctypes.create_string_buffer(8192)

        # Hash exactly the source sent to the compiler, even if the file is edited
        # while a prepared session remains alive.
        # Scene dimensions and estimator modes stay fixed for this session. Giving
        # the compiler those constants removes inactive paths and dynamic indexing.
        source = specialized_source(self.kernel, packed, nz, angles)
        bridge_hash = file_hash(self.library)
        self._native = native
        self._handle = native.monki_gpu_create(
            source,
            self.device,
            packed.ctypes.data_as(FLOAT_POINTER),
            packed.size,
            self._tally_count * self.tally_shards,
            self._round_capacity,
            error,
            len(error),
        )

        if not self._handle:
            self.close()
            raise RuntimeError(error.value.decode())

        self.metadata = metadata
        self._configuration = self._settings()
        self._round = 0
        self.provenance = {
            "backend": self.backend_name,
            "table_layout_version": LAYOUT_VERSION,
            "device": native.monki_gpu_device_name(self._handle).decode(),
            "preparation": preparation,
            "rng": "Philox4x32-10; seed/round/photon/draw-block",
            "seed": self.seed,
            "transport_precision": "float32; fast math disabled",
            "round_precision": "float64",
            "photons_per_dispatch": self.photons_per_dispatch,
            "tally_shards": self.tally_shards,
            "round_capacity": self._round_capacity,
            "max_steps": self.max_steps,
            "bridge_sha256": bridge_hash,
            "kernel_sha256": hashlib.sha256(source).hexdigest(),
            "platform": platform.platform(),
            "python_version": platform.python_version(),
            "gpu_dispatches": 0,
        }

        runtime_info = getattr(native, f"{self.native_prefix}_runtime_info", None)

        if runtime_info is not None:
            runtime_info.argtypes = [ctypes.c_void_p]
            runtime_info.restype = ctypes.c_char_p
            self.provenance["runtime"] = runtime_info(self._handle).decode()

    def run_batch(self, controls: object) -> RawTallies:
        """Execute exactly one complete statistical round."""

        return self._run_group([controls])[0]

    def run_batches(self, controls: list) -> list:
        """Schedule small equal-sized rounds together with independent accumulators."""

        # Validate the complete request before any GPU work, including the tail.
        self._validate_requests(controls)
        results = []
        first = 0
        _, nx, ny, _, _ = self._shape

        while first < len(controls):
            photons = controls[first].photons_per_column * nx * ny
            capacity = max(1, min(self._round_capacity, self.photons_per_dispatch // photons))
            last = first + 1

            while last < min(len(controls), first + capacity):
                if controls[last].photons_per_column != controls[first].photons_per_column:
                    break

                last += 1

            results.extend(self._run_group(controls[first:last]))
            first = last

        return results

    def _validate_requests(self, controls: list) -> None:
        if self._handle is None or self._native is None:
            raise RuntimeError(f"{self.backend_name} must be prepared before running a batch")

        if self._settings() != self._configuration:
            raise ValueError(
                f"{self.backend_name} settings changed; reprepare before running another batch"
            )

        _, nx, ny, _, _ = self._shape

        for offset, request in enumerate(controls, start=1):
            if not isinstance(request, BatchControls):
                raise TypeError("run_batch expects BatchControls")

            if request.round_index != self._round + offset:
                raise ValueError("round_index must be sequential, starting at one")

            if request.round_index >= 2**32:
                raise ValueError("round_index exceeds the 32-bit RNG identity range")

            if request.photons_per_column * nx * ny >= 2**32:
                raise ValueError(
                    f"{self.backend_name} photon identity exceeds its 32-bit round range"
                )

    def _run_group(self, controls: list) -> list[RawTallies]:
        """Reduce GPU dispatches into separate float64 totals for each requested round."""

        self._validate_requests(controls)
        assert self._native is not None
        nmat, nx, ny, nz, angles = self._shape
        photons_per_round = controls[0].photons_per_column * nx * ny
        photon_count = photons_per_round * len(controls)
        total = np.zeros((len(controls), self._tally_count), dtype=np.float64)
        scratch = np.empty((len(controls), self.tally_shards, self._tally_count), dtype=np.float32)
        error = ctypes.create_string_buffer(8192)

        try:
            for first in range(0, photon_count, self.photons_per_dispatch):
                count = min(self.photons_per_dispatch, photon_count - first)
                command = (ctypes.c_uint32 * 8)(
                    count,
                    first,
                    controls[0].round_index,
                    self.seed,
                    self.max_steps,
                    int(self.profile),
                    photons_per_round,
                    len(controls),
                )
                status = self._native.monki_gpu_run(
                    self._handle, command, scratch.ctypes.data_as(FLOAT_POINTER), error, len(error)
                )

                if status or not np.all(np.isfinite(scratch)):
                    raise RuntimeError(
                        error.value.decode() or f"Nonfinite {self.backend_name} tallies"
                    )

                total += scratch.sum(axis=1, dtype=np.float64)
                self.provenance["gpu_dispatches"] += 1

            if self.profile:
                self.provenance["gpu_seconds"] = self._native.monki_gpu_gpu_seconds(self._handle)

            results = []
            stokes_count = nmat * nx * ny * angles

            for request, values in zip(controls, total, strict=True):
                stokes = values[:stokes_count].reshape((nmat, nx, ny, angles), order="F").copy()
                paths = weights = None

                if self.save_path:
                    path_count = nx * ny * nz * angles
                    paths = (
                        values[stokes_count : stokes_count + path_count]
                        .reshape((nx, ny, nz, angles), order="F")
                        .copy()
                    )
                    weights = (
                        values[stokes_count + path_count :]
                        .reshape((nx, ny, angles), order="F")
                        .copy()
                    )

                results.append(RawTallies(request, stokes, paths, weights))

            self._round = controls[-1].round_index
            return results
        except BaseException:
            self.close()
            raise

    def close(self) -> None:
        """Release resident GPU buffers; owned returned arrays stay valid."""

        if self._handle is not None and self._native is not None:
            self._native.monki_gpu_close(self._handle)

        self._handle = None
        self._native = None
        self.metadata = None
        self._configuration = None

    def _settings(self) -> tuple:
        return (
            self.device,
            self.save_path,
            self.seed,
            self.photons_per_dispatch,
            self.max_steps,
            self.profile,
            self.tally_shards,
        )
