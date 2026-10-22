"""Persistent local Fortran sessions retaining Victor's transport and RNG state."""

import json
import os
import subprocess
from contextlib import suppress
from dataclasses import asdict, replace
from pathlib import Path
from queue import Empty, Queue
from tempfile import TemporaryDirectory
from threading import Thread

import numpy as np

from ..input import ReferenceCase
from ..output import SimulationResult
from ..reference import file_hash, fixture_text
from .backend import Backend
from .tallies import BatchControls, RawTallies

ROOT = Path(__file__).resolve().parents[3]


class CPU(Backend):
    """Run a local persistent worker; one caller may operate a session at a time.

    Build with scripts/build_cpu.py first. Native failures are isolated in the worker
    process, and protocol waits are bounded. Use a context manager to release resources.
    """

    def __init__(
        self,
        threads: int = 1,
        *,
        save_path: bool = True,
        seed: int = 23274,
        timeout_seconds: float = 120,
        executable: Path | None = None,
        upstream: Path = ROOT / "fortran/monki/upstream",
    ):
        if type(threads) is not int or threads < 1:
            raise ValueError("threads must be a positive integer")

        if type(save_path) is not bool:
            raise ValueError("save_path must be a bool")

        if not np.isfinite(timeout_seconds) or timeout_seconds <= 0:
            raise ValueError("timeout_seconds must be positive and finite")

        # Leave headroom for the Fortran per-thread seed-array offsets.
        if type(seed) is not int or not 0 <= seed <= 2**30:
            raise ValueError("seed must be an integer in [0, 2**30]")

        self.seed = seed
        profile = "paths" if save_path else "radiance"
        self.threads = threads
        self.save_path = save_path
        self.timeout_seconds = timeout_seconds
        self.executable = executable or ROOT / f"build/cpu/{profile}/monki_cpu"
        self.upstream = upstream
        self.provenance: dict = {}
        self.metadata: SimulationResult | None = None
        self._process: subprocess.Popen | None = None
        self._temporary: TemporaryDirectory | None = None
        self._reader: Thread | None = None
        self._messages: Queue[str | None] = Queue()
        self._round = 0
        self._shape = (0, 0, 0, 0, 0)
        self.log: list[str] = []

    def prepare(self, scene: object) -> None:
        """Load a validated monochromatic fixture once and retain its prepared tables."""

        self.close()

        if not isinstance(scene, ReferenceCase):
            raise TypeError("CPU.prepare expects a ReferenceCase")

        if not self.executable.is_file():
            raise FileNotFoundError("Build the CPU worker with uv run python scripts/build_cpu.py")

        source = self.upstream / "INPUT" / scene.input_file
        # rdsize pauses interactively below 30 configured rounds. The worker does
        # not use that loop count; reserve a safe metadata history and let the
        # driver choose the actual number of completed statistical rounds.
        text, (nx, ny, nz) = fixture_text(source, replace(scene, rounds=32))

        if nx > 30 or ny > 20 or nz + 1 > 35 or len(scene.view_zeniths) > 3:
            raise ValueError("case exceeds CPU build limits: 30x20x35 including vacuum, 3 angles")

        lines = text.splitlines()

        if "save the layer AMFs" not in lines[38] or "save the single" not in lines[39]:
            raise ValueError("Unsupported positional diagnostic settings")

        requested_paths = lines[38].split()[0].lower() == ".true."

        if requested_paths != self.save_path:
            raise ValueError("CPU save_path must match the fixture's layer AMF setting")

        if lines[39].split()[0].lower() != ".false.":
            raise ValueError("CPU scattering contribution diagnostics are not implemented")

        self._validate_count(scene.photons_per_column, nx, ny)
        self._temporary = TemporaryDirectory(prefix="monki-cpu-")
        directory = Path(self._temporary.name)
        self.log = []
        self._messages = Queue()

        try:
            (directory / "INPUT").mkdir()
            (directory / "OUTPUT").mkdir()
            (directory / "DATA").symlink_to((self.upstream / "DATA").resolve())
            (directory / "INPUT/case.in").write_text(text)
            command = [
                str(self.executable.resolve()),
                str(scene.surface_albedo),
                str(scene.solar_zenith),
                str(scene.solar_azimuth),
                str(scene.view_azimuth),
                str(len(scene.view_zeniths)),
                *map(str, scene.view_zeniths),
                str(int(scene.backward)),
                "case.in",
                "result",
            ]
            env = os.environ.copy()
            env.update(
                MONKI_HOME=str(directory) + "/",
                OMP_NUM_THREADS=str(self.threads),
                OMP_DYNAMIC="FALSE",
                MONKI_CPU_SEED=str(self.seed),
                GOMP_STACKSIZE="128M",
                OMP_STACKSIZE="128M",
            )
            self._process = subprocess.Popen(
                command,
                cwd=directory,
                env=env,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
            # Drain stdout continuously so verbose legacy output cannot fill a pipe
            # and deadlock transport. A queue also lets the caller enforce timeouts.
            self._reader = Thread(target=self._read_output, args=(self._process,), daemon=True)
            self._reader.start()
            ready = self._wait("MONKI_CPU_READY")
            version, x, y, z, nmat, angles, paths, threads = map(int, ready.split()[1:])

            if (version, x, y, z, angles, paths, threads) != (
                1,
                nx,
                ny,
                nz + 1,
                len(scene.view_zeniths),
                int(self.save_path),
                self.threads,
            ):
                raise RuntimeError(f"CPU worker configuration mismatch: {ready}")

            self._shape = (nmat, x, y, z, angles)
            self.metadata = SimulationResult.load(directory / "OUTPUT/scene.nc")
            self._round = 0
            self.provenance = {
                "backend": "CPU",
                "case": asdict(scene),
                "threads": self.threads,
                "save_path": self.save_path,
                "binary_sha256": file_hash(self.executable),
                "source_input_sha256": file_hash(source),
                "generated_input_sha256": file_hash(directory / "INPUT/case.in"),
                "rng": f"upstream seed={self.seed}; compiler RNG; thread-dependent streams",
                "seed": self.seed,
            }
            build_receipt = self.executable.parent / "build.json"

            if build_receipt.is_file():
                self.provenance["build"] = json.loads(build_receipt.read_text())
        except BaseException:
            self.close()
            raise

    def transport_tables(self) -> np.ndarray:
        """Copy the prepared float64 table export without running any CPU photons."""

        if self.metadata is None or self._temporary is None:
            raise RuntimeError("CPU must be prepared before exporting transport tables")

        return np.fromfile(Path(self._temporary.name) / "transport.bin", dtype=np.float64)

    def run_batch(self, controls: object) -> RawTallies:
        """Advance the retained random streams by one round and copy its raw sums."""

        if self.metadata is None or self._process is None or self._temporary is None:
            raise RuntimeError("CPU must be prepared before running a batch")

        if not isinstance(controls, BatchControls):
            raise TypeError("run_batch expects BatchControls")

        if controls.round_index != self._round + 1:
            raise ValueError("round_index must be sequential, starting at one")

        nmat, nx, ny, nz, angles = self._shape
        self._validate_count(controls.photons_per_column, nx, ny)
        stream = self._process.stdin
        assert stream is not None

        try:
            stream.write(f"{controls.round_index} {controls.photons_per_column}\n")
            stream.flush()
            reply = self._wait("MONKI_CPU_BATCH")

            if list(map(int, reply.split()[1:])) != [
                controls.round_index,
                controls.photons_per_column,
            ]:
                raise RuntimeError(f"CPU batch acknowledgement mismatch: {reply}")

            data = np.fromfile(Path(self._temporary.name) / "batch.bin", dtype=np.float64)
            stokes_count = nmat * nx * ny * angles
            path_count = nx * ny * nz * angles if self.save_path else 0
            weight_count = nx * ny * angles if self.save_path else 0

            if data.size != stokes_count + path_count + weight_count or not np.all(
                np.isfinite(data)
            ):
                raise RuntimeError("CPU worker returned malformed/nonfinite tallies")

            stokes = data[:stokes_count].reshape((nmat, nx, ny, angles), order="F").copy()
            paths = weights = None

            if self.save_path:
                paths = (
                    data[stokes_count : stokes_count + path_count]
                    .reshape((nx, ny, nz, angles), order="F")
                    .copy()
                )
                weights = data[-weight_count:].reshape((nx, ny, angles), order="F").copy()

            self._round = controls.round_index
            return RawTallies(controls, stokes, paths, weights)
        except BaseException:
            # A failed protocol has unknown RNG progress. Close the session to prevent
            # skipped photons or stale sums in subsequent calls.
            self.close()
            raise

    def close(self) -> None:
        """Release native memory and temporary files, including after a worker failure."""

        process = self._process
        self._process = None
        self.metadata = None

        if process is not None:
            if process.poll() is None:
                process.terminate()

            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()

            if self._reader is not None:
                self._reader.join(timeout=5)

            for stream in (process.stdin, process.stdout):
                if stream is not None:
                    # A broken pipe may surface again while buffered stdin closes;
                    # it must not prevent freeing the session directory.
                    with suppress(OSError):
                        stream.close()

        if self._temporary is not None:
            self._temporary.cleanup()
            self._temporary = None

    def _validate_count(self, count: int, nx: int, ny: int) -> None:
        if count % self.threads:
            raise ValueError("photons_per_column must be divisible by threads")

        if count * nx * ny // self.threads >= 2**31 - 1:
            raise ValueError("legacy photon loop exceeds its 32-bit index")

    def _read_output(self, process: subprocess.Popen) -> None:
        assert process is not None and process.stdout is not None

        for line in process.stdout:
            self.log.append(line.rstrip())

            if line.lstrip().startswith("MONKI_CPU_"):
                self._messages.put(line.strip())

        self._messages.put(None)

    def _wait(self, marker: str) -> str:
        try:
            message = self._messages.get(timeout=self.timeout_seconds)
        except Empty as exc:
            raise TimeoutError(f"CPU worker timed out waiting for {marker}") from exc

        if message is None or not message.startswith(marker + " "):
            raise RuntimeError(
                f"CPU worker failed waiting for {marker}:\n" + "\n".join(self.log[-20:])
            )

        return message
