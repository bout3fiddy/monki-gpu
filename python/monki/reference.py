"""Isolated execution of the supplied reference solver."""

import hashlib
import json
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path

from .environment import native_environment
from .input import ReferenceCase
from .output import SimulationResult


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture_text(source: Path, case: ReferenceCase) -> tuple[str, tuple[int, int, int]]:
    """Patch only known positional testrun fixtures, checking their original layout."""

    if source.parent.name != "testrun":
        raise ValueError("Only the supplied testrun input format is supported")

    # The original parser reads fixed line positions. Check the known layout before
    # changing counts so an older input format cannot silently shift the physics.
    lines = source.read_text().splitlines()
    checks = {
        4: "number of columns in horizontal x",
        5: "number of columns in horizontal y",
        6: "number of layers",
        16: "number of photons per column per round",
        17: "maximum number of rounds",
        18: "stop simulation when result has converged",
        19: "maximum allowed reflectance standard deviation",
        24: "size of Stokes vector",
    }

    for index, marker in checks.items():
        if index >= len(lines) or marker not in lines[index]:
            raise ValueError(f"Unsupported positional input at line {index + 1}: {source}")

    dimensions = (
        int(lines[4].split()[0]),
        int(lines[5].split()[0]),
        int(lines[6].split()[0]),
    )

    for index, value in ((16, case.photons_per_column), (17, case.rounds), (18, ".false.")):
        lines[index] = f"{str(value):<24}:" + lines[index].split(":", 1)[1]

    return "\n".join(lines) + "\n", dimensions


@dataclass(frozen=True)
class ReferenceRunner:
    """Run the legacy solver in its own directory and retain a reproducible receipt."""

    executable: Path
    upstream: Path
    threads: int = 1
    timeout_seconds: float = 120

    def run(self, case: ReferenceCase, directory: str | Path) -> SimulationResult:
        """Execute a fixed-round case and return owned arrays from a complete result.

        The caller chooses a new run directory. Inputs, logs, and hashes stay there
        after failure so a native error can be diagnosed without rerunning it.
        """

        executable = self.executable.resolve()
        upstream = self.upstream.resolve()
        directory = Path(directory).resolve()

        if type(self.threads) is not int or self.threads < 1:
            raise ValueError("threads must be a positive integer")

        if not executable.is_file():
            raise FileNotFoundError(f"Build the reference executable first: {executable}")

        source = upstream / "INPUT" / case.input_file
        text, (nx, ny, nz) = fixture_text(source, case)

        if nx > 30 or ny > 20 or nz > 35 or len(case.view_zeniths) > 3:
            raise ValueError("case exceeds reference harness build limits: 30x20x35, 3 angles")

        photons_per_round = case.photons_per_column * nx * ny

        # Requiring per-column divisibility also preserves even injection per thread.
        if case.photons_per_column % self.threads:
            raise ValueError("photons_per_column must be divisible by threads in legacy reference")

        if photons_per_round // self.threads >= 2**31:
            raise ValueError("legacy photon loop exceeds its 32-bit index")

        if len(str(directory)) > 240:
            raise ValueError("run directory is too long for legacy Fortran path buffers")

        # MONKI resolves both inputs and outputs beneath MONKI_HOME. Give each run
        # a fresh home to protect archived results and avoid accepting a stale file.
        directory.mkdir(parents=True, exist_ok=False)
        (directory / "INPUT").mkdir()
        (directory / "OUTPUT").mkdir()
        (directory / "DATA").symlink_to(upstream / "DATA", target_is_directory=True)
        input_path = directory / "INPUT/case.in"
        input_path.write_text(text)
        # Argument order is the legacy CLI's public interface, including the angle
        # count before its values. No shell interpolation is involved.
        command = [
            str(executable),
            str(case.surface_albedo),
            str(case.solar_zenith),
            str(case.solar_azimuth),
            str(case.view_azimuth),
            str(len(case.view_zeniths)),
            *map(str, case.view_zeniths),
            str(int(case.backward)),
            "case.in",
            "result",
        ]
        env = native_environment()
        env.update(
            MONKI_HOME=str(directory) + "/",
            OMP_NUM_THREADS=str(self.threads),
            OMP_DYNAMIC="FALSE",
            GOMP_STACKSIZE="128M",
            OMP_STACKSIZE="128M",
        )
        # Record actual inputs and executable identity for reproducibility.
        # Compiler and thread settings matter because upstream RNG streams depend on them.
        manifest = {
            "schema_version": 1,
            "case": asdict(case),
            "command": command,
            "threads": self.threads,
            "photons_per_round": photons_per_round,
            "requested_total_photons": photons_per_round * case.rounds,
            "source_input_sha256": file_hash(source),
            "generated_input_sha256": file_hash(input_path),
            "binary_sha256": file_hash(executable),
            "rng": "upstream seed=23274; compiler RNG; thread-dependent streams",
            "status": "running",
        }
        vendor_manifest = upstream.parent / "manifest.json"

        if vendor_manifest.is_file():
            manifest["vendor_manifest_sha256"] = file_hash(vendor_manifest)

        build_receipt = executable.parent / "build.json"

        if build_receipt.is_file():
            manifest["build"] = json.loads(build_receipt.read_text())

        receipt = directory / "run.json"
        receipt.write_text(json.dumps(manifest, indent=2) + "\n")
        started = time.monotonic()

        try:
            with (directory / "solver.log").open("w") as log:
                completed = subprocess.run(
                    command,
                    cwd=directory,
                    env=env,
                    stdin=subprocess.DEVNULL,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    timeout=self.timeout_seconds,
                    check=True,
                )

            # Fortran STOP can return success on a failed run. Require a new result
            # and matching completed counts to verify completion.
            manifest["returncode"] = completed.returncode
            outputs = list((directory / "OUTPUT").glob("*.nc"))

            if len(outputs) != 1:
                raise RuntimeError(
                    "Reference did not produce exactly one NetCDF file; see solver.log"
                )

            result = SimulationResult.load(outputs[0])

            for name, expected in (
                ("n_rounds", case.rounds),
                ("n_photons_per_column", case.photons_per_column),
            ):
                if name not in result.values or int(result.values[name].item()) != expected:
                    raise RuntimeError(f"Incomplete reference output: {name}")

            manifest.update(
                status="complete",
                output=str(outputs[0].relative_to(directory)),
                output_sha256=file_hash(outputs[0]),
            )
            return result
        except Exception as exc:
            manifest.update(status="failed", error=str(exc))
            raise
        finally:
            # Persist failures and elapsed time too; an interrupted or failed run must
            # never be indistinguishable from a scientifically accepted result.
            manifest["elapsed_seconds"] = time.monotonic() - started
            receipt.write_text(json.dumps(manifest, indent=2) + "\n")
