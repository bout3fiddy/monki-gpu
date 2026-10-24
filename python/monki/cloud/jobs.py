"""Versioned simulation requests and bounded, checksummed inline artifacts."""

import base64
import hashlib
import json
from dataclasses import asdict
from pathlib import Path
from tempfile import TemporaryDirectory

from ..convergence import AMFConvergence
from ..execution import CUDA
from ..input import ReferenceCase
from ..simulation import simulate

FIXTURES = {
    "testrun/MONKI_MLS_cloudfree_340_340_test.in",
    "testrun/MONKI_MLS_1layercloud_340_340_test.in",
    "testrun/MONKI_MLS_2boxclouds_340_340_test.in",
}
MAX_ARTIFACT_BYTES = 6 * 1024 * 1024


def simulation_request(
    case: ReferenceCase,
    *,
    save_path: bool = True,
    seed: int = 23274,
    convergence: AMFConvergence | None = None,
) -> dict:
    """Describe a full CUDA simulation using fixtures bundled in the worker image."""

    request = {
        "version": 1,
        "case": asdict(case),
        "save_path": save_path,
        "seed": seed,
        "convergence": asdict(convergence) if convergence else None,
    }
    validate_request(request)
    return request


def validate_request(request: dict) -> tuple[ReferenceCase, AMFConvergence | None]:
    """Reject unsupported fields and inputs before allocating a GPU session."""

    if (
        not isinstance(request, dict)
        or set(request) != {"version", "case", "save_path", "seed", "convergence"}
        or type(request["version"]) is not int
        or request["version"] != 1
    ):
        raise ValueError("Unsupported MONKI job schema")

    case = ReferenceCase(**request["case"])

    if case.input_file not in FIXTURES or case.backward or len(case.view_zeniths) > 3:
        raise ValueError(
            "Remote CUDA jobs require a supported forward fixture and up to three views"
        )

    if type(request["save_path"]) is not bool:
        raise ValueError("save_path must be a bool")

    if type(request["seed"]) is not int or not 0 <= request["seed"] < 2**32:
        raise ValueError("seed must be an unsigned 32-bit integer")

    box = "2boxclouds" in case.input_file

    if request["save_path"] == box:
        raise ValueError("Box clouds require save_path=False; clear/layered fixtures require True")

    policy = (
        AMFConvergence(**request["convergence"]) if request["convergence"] is not None else None
    )

    if policy is not None and (not request["save_path"] or case.rounds < policy.minimum_rounds):
        raise ValueError("Convergence requires paths and a sufficient round budget")

    # Fixed execution retains round arrays; bound memory and result size for this
    # initial fixture-only serverless protocol. Larger workloads need streaming/storage.
    if case.rounds > 8192 or (box and case.rounds > 256):
        raise ValueError("Remote round budget exceeds validated inline-job limits")

    return case, policy


def execute(request: dict) -> dict:
    """Run one isolated job and return NetCDF plus its unchanged scientific receipt."""

    case, policy = validate_request(request)

    with TemporaryDirectory(prefix="monki-job-") as temporary:
        output = Path(temporary) / "result.nc"
        simulate(
            case,
            CUDA(save_path=request["save_path"], seed=request["seed"]),
            output,
            convergence=policy,
        )
        artifacts = {}
        total = 0

        for path in (output, output.with_suffix(".run.json")):
            data = path.read_bytes()
            total += len(data)

            if total > MAX_ARTIFACT_BYTES:
                raise ValueError("Result exceeds inline artifact limit")

            artifacts[path.name] = {
                "sha256": hashlib.sha256(data).hexdigest(),
                "base64": base64.b64encode(data).decode("ascii"),
            }

        return {"version": 1, "artifacts": artifacts}


def save_artifacts(result: dict, destination: Path) -> None:
    """Verify the complete response before creating a fresh local result directory."""

    if result.get("version") != 1 or set(result.get("artifacts", {})) != {
        "result.nc",
        "result.run.json",
    }:
        raise ValueError("Invalid MONKI result schema")

    decoded = {}
    total = 0

    for name, artifact in result["artifacts"].items():
        encoded = artifact["base64"]

        if not isinstance(encoded, str) or len(encoded) > MAX_ARTIFACT_BYTES * 4 // 3 + 4:
            raise ValueError("Encoded artifact exceeds size limit")

        data = base64.b64decode(encoded, validate=True)
        total += len(data)

        if total > MAX_ARTIFACT_BYTES or hashlib.sha256(data).hexdigest() != artifact["sha256"]:
            raise ValueError("Artifact size/hash verification failed")

        decoded[name] = data

    json.loads(decoded["result.run.json"])
    destination.mkdir(parents=True, exist_ok=False)

    for name, data in decoded.items():
        (destination / name).write_bytes(data)
