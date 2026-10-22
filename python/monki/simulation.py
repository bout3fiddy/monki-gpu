"""Shared round statistics following Victor's main.f90 and accumulate routine."""

import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from .convergence import AMFConvergence, AMFStopping
from .execution.backend import Backend
from .execution.tallies import BatchControls, RawTallies
from .input import ReferenceCase
from .output import SimulationResult
from .reference import file_hash


def simulate(
    case: ReferenceCase,
    backend: Backend,
    output: str | Path,
    *,
    convergence: AMFConvergence | None = None,
) -> SimulationResult:
    """Run fixed, equally sized rounds and save a complete NetCDF product.

    With AMFConvergence, case.rounds sets the maximum round budget.
    The saved receipt distinguishes convergence from exhausting that budget.

    This helper owns the prepared session and closes it on success or failure. Each
    batch must contain one complete statistical round to preserve MONKI's
    uncertainty estimator. The same path estimator and
    stopping rule apply to every backend.
    """

    if convergence is not None and case.rounds < convergence.minimum_rounds:
        raise ValueError("round cap must be at least minimum_rounds")

    output = Path(output).resolve()

    receipt = output.with_suffix(".run.json")

    if output.exists() or receipt.exists():
        raise FileExistsError(f"Refusing to overwrite an existing result: {output}")

    with backend:
        backend.prepare(case)
        metadata = backend.metadata

        if metadata is None:
            raise RuntimeError("Backend did not supply prepared scene metadata")

        total = rounds_sum = rounds_squared = paths = weights = None
        history = []

        stopping = None

        if convergence is not None:
            if "AMF_layer" not in metadata.values:
                raise ValueError("AMF stopping requires layer path diagnostics")

            thickness = np.diff(metadata.values["alt"]) / 1000
            stopping = AMFStopping(convergence, thickness)
            batches = stopping.batches(backend, case)
        else:
            requests = [
                BatchControls(index, case.photons_per_column) for index in range(1, case.rounds + 1)
            ]
            batches = backend.run_batches(requests)

            if len(batches) != case.rounds:
                raise ValueError("Backend returned the wrong number of statistical rounds")

        for index, batch in enumerate(batches, start=1):
            if not isinstance(batch, RawTallies):
                raise TypeError("Backend must return RawTallies")

            if batch.controls != BatchControls(index, case.photons_per_column):
                raise ValueError("Backend returned tallies for a different round")

            if total is None:
                total = np.zeros_like(batch.stokes)
                rounds_sum = np.zeros_like(batch.stokes)
                rounds_squared = np.zeros_like(batch.stokes)

                if batch.weighted_paths is not None:
                    paths = np.zeros_like(batch.weighted_paths)
                    weights = np.zeros_like(batch.weights)

            assert rounds_sum is not None and rounds_squared is not None
            total += batch.stokes
            normalized = batch.stokes / case.photons_per_column
            rounds_sum += normalized
            rounds_squared += normalized**2
            # Match the original successive spatial sums and divisions, retaining
            # all Stokes components in the convergence history.
            _, nx, ny, angles = total.shape
            history.append(
                rounds_sum.sum(axis=1).sum(axis=1).sum(axis=1) / nx / ny / angles / index
            )

            if paths is not None and weights is not None:
                if batch.weighted_paths is None or batch.weights is None:
                    raise ValueError("Backend dropped path diagnostics during a session")

                paths += batch.weighted_paths
                weights += batch.weights

        completed_rounds = len(history)
        assert total is not None and rounds_sum is not None and rounds_squared is not None
        radiance = total / completed_rounds / case.photons_per_column
        uncertainty = np.full_like(radiance, 999.0)

        # Upstream uses 999 as 'uncertainty unavailable' through round 30. This
        # is a standard error of the mean across rounds.
        if completed_rounds > 30:
            variance_sum = rounds_squared - (1.0 / completed_rounds) * rounds_sum**2
            tolerance = 64 * np.finfo(float).eps * np.maximum(rounds_squared, 1e-300)

            if np.any(variance_sum < -tolerance):
                raise ArithmeticError("Negative round variance exceeds floating-point roundoff")

            uncertainty = np.sqrt(np.maximum(variance_sum, 0) / (completed_rounds - 1)) / np.sqrt(
                completed_rounds
            )

        radiance_history = np.asarray(history)

        # The backward transport convention reverses U only after accumulation.
        if case.backward and radiance.shape[0] >= 3:
            radiance[2] *= -1
            radiance_history[:, 2] *= -1

        intensity = radiance[0]
        polarization = np.sqrt(np.sum(radiance[1:] ** 2, axis=0))
        dop = np.divide(
            polarization, intensity, out=np.zeros_like(intensity), where=intensity > 1e-16
        )
        values = {name: value.copy() for name, value in metadata.values.items()}
        values["radiance"] = radiance.transpose(3, 2, 1, 0)[..., None]
        values["std_radiance"] = uncertainty.transpose(3, 2, 1, 0)[..., None]
        values["DOP"] = dop.transpose(2, 1, 0)[..., None]
        values["radiance_convergence"] = radiance_history[..., None]
        values["n_rounds"] = np.full_like(values["n_rounds"], completed_rounds)
        values["n_photons_per_column"] = np.full_like(
            values["n_photons_per_column"], case.photons_per_column
        )

        if paths is not None and weights is not None:
            if np.any(weights <= 0):
                raise ValueError("Layer AMF is undefined for a pixel with no detected weight")

            # Transport lengths are kilometres; altitude metadata is in metres.
            # The final vacuum layer is deliberately omitted from the public output.
            thickness = np.diff(values["alt"]) / 1000.0
            amf = paths[:, :, :-1, :] / weights[:, :, None, :] / thickness[None, None, :, None]
            values["AMF_layer"] = amf.transpose(3, 2, 1, 0)[..., None]

        dimensions = metadata.dimensions.copy()
        attributes = metadata.attributes.copy()

        if stopping is not None:
            _, amf_error = stopping.moments.estimate(stopping.thickness)
            values["std_AMF_layer"] = amf_error.transpose(3, 2, 1, 0)[..., None]
            dimensions["std_AMF_layer"] = dimensions["AMF_layer"]
            attributes["std_AMF_layer"] = {"long_name": "Layer AMF standard error", "units": "-"}

        result = SimulationResult(output, values, dimensions, attributes)
        result.save()
        provenance = {
            **backend.provenance,
            "status": "complete",
            "completed_rounds": completed_rounds,
            "total_photons": completed_rounds * case.photons_per_column * nx * ny,
            "output_sha256": file_hash(output),
        }

        if stopping is not None:
            provenance["status"] = "converged" if stopping.converged else "max_rounds_reached"
            provenance["amf_convergence"] = {
                "policy": asdict(stopping.policy),
                "converged": stopping.converged,
                "checks": stopping.history,
                "uncertainty": "pointwise relative standard error",
            }

        with receipt.open("x") as stream:
            stream.write(json.dumps(provenance, indent=2) + "\n")

        return result
