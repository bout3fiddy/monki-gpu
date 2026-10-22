"""Round-based uncertainty for CPU/Metal comparisons of the layer AMF ratio."""

from statistics import NormalDist

import numpy as np
from monki.execution import Backend, RawTallies
from monki.output import SimulationResult


class PathRecorder(Backend):
    """Retain per-round path sums for validation."""

    def __init__(self, backend: Backend):
        self.backend = backend
        self.rounds: list[RawTallies] = []

    def prepare(self, scene: object) -> None:
        self.rounds.clear()
        self.backend.prepare(scene)
        self.metadata = self.backend.metadata
        self.provenance = self.backend.provenance

    def run_batch(self, controls: object) -> RawTallies:
        batch = self.backend.run_batch(controls)

        if not isinstance(batch, RawTallies):
            raise TypeError("Validation requires raw round tallies")

        # Box-cloud radiance-only runs need no retained per-round path buffers.
        if batch.weighted_paths is not None:
            self.rounds.append(batch)

        return batch

    def run_batches(self, controls: list) -> list:
        batches = self.backend.run_batches(controls)

        for batch in batches:
            if not isinstance(batch, RawTallies):
                raise TypeError("Validation requires raw round tallies")

            if batch.weighted_paths is not None:
                self.rounds.append(batch)

        return batches

    def close(self) -> None:
        self.backend.close()

    def standard_error(self, result: SimulationResult) -> np.ndarray:
        """Delta-method standard error of sum(weighted paths)/sum(weights)/thickness.

        The numerator and denominator come from the same photons and are correlated.
        Residuals P_round - ratio*W_round retain that covariance. Rare paths can make
        this finite-round uncertainty estimate optimistic.
        """

        if len(self.rounds) <= 30:
            raise ValueError("AMF uncertainty requires more than 30 complete rounds")

        paths = np.stack([batch.weighted_paths for batch in self.rounds])[:, :, :, :-1, :]
        weights = np.stack([batch.weights for batch in self.rounds])[:, :, :, None, :]
        mean_weight = weights.mean(axis=0)

        if np.any(mean_weight <= 0):
            raise ValueError("AMF uncertainty is undefined without detected weight")

        ratio = paths.mean(axis=0) / mean_weight
        residual = paths - ratio * weights
        error = residual.std(axis=0, ddof=1) / np.sqrt(len(self.rounds)) / mean_weight
        thickness = np.diff(result.values["alt"]) / 1000
        error /= thickness[None, None, :, None]
        return error.transpose(3, 2, 1, 0)[..., None]


def compare_paths(
    reference: SimulationResult,
    candidate: SimulationResult,
    reference_rounds: PathRecorder,
    candidate_rounds: PathRecorder,
) -> dict:
    """Screen layer AMFs with independent-run errors and a fixed family error rate."""

    if "AMF_layer" not in reference.values:
        return {"applicable": False, "passed": "AMF_layer" not in candidate.values}

    reference_error = reference_rounds.standard_error(reference)
    candidate_error = candidate_rounds.standard_error(candidate)
    delta = np.abs(candidate.values["AMF_layer"] - reference.values["AMF_layer"])
    family_error_rate = 1e-6
    z = NormalDist().inv_cdf(1 - family_error_rate / (2 * delta.size))
    bound = 1e-10 + z * np.hypot(reference_error, candidate_error)
    finite = np.all(np.isfinite(delta)) and np.all(np.isfinite(bound))
    failed = int(np.count_nonzero(delta > bound))
    return {
        "applicable": True,
        "passed": bool(finite and failed == 0),
        "estimator": "correlated ratio delta-method across rounds",
        "family_error_rate": family_error_rate,
        "absolute_tolerance": 1e-10,
        "compared_values": delta.size,
        "failed_values": failed,
        "max_absolute_error": float(np.max(delta)),
        "max_normalized_error": float(np.max(delta / bound)),
    }
