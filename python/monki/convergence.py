"""Shared AMF precision stopping, independent of the transport backend."""

from dataclasses import dataclass

import numpy as np

from .execution.tallies import BatchControls, RawTallies


@dataclass(frozen=True)
class AMFConvergence:
    """Stop when every physical layer/pixel/view has relative SE below the target.

    Repeated checking and rare-path tails can make the delta-method precision
    estimate optimistic. ReferenceCase.rounds remains a hard cap.
    """

    relative_error: float = 0.01
    minimum_rounds: int = 128
    check_every: int = 32
    consecutive_checks: int = 3

    def __post_init__(self):
        if not np.isfinite(self.relative_error) or not 0 < self.relative_error < 1:
            raise ValueError("relative_error must be finite and in (0, 1)")

        for name in ("minimum_rounds", "check_every", "consecutive_checks"):
            value = getattr(self, name)

            if type(value) is not int or value < 1:
                raise ValueError(f"{name} must be a positive integer")

        if self.minimum_rounds < 32:
            raise ValueError("minimum_rounds must be at least 32")


class AMFMoments:
    """Online covariance of round path/weight sums; memory does not grow with rounds."""

    def __init__(self):
        self.count = 0
        self.path_mean = self.weight_mean = None
        self.path_m2 = self.weight_m2 = self.covariance = None

    def add(self, batch: RawTallies):
        """Accumulate a complete round, excluding the nonphysical vacuum layer."""

        if batch.weighted_paths is None or batch.weights is None:
            raise ValueError("AMF stopping requires layer path diagnostics")

        paths = batch.weighted_paths[:, :, :-1, :]
        weights = batch.weights[:, :, None, :]

        if not np.all(np.isfinite(paths)) or not np.all(np.isfinite(weights)):
            raise ValueError("Nonfinite AMF tallies")

        if self.count == 0:
            self.path_mean = np.zeros_like(paths)
            self.weight_mean = np.zeros_like(weights)
            self.path_m2 = np.zeros_like(paths)
            self.weight_m2 = np.zeros_like(weights)
            self.covariance = np.zeros_like(paths)

        assert self.path_mean is not None and self.weight_mean is not None
        assert self.path_m2 is not None and self.weight_m2 is not None
        assert self.covariance is not None
        self.count += 1
        path_delta = paths - self.path_mean
        weight_delta = weights - self.weight_mean
        self.path_mean += path_delta / self.count
        self.weight_mean += weight_delta / self.count
        self.path_m2 += path_delta * (paths - self.path_mean)
        self.weight_m2 += weight_delta * (weights - self.weight_mean)
        self.covariance += path_delta * (weights - self.weight_mean)

    def estimate(self, thickness):
        """Return AMF and its SE, retaining numerator/denominator covariance."""

        if self.count < 2:
            raise ValueError("At least two rounds are needed for AMF uncertainty")

        assert self.path_mean is not None and self.weight_mean is not None
        assert self.path_m2 is not None and self.weight_m2 is not None
        assert self.covariance is not None

        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = self.path_mean / self.weight_mean
            # This equals sum((P_round - ratio*W_round)**2), but uses stable
            # centred moments to bound memory use across rounds.
            residual_m2 = self.path_m2 - 2 * ratio * self.covariance + ratio**2 * self.weight_m2
            scale = self.path_m2 + 2 * np.abs(ratio * self.covariance) + ratio**2 * self.weight_m2

            if np.any(residual_m2 < -128 * np.finfo(float).eps * np.maximum(scale, 1e-300)):
                raise ArithmeticError("Negative AMF variance exceeds roundoff")

            error = np.sqrt(np.maximum(residual_m2, 0) / (self.count - 1) / self.count)
            error /= self.weight_mean
            shape = thickness[None, None, :, None]
            return ratio / shape, error / shape


class AMFStopping:
    """Record precision checks and stop only after sustained estimated precision."""

    def __init__(self, policy: AMFConvergence, thickness):
        self.policy = policy
        self.thickness = thickness
        self.moments = AMFMoments()
        self.history = []
        self.streak = 0
        self.converged = False

    def batches(self, backend, case):
        """Schedule equally sized rounds in bounded groups on either backend."""

        for start in range(1, case.rounds + 1, self.policy.check_every):
            end = min(start + self.policy.check_every - 1, case.rounds)
            requests = [BatchControls(i, case.photons_per_column) for i in range(start, end + 1)]
            batches = backend.run_batches(requests)

            if len(batches) != len(requests):
                raise ValueError("Backend returned the wrong number of statistical rounds")

            for request, batch in zip(requests, batches, strict=True):
                if not isinstance(batch, RawTallies) or batch.controls != request:
                    raise ValueError("Backend returned invalid round tallies")

                self.moments.add(batch)
                yield batch

            if end < self.policy.minimum_rounds:
                continue

            amf, error = self.moments.estimate(self.thickness)
            relative = np.divide(
                error, np.abs(amf), out=np.full_like(error, np.inf), where=amf != 0
            )
            maximum = float(np.max(relative))
            passed = bool(np.all(np.isfinite(relative)) and maximum < self.policy.relative_error)
            self.streak = self.streak + 1 if passed else 0
            self.converged = self.streak >= self.policy.consecutive_checks
            self.history.append(
                {
                    "rounds": end,
                    "max_relative_se": maximum if np.isfinite(maximum) else None,
                    "below_target": passed,
                    "streak": self.streak,
                }
            )

            if self.converged:
                return
