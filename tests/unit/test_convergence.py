"""Check AMF uncertainty against explicit correlated round residuals."""

import numpy as np
import pytest
from monki import AMFConvergence
from monki.convergence import AMFMoments, AMFStopping
from monki.execution.tallies import BatchControls, RawTallies


def test_online_amf_error_matches_explicit_correlated_ratio():
    rng = np.random.default_rng(123)
    weights = rng.uniform(1, 4, (100, 1, 1, 2))
    paths = weights[:, :, :, None, :] * rng.uniform(1, 5, (100, 1, 1, 4, 2))
    moments = AMFMoments()

    for index in range(100):
        moments.add(
            RawTallies(
                BatchControls(index + 1, 100), np.zeros((3, 1, 1, 2)), paths[index], weights[index]
            )
        )

    thickness = np.array([0.5, 2.0, 5.0])
    actual, error = moments.estimate(thickness)
    ratio = paths[:, :, :, :-1, :].mean(axis=0) / weights.mean(axis=0)[:, :, None, :]
    residual = paths[:, :, :, :-1, :] - ratio * weights[:, :, :, None, :]
    expected_error = residual.std(axis=0, ddof=1) / 10 / weights.mean(axis=0)[:, :, None, :]
    np.testing.assert_allclose(actual, ratio / thickness[None, None, :, None])
    np.testing.assert_allclose(error, expected_error / thickness[None, None, :, None], rtol=1e-12)


@pytest.mark.parametrize(
    "settings",
    [
        {"relative_error": 0},
        {"relative_error": np.nan},
        {"minimum_rounds": 31},
        {"check_every": 0},
        {"consecutive_checks": True},
    ],
)
def test_invalid_stopping_policy(settings):
    with pytest.raises(ValueError):
        AMFConvergence(**settings)


def test_late_outlier_resets_precision_streak():
    class Synthetic:
        def run_batches(self, requests):
            return [
                RawTallies(
                    request,
                    np.ones((3, 1, 1, 1)),
                    np.full((1, 1, 2, 1), 100.0 if request.round_index == 48 else 1.0),
                    np.ones((1, 1, 1)),
                )
                for request in requests
            ]

    class Case:
        rounds = 64
        photons_per_column = 100

    stopping = AMFStopping(AMFConvergence(minimum_rounds=32, check_every=16), np.array([1.0]))
    batches = list(stopping.batches(Synthetic(), Case()))
    assert len(batches) == 64
    assert not stopping.converged
    assert stopping.history[0]["streak"] == 1
    assert stopping.history[1]["streak"] == 0


def test_zero_amf_cannot_satisfy_relative_precision():
    class Synthetic:
        def run_batches(self, requests):
            return [
                RawTallies(
                    request, np.ones((3, 1, 1, 1)), np.zeros((1, 1, 2, 1)), np.ones((1, 1, 1))
                )
                for request in requests
            ]

    class Case:
        rounds = 32
        photons_per_column = 100

    stopping = AMFStopping(AMFConvergence(minimum_rounds=32, consecutive_checks=1), np.array([1.0]))
    list(stopping.batches(Synthetic(), Case()))
    assert not stopping.converged
    assert stopping.history[-1]["max_relative_se"] is None
