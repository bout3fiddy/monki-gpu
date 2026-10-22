"""Execute the CPU session against the independently built original program."""

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from monki import ReferenceCase, ReferenceRunner, SimulationResult, simulate
from monki.execution import CPU
from monki.execution.tallies import BatchControls
from monki.validation import check_result, compare

ROOT = Path(__file__).resolve().parents[2]
UPSTREAM = ROOT / "fortran/monki/upstream"
CASE = ReferenceCase("clear", "testrun/MONKI_MLS_cloudfree_340_340_test.in")
pytestmark = pytest.mark.integration


@pytest.mark.parametrize(
    "threads,backward,paths",
    [(1, False, True), (2, False, True), (1, True, True), (2, True, True)],
)
def test_cpu_matches_original_round_statistics(tmp_path, threads, backward, paths):
    case = replace(CASE, backward=backward)
    profile = "paths" if paths else "radiance"
    reference = ReferenceRunner(ROOT / f"build/reference/{profile}/monki", UPSTREAM, threads).run(
        case, tmp_path / "reference"
    )
    backend = CPU(threads=threads, save_path=paths)
    actual = simulate(case, backend, tmp_path / "cpu.nc")
    comparison = compare(reference, actual, mode="deterministic", relative_tolerance=1e-10)
    assert comparison.passed, comparison.details
    assert check_result(actual) == []
    assert backend.metadata is None
    reloaded = SimulationResult.load(actual.path)
    assert compare(
        actual, reloaded, mode="deterministic", absolute_tolerance=0, relative_tolerance=0
    ).passed


def test_sessions_keep_rng_state_and_owned_results():
    with CPU() as first, CPU() as second:
        first.prepare(CASE)
        second.prepare(CASE)
        first_round = first.run_batch(BatchControls(1, 200))
        snapshot = first_round.stokes.copy()
        second_round = first.run_batch(BatchControls(2, 200))
        other_first = second.run_batch(BatchControls(1, 200))
        other_second = second.run_batch(BatchControls(2, 200))
        np.testing.assert_array_equal(snapshot, other_first.stokes)
        np.testing.assert_array_equal(second_round.stokes, other_second.stokes)
        assert not np.array_equal(snapshot, second_round.stokes)
        first.prepare(CASE)
        replay = first.run_batch(BatchControls(1, 200))
        np.testing.assert_array_equal(snapshot, replay.stokes)

    np.testing.assert_array_equal(first_round.stokes, snapshot)
    first.close()

    with pytest.raises(RuntimeError, match="prepared"):
        first.run_batch(BatchControls(1, 200))


def test_bad_requests_do_not_advance_random_stream():
    with CPU(threads=2) as backend:
        backend.prepare(CASE)

        with pytest.raises(ValueError, match="sequential"):
            backend.run_batch(BatchControls(2, 200))

        with pytest.raises(ValueError, match="divisible"):
            backend.run_batch(BatchControls(1, 201))

        assert backend.run_batch(BatchControls(1, 200)).controls.round_index == 1


def test_native_failure_invalidates_session():
    with CPU() as backend:
        backend.prepare(CASE)
        process = backend._process
        assert process is not None
        process.kill()
        process.wait()

        with pytest.raises((RuntimeError, BrokenPipeError)):
            backend.run_batch(BatchControls(1, 200))

        assert backend.metadata is None
        assert backend._process is None
        assert backend._temporary is None


def test_timeout_closes_native_worker(monkeypatch):
    backend = CPU()

    def timeout(marker):
        raise TimeoutError("test stalled native worker")

    monkeypatch.setattr(backend, "_wait", timeout)

    with pytest.raises(TimeoutError):
        backend.prepare(CASE)

    assert backend._process is None
    assert backend._temporary is None


def test_short_run_preserves_uncertainty_sentinel(tmp_path):
    result = simulate(replace(CASE, rounds=2), CPU(), tmp_path / "short.nc")
    assert np.all(result.standard_error == 999)
    assert result.values["radiance_convergence"].shape[0] == 2
    assert result.values["n_rounds"].item() == 2

    with pytest.raises(FileExistsError):
        simulate(CASE, CPU(), result.path)


def test_explicit_seed_replays_and_changes_photon_stream():
    """Independent-seed studies require a seed to affect transport reproducibly."""

    tallies = []

    for seed in (714103, 714103, 23274):
        with CPU(seed=seed) as backend:
            backend.prepare(CASE)
            batch = backend.run_batch(BatchControls(1, 200))
            tallies.append(batch)
            assert backend.provenance["seed"] == seed

    np.testing.assert_array_equal(tallies[0].stokes, tallies[1].stokes)
    np.testing.assert_array_equal(tallies[0].weighted_paths, tallies[1].weighted_paths)
    assert not np.array_equal(tallies[0].stokes, tallies[2].stokes)


@pytest.mark.parametrize("seed", [-1, 2**30 + 1, True, 1.5])
def test_invalid_cpu_seed(seed):
    with pytest.raises(ValueError, match="seed"):
        CPU(seed=seed)


def test_adaptive_run_matches_fixed_run_at_completed_count(tmp_path):
    from monki import AMFConvergence

    case = replace(CASE, photons_per_column=2000, rounds=96)
    adaptive = simulate(
        case,
        CPU(),
        tmp_path / "adaptive.nc",
        convergence=AMFConvergence(
            relative_error=0.9, minimum_rounds=32, check_every=16, consecutive_checks=2
        ),
    )
    count = int(adaptive.values["n_rounds"].item())
    assert count == 48
    fixed = simulate(replace(case, rounds=count), CPU(), tmp_path / "fixed.nc")

    for name in ("radiance", "AMF_layer", "std_radiance"):
        np.testing.assert_array_equal(adaptive.values[name], fixed.values[name])

    assert np.all(np.isfinite(adaptive.values["std_AMF_layer"]))
    assert np.all(adaptive.values["std_AMF_layer"] >= 0)


def test_adaptive_budget_exhaustion_is_not_convergence(tmp_path):
    import json

    from monki import AMFConvergence

    simulate(
        replace(CASE, rounds=32),
        CPU(),
        tmp_path / "limited.nc",
        convergence=AMFConvergence(relative_error=1e-12, minimum_rounds=32),
    )
    receipt = json.loads((tmp_path / "limited.run.json").read_text())
    assert receipt["status"] == "max_rounds_reached"
    assert receipt["completed_rounds"] == 32
    assert not receipt["amf_convergence"]["converged"]
