"""Executed reference checks; fail with a useful message if the build is missing."""

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from monki import ReferenceCase, ReferenceRunner
from monki.validation import check_result, compare

ROOT = Path(__file__).resolve().parents[2]
UPSTREAM = ROOT / "fortran/monki/upstream"

pytestmark = pytest.mark.integration


def test_repeated_seeded_runs_reproduce_all_outputs(tmp_path):
    runner = ReferenceRunner(ROOT / "build/reference/paths/monki", UPSTREAM)
    case = ReferenceCase("clear", "testrun/MONKI_MLS_cloudfree_340_340_test.in")
    first = runner.run(case, tmp_path / "first")
    second = runner.run(case, tmp_path / "second")
    assert check_result(first) == []
    result = compare(
        first, second, mode="deterministic", absolute_tolerance=0, relative_tolerance=0
    )
    assert result.passed, result.details
    assert first.radiance.shape == (3, 1, 1, 3, 1)
    assert "AMF_layer" in first.values
    assert np.all(first.values["AMF_layer"] >= 0)


def test_backward_path_produces_valid_solar_angle_axis(tmp_path):
    runner = ReferenceRunner(ROOT / "build/reference/paths/monki", UPSTREAM)
    case = ReferenceCase("backward", "testrun/MONKI_MLS_cloudfree_340_340_test.in", backward=True)
    result = runner.run(case, tmp_path / "backward")
    assert check_result(result) == []
    assert result.dimensions["radiance"][0] == "nsza"


def test_threaded_rounds_are_statistically_consistent(tmp_path):
    runner = ReferenceRunner(ROOT / "build/reference/paths/monki", UPSTREAM)
    case = ReferenceCase(
        "clear", "testrun/MONKI_MLS_cloudfree_340_340_test.in", photons_per_column=2000
    )
    serial = runner.run(case, tmp_path / "serial")
    parallel = replace(runner, threads=2).run(case, tmp_path / "parallel")
    result = compare(serial, parallel)
    assert result.passed, result.details
