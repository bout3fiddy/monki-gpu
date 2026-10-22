"""Input and executable failure handling."""

from pathlib import Path

import pytest
from monki import ReferenceCase, ReferenceRunner
from monki.reference import fixture_text

ROOT = Path(__file__).resolve().parents[2]
UPSTREAM = ROOT / "fortran/monki/upstream"
INPUT = "testrun/MONKI_MLS_cloudfree_340_340_test.in"


@pytest.mark.parametrize(
    "arguments",
    [
        {"rounds": 0},
        {"rounds": True},
        {"photons_per_column": -1},
        {"view_zeniths": ()},
        {"view_zeniths": (90,)},
        {"surface_albedo": 1.1},
        {"solar_zenith": float("nan")},
        {"input_file": "../escape.in"},
        {"name": "../escape"},
        {"backward": 1},
    ],
)
def test_invalid_case_rejected(arguments):
    with pytest.raises(ValueError):
        ReferenceCase(**({"name": "clear", "input_file": INPUT} | arguments))


def test_only_intended_positional_fields_change():
    source = UPSTREAM / "INPUT" / INPUT
    original = source.read_text().splitlines()
    generated, dimensions = fixture_text(source, ReferenceCase("clear", INPUT))
    differences = {
        i for i, (a, b) in enumerate(zip(original, generated.splitlines(), strict=True)) if a != b
    }
    assert differences <= {16, 17, 18}
    assert dimensions == (1, 1, 33)


def test_unsupported_input_version_rejected():
    with pytest.raises(ValueError, match="Only the supplied testrun"):
        fixture_text(
            Path("unsupported/case.in"),
            ReferenceCase("unsupported", "unsupported/case.in"),
        )


def test_zero_exit_without_output_is_failure(tmp_path):
    executable = tmp_path / "fake-solver"
    executable.write_text("#!/bin/sh\necho 'Dimension error'\nexit 0\n")
    executable.chmod(0o755)
    runner = ReferenceRunner(executable, UPSTREAM)

    with pytest.raises(RuntimeError, match="exactly one NetCDF"):
        runner.run(ReferenceCase("clear", INPUT), tmp_path / "run")

    assert '"status": "failed"' in (tmp_path / "run/run.json").read_text()


def test_non_divisible_thread_partition_rejected(tmp_path):
    executable = tmp_path / "fake"
    executable.touch()

    with pytest.raises(ValueError, match="divisible"):
        ReferenceRunner(executable, UPSTREAM, threads=3).run(
            ReferenceCase("clear", INPUT, photons_per_column=200), tmp_path / "run"
        )

    assert not (tmp_path / "run").exists()
