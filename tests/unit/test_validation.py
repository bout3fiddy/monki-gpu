"""Comparator coverage: plausible wrong results must not pass the gate."""

from dataclasses import replace
from pathlib import Path

import netCDF4
import numpy as np
import pytest
from monki import SimulationResult
from monki.validation import check_result, compare, select_views


@pytest.fixture
def result():
    values = {
        "radiance": np.array([[[[[0.5], [0.1], [-0.05]]]]]),
        "std_radiance": np.full((1, 1, 1, 3, 1), 0.001),
        "DOP": np.array([[[[np.hypot(0.1, 0.05) / 0.5]]]]),
        "wav": np.array([340.0]),
        "theta": np.array([30.0]),
        "n_rounds": np.array([32]),
        "n_photons_per_column": np.array([2000]),
    }
    dimensions = {
        "radiance": ("nvza", "ny", "nx", "nmat", "nwavs"),
        "std_radiance": ("nvza", "ny", "nx", "nmat", "nwavs"),
        "DOP": ("nvza", "ny", "nx", "nwavs"),
        "wav": ("nwavs",),
        "theta": ("nvza",),
        "n_rounds": ("single",),
        "n_photons_per_column": ("single",),
    }
    return SimulationResult(Path("synthetic.nc"), values, dimensions, {})


def changed(result, name, value):
    return replace(result, values={**result.values, name: value})


def test_identical_stokes_and_uncertainty_pass_both_modes(result):
    assert check_result(result) == []
    assert compare(result, result, mode="deterministic").passed
    assert compare(result, result).passed


def test_small_signed_stokes_variation_is_statistically_allowed(result):
    values = result.radiance.copy()
    values[..., 1, :] += 0.001
    candidate = changed(result, "radiance", values)
    candidate = changed(
        candidate, "DOP", np.linalg.norm(values[..., 1:, :], axis=3) / values[..., 0, :]
    )
    assert compare(result, candidate).passed
    assert not compare(result, candidate, mode="deterministic").passed


def test_wrong_but_physically_valid_radiance_fails(result):
    candidate = changed(result, "radiance", result.radiance * 1.5)
    assert check_result(candidate) == []  # Scaling preserves DOP and polarization bounds.
    comparison = compare(result, candidate)
    assert not comparison.passed
    assert comparison.failed_values > 0
    assert comparison.max_normalized_error > 1


@pytest.mark.parametrize(
    "name,value,message",
    [
        ("theta", np.array([31.0]), "coordinate mismatch"),
        ("wav", np.array([550.0]), "coordinate mismatch"),
        ("std_radiance", np.full((1, 1, 1, 3, 1), 999.0), "uncertainty unavailable"),
        ("std_radiance", np.full((1, 1, 1, 3, 1), -0.1), "negative standard error"),
        ("radiance", np.full((1, 1, 1, 3, 1), np.nan), "non-finite"),
        ("n_rounds", np.array([30]), "uncertainty unavailable"),
        ("DOP", np.array([[[[0.9]]]]), "DOP is inconsistent"),
    ],
)
def test_invalid_or_misaligned_results_fail(result, name, value, message):
    report = compare(result, changed(result, name, value))
    assert not report.passed
    assert any(message in item for item in report.details)


def test_same_shape_wrong_axis_names_fail(result):
    candidate = replace(
        result, dimensions={**result.dimensions, "radiance": ("nvza", "nx", "ny", "nmat", "nwavs")}
    )
    assert not compare(result, candidate).passed


def test_units_mismatch_fails_even_when_numbers_match(result):
    candidate = replace(result, attributes={"radiance": {"units": "W m-2 sr-1"}})
    comparison = compare(result, candidate)
    assert not comparison.passed
    assert "units mismatch: radiance" in comparison.details


def test_changed_optical_scene_cannot_pass_as_same_case(result):
    reference = changed(result, "bmsca", np.array([0.1]))
    reference = replace(reference, dimensions={**reference.dimensions, "bmsca": ("nz",)})
    candidate = changed(reference, "bmsca", np.array([0.2]))
    assert not compare(reference, candidate).passed


def test_single_plus_multiple_must_match_total(result):
    candidate = changed(result, "radiance_sing", result.radiance * 0.4)
    candidate = changed(candidate, "radiance_mult", result.radiance * 0.7)
    assert "single + multiple scattering does not equal total radiance" in check_result(candidate)


def test_negative_intensity_and_excess_polarization_rejected(result):
    values = result.radiance.copy()
    values[..., 0, :] = -0.1
    failures = check_result(changed(result, "radiance", values))
    assert "negative intensity" in failures
    assert "polarized magnitude exceeds intensity" in failures


def test_views_selected_by_coordinate_not_position(result):
    assert select_views(result, (30.0,)).values["theta"].item() == 30

    with pytest.raises(ValueError, match="absent"):
        select_views(result, (0.0,))


def test_masked_netcdf_is_rejected(tmp_path):
    path = tmp_path / "missing.nc"

    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("n", 2)
        variable = dataset.createVariable("radiance", "f8", ("n",), fill_value=-999)
        variable[0] = 1

    with pytest.raises(ValueError, match="masked/missing"):
        SimulationResult.load(path)


@pytest.mark.parametrize("rate", [0, 1, -1, float("nan")])
def test_invalid_statistical_policy_rejected(result, rate):
    with pytest.raises(ValueError):
        compare(result, result, family_error_rate=rate)
