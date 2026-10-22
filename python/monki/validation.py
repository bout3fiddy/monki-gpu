"""Scientific output checks and explicit deterministic/statistical comparisons."""

from dataclasses import asdict, dataclass
from statistics import NormalDist

import numpy as np

from .output import SimulationResult

COORDINATES = ("wav", "theta", "theta0", "phi", "phi0", "asurf", "x_unit", "y_unit", "z_unit")
SCENE_VARIABLES = ("pres", "temp", "alt", "bmabs", "bmsca", "baer", "vmrgas")

# Physical polarization guard from mc.f90's stokes_criterium_error.
STOKES_PHYSICAL_TOLERANCE = 1e-2
COORDINATE_TOLERANCE = 1e-12
UNCERTAINTY_SENTINEL = 999
MINIMUM_UNCERTAINTY_ROUNDS = 31


def select_views(result: SimulationResult, angles: tuple[float, ...]) -> SimulationResult:
    """Explicitly select exact forward viewing coordinates from a retained result."""

    if "theta" not in result.values or result.dimensions["theta"] != ("nvza",):
        raise ValueError("view selection requires forward theta(nvza)")

    indices = []

    for angle in angles:
        matches = np.flatnonzero(np.isclose(result.values["theta"], angle, rtol=0, atol=1e-12))

        if len(matches) != 1:
            raise ValueError(f"view {angle} is absent or ambiguous")

        indices.append(int(matches[0]))

    values = {
        name: np.take(data, indices, axis=result.dimensions[name].index("nvza"))
        if "nvza" in result.dimensions[name]
        else data.copy()
        for name, data in result.values.items()
    }
    return SimulationResult(result.path, values, result.dimensions.copy(), result.attributes.copy())


@dataclass(frozen=True)
class Comparison:
    passed: bool
    mode: str
    compared_values: int
    failed_values: int
    max_absolute_error: float
    max_normalized_error: float
    details: tuple[str, ...]

    def to_dict(self):
        return asdict(self)


def check_result(result: SimulationResult, *, require_uncertainty: bool = True) -> list[str]:
    """Check stored Stokes geometry, positivity, polarization, and available uncertainty."""

    failures = []

    for name in ("radiance", "std_radiance", "DOP", "wav", "n_rounds", "n_photons_per_column"):
        if name not in result.values:
            failures.append(f"missing variable: {name}")

    if failures:
        return failures

    for name, data in result.values.items():
        if np.issubdtype(data.dtype, np.number) and not np.all(np.isfinite(data)):
            failures.append(f"non-finite values: {name}")

    if "nmat" not in result.dimensions["radiance"]:
        return failures + ["radiance has no nmat dimension"]

    # NetCDF's Python axis order differs from the Fortran declarations. Resolve the
    # Stokes axis by name.
    axis = result.dimensions["radiance"].index("nmat")
    stokes = np.moveaxis(result.radiance, axis, -1)

    if stokes.shape[-1] not in (1, 3, 4):
        failures.append("Stokes vector must contain 1, 3, or 4 elements")

    intensity = stokes[..., 0]

    if np.any(intensity < -1e-12):
        failures.append("negative intensity")

    polarized = np.linalg.norm(stokes[..., 1:], axis=-1)

    # Match upstream's absolute 1e-2 Stokes acceptance tolerance; this is not a fit tolerance.
    if np.any(polarized > intensity + STOKES_PHYSICAL_TOLERANCE):
        failures.append("polarized magnitude exceeds intensity")

    # DOP is derived from I/Q/U/V. Checking this catches inconsistent products even
    # when each array separately looks finite and plausible.
    expected_dop = np.divide(
        polarized, intensity, out=np.zeros_like(intensity), where=intensity > 0
    )

    if result.values["DOP"].shape != intensity.shape or not np.allclose(
        result.values["DOP"], expected_dop, rtol=1e-9, atol=1e-10
    ):
        failures.append("DOP is inconsistent with Stokes components")

    error = result.standard_error

    if (
        error.shape != result.radiance.shape
        or result.dimensions["std_radiance"] != result.dimensions["radiance"]
    ):
        failures.append("standard-error shape/dimensions differ from radiance")

    if np.any(error < 0):
        failures.append("negative standard error")

    if require_uncertainty and (
        int(result.values["n_rounds"].item()) < MINIMUM_UNCERTAINTY_ROUNDS
        or np.any(error == UNCERTAINTY_SENTINEL)
    ):
        failures.append("uncertainty unavailable: need more than 30 completed rounds")

    if np.any(result.values["wav"] <= 0):
        failures.append("non-positive wavelength")

    for name in ("n_rounds", "n_photons_per_column"):
        if np.any(result.values[name] <= 0):
            failures.append(f"non-positive {name}")

    if "radiance_sing" in result.values or "radiance_mult" in result.values:
        if not {"radiance_sing", "radiance_mult"} <= result.values.keys():
            failures.append("incomplete scattering decomposition")
        elif not np.allclose(
            result.radiance,
            result.values["radiance_sing"] + result.values["radiance_mult"],
            rtol=1e-9,
            atol=1e-10,
        ):
            failures.append("single + multiple scattering does not equal total radiance")

    return failures


def compare(
    reference: SimulationResult,
    candidate: SimulationResult,
    *,
    mode: str = "statistical",
    absolute_tolerance: float = 1e-10,
    relative_tolerance: float = 1e-10,
    family_error_rate: float = 1e-6,
) -> Comparison:
    """Compare aligned outputs. Statistical mode checks radiance.

    Statistical bounds use combined standard errors and a Bonferroni normal threshold.
    Monte Carlo tails require multi-seed study to assess equivalence.
    Deterministic mode compares all numeric variables and requires the same run configuration.
    """

    if mode not in {"statistical", "deterministic"}:
        raise ValueError("mode must be statistical or deterministic")

    if not all(np.isfinite(x) and x >= 0 for x in (absolute_tolerance, relative_tolerance)):
        raise ValueError("tolerances must be finite and non-negative")

    if not 0 < family_error_rate < 1:
        raise ValueError("family_error_rate must be in (0, 1)")

    # Reject malformed or misaligned products before doing numerical comparisons.
    # Equal array sizes alone do not establish equal wavelengths, geometry, or units.
    problems = [
        f"reference: {x}"
        for x in check_result(reference, require_uncertainty=mode == "statistical")
    ]
    problems += [
        f"candidate: {x}"
        for x in check_result(candidate, require_uncertainty=mode == "statistical")
    ]

    if set(reference.values) != set(candidate.values):
        problems.append("variable inventories differ")

    for name in reference.values.keys() & candidate.values.keys():
        different_shape = (
            reference.dimensions[name] != candidate.dimensions[name]
            or reference.values[name].shape != candidate.values[name].shape
        )
        # Round convergence histories may differ in statistical mode.
        allow_round_history_difference = mode == "statistical" and name == "radiance_convergence"

        if different_shape and not allow_round_history_difference:
            problems.append(f"shape/dimension mismatch: {name}")

        reference_units = reference.attributes.get(name, {}).get("units")
        candidate_units = candidate.attributes.get(name, {}).get("units")

        if reference_units != candidate_units:
            problems.append(f"units mismatch: {name}")

    for name in (*COORDINATES, *SCENE_VARIABLES):
        if name in reference.values and name in candidate.values:
            a, b = reference.values[name], candidate.values[name]

            if a.shape != b.shape or not np.allclose(a, b, rtol=0, atol=COORDINATE_TOLERANCE):
                problems.append(f"coordinate mismatch: {name}")

    if problems:
        return Comparison(False, mode, 0, 0, 0.0, 0.0, tuple(problems))

    # Upstream supplies the radiance standard errors used for this comparison.
    # Statistical acceptance is deliberately limited to that measured quantity.
    names = ["radiance"] if mode == "statistical" else sorted(reference.values)
    count = failed = 0
    max_absolute = max_normalized = 0.0

    for name in names:
        a, b = reference.values[name], candidate.values[name]

        if not np.issubdtype(a.dtype, np.number):
            if not np.array_equal(a, b):
                problems.append(f"non-numeric values differ: {name}")

            continue

        delta = np.abs(a.astype(float) - b.astype(float))

        if mode == "statistical":
            # Independent estimates add variances. The Bonferroni tail accounts for
            # checking many pixels/components; it does not assume independent pixels.
            # Normal-tail calibration is approximate for finite Monte Carlo rounds.
            tail = max(family_error_rate / (2 * a.size), np.finfo(float).eps)
            z = NormalDist().inv_cdf(1 - tail)
            bound = absolute_tolerance + z * np.hypot(
                reference.standard_error, candidate.standard_error
            )
        else:
            bound = absolute_tolerance + relative_tolerance * np.abs(a)

        bad = delta > bound
        count += a.size
        failed += int(np.count_nonzero(bad))
        max_absolute = max(max_absolute, float(np.max(delta, initial=0)))
        normalized = delta / np.maximum(bound, np.finfo(float).tiny)
        max_normalized = max(max_normalized, float(np.max(normalized, initial=0)))

        if np.any(bad):
            problems.append(f"{name}: {np.count_nonzero(bad)}/{a.size} values exceed tolerance")

    return Comparison(
        not problems, mode, count, failed, max_absolute, max_normalized, tuple(problems)
    )
