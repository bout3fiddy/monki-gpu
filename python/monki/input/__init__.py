"""Validated reference-run inputs; the complete Scene API is not implemented yet."""

from dataclasses import dataclass
from math import isfinite
from pathlib import Path


@dataclass(frozen=True)
class ReferenceCase:
    """Specify a supported legacy fixture plus explicit Monte Carlo and geometry settings."""

    name: str
    input_file: str
    photons_per_column: int = 200
    rounds: int = 32
    surface_albedo: float = 0.05
    solar_zenith: float = 45.0
    solar_azimuth: float = 0.0
    view_azimuth: float = 0.0
    view_zeniths: tuple[float, ...] = (0.0, 30.0, 60.0)
    backward: bool = False

    def __post_init__(self):
        """Reject invalid requests before they reach Fortran's positional parser."""

        if not self.name or not self.name.replace("-", "").replace("_", "").isalnum():
            raise ValueError("case name must contain only letters, numbers, hyphens, underscores")

        path = Path(self.input_file)

        if path.is_absolute() or ".." in path.parts:
            raise ValueError("input_file must be relative to the vendored INPUT directory")

        # The legacy input/output count fields are 32-bit even though some internal
        # counters are wider. Reject oversized requests to prevent overflow.
        for name in ("photons_per_column", "rounds"):
            value = getattr(self, name)

            if type(value) is not int or not 0 < value < 2**31:
                raise ValueError(f"{name} must be a positive 32-bit integer")

        if not isfinite(self.surface_albedo) or not 0 <= self.surface_albedo <= 1:
            raise ValueError("surface_albedo must be in [0, 1]")

        if not self.view_zeniths or len(self.view_zeniths) > 91:
            raise ValueError("supply 1 to 91 viewing zenith angles")

        for angle in (self.solar_zenith, *self.view_zeniths):
            if not isfinite(angle) or not 0 <= angle < 90:
                raise ValueError("reference harness currently supports zeniths in [0, 90)")

        for angle in (self.solar_azimuth, self.view_azimuth):
            if not isfinite(angle) or not 0 <= angle < 360:
                raise ValueError("azimuths must be in [0, 360)")

        if type(self.backward) is not bool:
            raise ValueError("backward must be a bool")
