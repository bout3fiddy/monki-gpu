"""Owned NetCDF result arrays, independent of executable lifetime."""

from dataclasses import dataclass
from pathlib import Path

import netCDF4
import numpy as np


@dataclass(frozen=True)
class SimulationResult:
    """Own data, dimension names, and units read from one MONKI NetCDF product."""

    path: Path
    values: dict[str, np.ndarray]
    dimensions: dict[str, tuple[str, ...]]
    attributes: dict[str, dict]

    @classmethod
    def load(cls, path: str | Path):
        """Copy arrays before closing the file; reject missing values immediately."""

        path = Path(path).resolve()

        with netCDF4.Dataset(path) as dataset:
            values = {}
            dimensions = {}
            attributes = {}

            for name, variable in dataset.variables.items():
                # Reject masked NetCDF values because they represent missing output.
                data = variable[:]

                if np.ma.isMaskedArray(data) and np.any(np.ma.getmaskarray(data)):
                    raise ValueError(f"{name}: masked/missing data in {path}")

                values[name] = np.array(data, copy=True)
                dimensions[name] = tuple(variable.dimensions)
                attributes[name] = {key: variable.getncattr(key) for key in variable.ncattrs()}

        return cls(path, values, dimensions, attributes)

    def save(self) -> None:
        """Write owned arrays without replacing an existing scientific result."""

        self.path.parent.mkdir(parents=True, exist_ok=True)

        # Exclusive creation protects supplied outputs, even if another writer
        # races the caller's initial path check.
        with netCDF4.Dataset(self.path, "w", clobber=False) as dataset:
            sizes = {}

            for name, value in self.values.items():
                for dimension, size in zip(self.dimensions[name], value.shape, strict=True):
                    if dimension in sizes and sizes[dimension] != size:
                        raise ValueError(f"Inconsistent dimension {dimension}")

                    sizes[dimension] = size

            for name, size in sizes.items():
                dataset.createDimension(name, size)

            for name, value in self.values.items():
                variable = dataset.createVariable(name, value.dtype, self.dimensions[name])
                variable.setncatts(self.attributes.get(name, {}))
                variable[:] = value

    @property
    def radiance(self):
        return self.values["radiance"]

    @property
    def standard_error(self):
        return self.values["std_radiance"]
