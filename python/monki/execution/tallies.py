"""Named batch requests and owned transport buffers shared by backend implementations."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class BatchControls:
    """Request one complete statistical round, numbered from one within a session."""

    round_index: int
    photons_per_column: int

    def __post_init__(self):
        for value in (self.round_index, self.photons_per_column):
            if type(value) is not int or not 0 < value < 2**31:
                raise ValueError("batch counts must be positive 32-bit integers")


@dataclass(frozen=True)
class RawTallies:
    """Own unnormalized Fortran-order Stokes and optional weighted-path sums.

    Stokes axes are (component, x, y, angle); paths are (x, y, layer, angle),
    including the vacuum layer; weights are (x, y, angle). No buffer aliases the worker.
    """

    controls: BatchControls
    stokes: np.ndarray
    weighted_paths: np.ndarray | None = None
    weights: np.ndarray | None = None
