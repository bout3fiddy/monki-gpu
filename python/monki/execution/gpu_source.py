"""Assemble local shader includes so runtime compilation hashes the complete source."""

import re
from pathlib import Path

import numpy as np

from .gpu_tables import HEADER_SIZE, SPECIALIZATION_FIELDS, UNIFORM_TOP


def kernel_source(path: Path, ancestors: tuple[Path, ...] = ()) -> bytes:
    """Inline quoted includes; leave compiler-provided system headers intact."""

    path = path.resolve()

    if path in ancestors:
        raise ValueError(f"Cyclic GPU source include: {path.name}")

    def include(match: re.Match) -> str:
        child = path.parent / match.group(1)
        return kernel_source(child, (*ancestors, path)).decode()

    return re.sub(
        r'^\s*#include "([^"]+)"\s*$', include, path.read_text(), flags=re.MULTILINE
    ).encode()


def specialized_source(kernel: Path, packed: np.ndarray, nz: int, angles: int) -> bytes:
    """Fix dimensions and view geometry for a prepared GPU session."""

    defines = "".join(
        f"#define MONKI_{name} {int(packed[index])}\n"
        for index, name in enumerate(SPECIALIZATION_FIELDS)
    )
    defines += f"#define MONKI_UNIFORM_TOP {int(packed[UNIFORM_TOP])}\n"

    for name, value in zip(("THETA0", "PHI0", "PHI"), packed[10:13], strict=True):
        defines += f"#define MONKI_{name} {float(value)!r}f\n"

    for index, value in enumerate(packed[HEADER_SIZE + nz : HEADER_SIZE + nz + angles]):
        defines += f"#define MONKI_VIEW_{index} {float(value)!r}f\n"

    return defines.encode() + kernel_source(kernel)
