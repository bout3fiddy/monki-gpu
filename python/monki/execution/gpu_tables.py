"""Validate and pack Victor's prepared transport tables for the GPU kernels."""

import numpy as np

# Version 1 retains the original Metal layout; the ABI fixes float32 and Fortran order.
LAYOUT_VERSION = 1
HEADER_SIZE = 16
TALLY_SHARDS = 14
UNIFORM_TOP = 15
SPECIALIZATION_FIELDS = ("NX", "NY", "NZ", "NMAT", "ANGLES", "LUT_SIZE", "PATHS")


def pack_tables(raw: np.ndarray) -> tuple[np.ndarray, tuple[int, int, int, int, int]]:
    """Convert the versioned float64 Fortran export to an explicit float32 layout.

    Preparation stays in Fortran. This adapter only validates supported physics,
    converts units for angles, and lays out buffers. Uniform-tail prefix sums
    integrate the supplied extinction; no optical properties are recomputed.
    """

    if raw.ndim != 1 or raw.size < 18 or not np.all(np.isfinite(raw)):
        raise ValueError("Invalid/nonfinite prepared transport export")

    version, nx, ny, nz, nmat, angles, species, source, detector, register, paths, backward = raw[
        :12
    ]

    if version != 1 or np.any(raw[:12] != raw[:12].astype(int)):
        raise ValueError("Unsupported transport export version or dimensions")

    nx, ny, nz, nmat, angles = map(int, (nx, ny, nz, nmat, angles))

    if not (1 <= nx <= 30 and 1 <= ny <= 20 and 2 <= nz <= 35 and 1 <= angles <= 3):
        raise ValueError("GPU scene exceeds validated grid limits")

    if (
        nmat not in (1, 3)
        or species != 1
        or backward
        or (source, detector, register) != (nz, nz, nz)
    ):
        raise NotImplementedError(
            "GPU supports forward I/IQU, one aerosol type, TOA source/detector/grid"
        )

    offset = 18

    def take(count):
        nonlocal offset
        result = raw[offset : offset + count]
        offset += count

        if result.size != count:
            raise ValueError("Truncated transport export")

        return result

    zunit = take(nz)
    views = take(angles)
    surface = take(nx * ny)
    albedo = take(nx * ny)
    elevation = take(nx * ny)
    beta = take(nx * ny * nz)
    fractions = take(2 * nx * ny * nz).reshape((-1, 2))
    ssalb = take(nz)
    lut_size, aerosol_albedo = take(2)

    if lut_size != int(lut_size) or not 2 <= lut_size <= 100_000:
        raise ValueError("Invalid scattering LUT size")

    lut_size = int(lut_size)
    tables = take(6 * lut_size).reshape((6, lut_size))

    if offset != raw.size or not np.all(surface == -1) or not np.all(elevation == 1):
        raise NotImplementedError("GPU supports only flat Lambertian surfaces and version-1 tables")

    if np.any(zunit <= 0) or np.any(beta < 0) or np.any(np.diff(tables[0]) <= 0):
        raise ValueError("Invalid layer dimensions, extinction, or scattering angle grid")

    if np.any(np.diff(tables[1]) < 0) or not np.allclose(tables[1, [0, -1]], [0, 1]):
        raise ValueError("Invalid scattering CDF")

    header = np.zeros(HEADER_SIZE)
    header[:14] = [
        nx,
        ny,
        nz,
        nmat,
        angles,
        lut_size,
        paths,
        raw[12],
        raw[13],
        raw[14],
        *np.deg2rad(raw[15:18]),
        aerosol_albedo,
    ]
    # Above the highest horizontally varying layer, detector attenuation depends
    # only on altitude. Integrate the same float32 extinction/thickness that the
    # kernel uses; this replaces repeated traversal without averaging any cells.
    beta32 = beta.astype(np.float32).reshape(nz, nx * ny)
    z32 = zunit.astype(np.float32)
    varying = np.flatnonzero(np.any(beta32 != beta32[:, :1], axis=1))
    header[UNIFORM_TOP] = int(varying[-1] + 1) if varying.size else 0
    height_above = np.zeros(nz)
    tau_above = np.zeros(nz)

    for layer in range(nz - 3, -1, -1):
        height_above[layer] = height_above[layer + 1] + float(z32[layer + 1])
        tau_above[layer] = tau_above[layer + 1] + float(z32[layer + 1]) * float(
            beta32[layer + 1, 0]
        )

    packed = np.concatenate(
        [
            header,
            zunit,
            np.deg2rad(views),
            albedo,
            beta,
            fractions[:, 0],
            ssalb,
            tables.ravel(),
            height_above,
            tau_above,
        ]
    ).astype(np.float32)

    if not np.all(np.isfinite(packed)):
        raise ValueError("Transport tables overflow GPU float32 precision")

    if np.any(np.diff(tables[0].astype(np.float32)) <= 0):
        raise ValueError("Scattering angles are not distinct at GPU float32 precision")

    if np.any(np.cos(np.deg2rad(views).astype(np.float32)) <= 0):
        raise ValueError("GPU viewing angles must remain above the horizon in float32")

    return packed, (nmat, nx, ny, nz, angles)
