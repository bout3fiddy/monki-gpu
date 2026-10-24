# Metal backend

Fortran prepares the optical tables; `transport.metal` dispatches photon histories from `kernels/common/transport.h`.
`bridge.m` owns Metal resources and dispatches. Python owns statistics and output.
CPU photon transport is never used as a fallback.

Requires macOS Command Line Tools and Metal 3 hardware. Build with
`uv run python scripts/build_metal.py`; shaders compile at runtime.

## Routine mapping

| Victor's Fortran | Metal |
|---|---|
| `monte_carlo` | `monte_carlo` |
| `ADVANCE_internal` | `advance` |
| `ray_scat_matrix`, `mie_scat_matrix` | `scattering_elements` |
| `stokes_rot_mat`, `matmul_stokes_rot` | `rotation`, `scatter_stokes` |
| `rayl_scatter_pol`, `mie_scatter_pol_bound` | `sample_direction` |
| `GROUND_LAMB2` | Lambertian step in `monte_carlo` |
| `calc_radiances_pol`, `calc_radiances_and_paths_pol` | `detector` |

Supports forward I/IQU, one aerosol species, flat Lambertian surfaces, periodic
3D grids and TOA source/detector/output. Backward tracing, Stokes V, multiple
species, other detector grids and scattering-order diagnostics are unsupported.

## Precision and memory

Transport and atomic tallies use float32 with fast math disabled; host reductions
and statistics use float64. Atomic addition order can vary. Philox4x32-10 uses
seed, round, photon ID and draw block; equal CPU/Metal seeds do not match paths.

Shared `MTLBuffer` allocations persist across rounds in Apple unified memory.
Completion precedes copying owned host results. Each session permits one caller.
Defaults are 65,536 photons per dispatch and 32 tally banks (allowed: 1–64).
Grouped rounds keep separate RNG identities and tallies, up to 32 slots and a
64 MiB reservation budget; one larger round remains supported.

Prepared dimensions specialize the shader. Uniform upper-atmosphere detector paths
use optical-depth prefixes; periodic 1 × 1 scenes skip horizontal boundaries.
Both retain physical slant lengths. `profile=True` records completed GPU command time.

## Native interfaces

`transport.bin` v1 is a native float64 stream in Fortran order:

1. Header: version, nx, ny, nz, nmat, angles, aerosol count, source layer,
   detector layer, output-grid layer, save_path, backward.
2. Geometry: xunit, yunit, depol, solar propagation zenith, solar azimuth, detector azimuth.
3. zunit, views, surface types, albedos, elevations, BETA, species fractions, gas albedos.
4. Per aerosol: LUT size, albedo, cosines, bound CDF, F11, F22, F33, F12.

`gpu_tables.py` validates and packs float32 buffers. Bridge ABI v2 takes a round
capacity and eight control words: photon count, first flattened photon, first
statistical round, seed, traversal limit, profiling flag, photons per round,
active round count. A stale bridge requires rebuilding.

Error flags: 1 nonfinite probability; 2 detector traversal limit; 4 sampling limit;
8 invalid Stokes; 16 history limit. Any failure rejects the round/group and closes
the session; histories are never silently truncated.

See [validation](../../validation/README.md) for GPU tests and CPU comparisons.
