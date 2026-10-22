# CPU worker

`worker.f90` adapts Victor J.H. Trees's `main.f90` (MIT). Transport and helper
routines link directly from the unchanged reference sources. The worker retains
scene preparation, tally resets, OpenMP reduction and `monte_carlo`; Python owns
round scheduling, statistics and final output.

Build both profiles with `uv run python scripts/build_cpu.py` and the same command
with `--no-path` for radiance-only runs.

## Protocol v1

- Ready: `MONKI_CPU_READY 1 nx ny nz nmat nangles save_path threads`.
- Request on stdin: `round_index photons_per_column` (round indices start at one).
- Completed: `MONKI_CPU_BATCH round_index photons_per_column`.
- EOF or zero photons exits. Python copies `batch.bin` before another request.

`batch.bin` is native float64 without record markers, in Fortran order:

| Array | Shape |
|---|---|
| Stokes | `(nmat, nx, ny, nangles)` |
| Weighted paths, if enabled | `(nx, ny, nz, nangles)` |
| Detected weights, if enabled | `(nx, ny, nangles)` |

`nz` includes vacuum. Tallies include the forward geometry factor but remain raw
sums; normalization and AMF division happen after accumulation.
`OUTPUT/scene.nc` supplies metadata only: its zero/999 outputs are placeholders.
Its configured 32 rounds avoid the upstream interactive prompt; Python sets the
actual round count and history in the final output.

Preparation resets the configured seed. RNG streams depend on compiler and thread
count. Each session permits one caller; errors invalidate it and retain diagnostic
logs. Separate processes isolate legacy module state and Fortran `STOP`.
See [validation](../../validation/README.md) for reference checks.
