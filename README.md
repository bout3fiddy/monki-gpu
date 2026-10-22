# MONKI

Monte Carlo radiative transfer by Victor J.H. Trees, with Python interfaces,
local CPU and Apple Metal backends, NetCDF output and reference validation.
CPU transport reuses Victor's Fortran; Metal translates the transport kernels.
CUDA and cloud execution are not implemented.

## Setup

Requires Python 3.14+, uv, GNU Fortran and NetCDF Fortran. Metal also requires
macOS Command Line Tools and Metal 3 hardware.

```sh
brew install gcc netcdf-fortran
uv sync --locked
uv run python scripts/build_cpu.py
uv run python scripts/build_cpu.py --no-path
uv run python scripts/build_metal.py
```

On Debian/Ubuntu, install `gfortran libnetcdff-dev libnetcdf-dev` and omit Metal.
Run from this checkout; wheels do not bundle native binaries or scientific data.

## Run

Run Python scripts with `uv run python`:

```python
from monki import AMFConvergence, ReferenceCase, simulate
from monki.execution import CPU, Metal

case = ReferenceCase(
    name="clear",
    input_file="testrun/MONKI_MLS_cloudfree_340_340_test.in",
    photons_per_column=2_000,
    rounds=32,
)

result = simulate(case, backend=CPU(threads=4), output="results/clear-cpu.nc")
result = simulate(case, backend=Metal(), output="results/clear-metal.nc")
print(result.radiance)
print(result.standard_error)
```

Each call prepares one session, completes statistical rounds and writes NetCDF
plus a `.run.json` receipt. Existing outputs are never overwritten.
CPU photon counts per column must be divisible by the thread count.
Use `save_path=False` on either backend for the box-cloud fixture.

Supported inputs are the supplied 340 nm clear, layered-cloud and box-cloud
fixtures. CPU supports forward/backward tracing; Metal supports forward I or
I/Q/U, one aerosol species, flat Lambertian surfaces, periodic grids and a
source/detector/output grid at top of atmosphere. Unsupported options fail.
Build limits: 30 × 20 pixels, 35 layers including vacuum, three viewing angles.

## AMF stopping

Set `case.rounds` to the maximum budget and pass
`convergence=AMFConvergence(relative_error=0.01)` to `simulate`.
This requires AMF diagnostics (`save_path=True`).

The driver requires estimated SE / |AMF| < 1% in every physical layer, pixel and
view. Defaults: at least 128 rounds, checks every 32 rounds, three consecutive
passes. Zero AMFs cannot meet a relative-only criterion. NetCDF includes
`std_AMF_layer`; the receipt reports `converged` or `max_rounds_reached`.
Reaching the cap saves the result without claiming convergence.
This estimates sampling precision. Physical accuracy requires independent validation;
rare paths and repeated checks can underestimate uncertainty. Omit the policy for fixed-round execution.

## Checks and timings

```sh
uv run python scripts/verify.py
uv run python validation/benchmark.py --threads 1 4 --repeats 3
uv run python validation/benchmark_amf_convergence.py
```

The verification gate includes reference simulations, but no benchmarks.
Benchmark reports include correctness checks, timings and hardware details.
Add `--profile-metal` to the fixed-work benchmark for GPU command timings.

[Architecture](ARCHITECTURE.md) · [Code style](STYLEGUIDE.md) ·
[Validation and plots](validation/README.md) · [Reference and licenses](fortran/monki/README.md)
