# Validation

```sh
uv run python scripts/verify.py
uv run python validation/run_suite.py --case clear --metal
uv run python validation/compare.py reference.nc candidate.nc --mode statistical
```

The gate checks hashes, formatting, lint, types, native builds, tests and all three
scientific cases. Metal runs on macOS; other systems explicitly skip it.
Benchmarks are separate.

## Cases and thresholds

| Case | Grid | Photons | Layer AMFs |
|---|---|---:|---|
| Clear sky | 1 × 1 | 64,000 | Yes |
| Layered cloud | 1 × 1 | 64,000 | Yes |
| Box clouds | 30 × 20 | 3,840,000 | No |

All use 32 rounds, 340 nm, albedo 0.05, solar zenith 45°, azimuths 0° and viewing
zeniths 0°, 30°, 60°. Definitions: `cases/`.

- CPU versus independently built original: numeric outputs, `atol=rtol=1e-10`.
- CPU versus supplied archives: statistical radiance comparison.
- Metal versus CPU: statistical radiance and enabled AMFs.
- Invariants: dimensions, coordinates, units, finite values, nonnegative intensity,
  polarization bounds, DOP consistency and valid standard errors.
- Tests also cover native helpers, RNG, threaded/backward CPU runs, ownership,
  failure handling and adaptive covariance/stopping.

Statistical radiance bound:

```text
abs(candidate - reference) <= 1e-10 + z * hypot(SE_candidate, SE_reference)
```

`z` uses a two-sided normal threshold with Bonferroni family error 1e-6.
`std_radiance` is SE of the round mean; the `999` sentinel through round 30 is rejected.
AMF SE uses round residuals `P - ratio*W`, retaining path/weight covariance,
with the same combined-error policy. These checks estimate finite-sample agreement. Physical accuracy requires
independent validation. Rare paths and repeated stopping can underestimate uncertainty.
Archived compiler/RNG details are unknown. Tolerances are not fitted to failures.

## Reports and plots

Unique `runs/` directories hold outputs, hashes, settings and photon counts.
Reference outputs are never overwritten. Adaptive receipts record every check
and distinguish `converged` from `max_rounds_reached`.

| Retained directory | Evidence |
|---|---|
| `20261022T120814640769Z` | CPU/reference/Metal comparisons and plots |
| `benchmark-20261022T120952180298Z` | Three-repeat fixed-work timings |
| `amf-study-20261022T130245Z` | Independent seeds and raw AMF round tallies |
| `amf-timing-20261022T135836Z` | Three-seed adaptive timings and archive checks |
| `20261022T142231760983Z` | Validation after repository cleanup |

```sh
uv run python validation/plot_results.py --run validation/runs/20261022T120814640769Z
uv run python validation/investigate_amf.py --plot-run validation/runs/amf-study-20261022T130245Z --benchmark validation/runs/benchmark-20261022T120952180298Z/report.json
```

These read saved results and write PNG/SVG. Reflectance is stored radiance divided
by cos(solar zenith), using MONKI's E0 = pi normalization. Difference maps use
combined SE units; AMF bands use retained round statistics.

## Benchmarks

```sh
uv run python validation/benchmark.py --threads 1 4 --repeats 3
uv run python validation/investigate_amf.py
uv run python validation/benchmark_amf_convergence.py
```

Timings include preparation, transport, statistics, output and cleanup; compilation
and warmup are excluded. Fixed-work repeats reuse streams; independent-seed studies
test sampling variation. Add `--profile-metal` for GPU command timings.

Adaptive timing uses four CPU threads, 8,000 photons/round, a 1% relative-SE target,
128 minimum rounds, checks every 32 rounds, three passes and an 8,192-round cap.
All twelve retained Apple M4 runs converged and passed archived-radiance checks.
Median CPU/Metal times: clear 2.074/0.466 s; layered cloud 232.839/26.160 s.
Maximum paired AMF difference: 3.493 combined estimated SE. Per-seed counts and
results are in the retained reports.
