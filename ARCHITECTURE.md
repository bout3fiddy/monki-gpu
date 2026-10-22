# Architecture

Python owns inputs, session lifetime, statistics, convergence and NetCDF output.
Victor's retained Fortran sources remain unchanged in `fortran/monki/upstream/`.

## Execution

`Backend` defines `prepare`, `run_batch` and `close`. `run_batches` defaults to
sequential rounds; Metal can group rounds with separate tallies and random streams.
`CPU` and `Metal` implement this interface. The driver does not switch on backend type.

1. `ReferenceCase` validates inputs and controls.
2. `prepare` loads the scene and optical tables once.
3. `run_batches` returns owned raw Stokes, weighted-path and weight arrays.
4. `simulate` computes radiance, standard errors, polarization and layer AMFs;
   optional `AMFConvergence` controls stopping.
5. NetCDF and a JSON receipt record results, counts, settings and hashes.
6. `close` releases resources, including after failure.

[CPU](kernels/cpu/README.md) uses a persistent Fortran process and OpenMP.
[Metal](kernels/metal/README.md) reuses Fortran preparation, retains GPU buffers
and dispatches translated transport through an Objective-C bridge.
Each session allows one caller; returned arrays survive later batches and closure.

CPU transport is float64 with the compiler RNG. Metal uses float32 with Philox;
host statistics use float64. Equal seeds do not produce matching photon paths.
AMFs divide accumulated weighted paths by weights and layer thickness;
uncertainty retains path/weight covariance. See [validation](validation/README.md).
