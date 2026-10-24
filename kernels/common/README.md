# Shared GPU transport

`transport.h` is the scientific source for Metal and CUDA. It retains the routine
mapping, float32 assumptions, Philox counter mapping, event estimator and failure
flags from the validated Metal implementation. Fortran prepares all optical tables.

The small platform adapters supply vectors, address spaces, math and atomics.
Each entry point launches one complete photon history per thread. Host bridges
own allocation, compilation, dispatch, synchronization and error reporting.

Python expands local includes before runtime compilation and hashes the complete
specialized source. Dimensions, view geometry and path diagnostics specialize each
session. Changes here require validation on both GPU platforms.

Philox integer words must match known answers exactly. Transcendental functions,
floating-point expressions and atomic order can differ across devices; transport
results require statistical validation against the independent Fortran reference.
