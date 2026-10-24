# NVIDIA CUDA backend

`CUDA` uses the shared Python GPU lifecycle and `kernels/common/transport.h`.
`platform.cuh` adapts vector math, address spaces and atomics. `transport.cu`
provides the launch entry point. `bridge.cpp` owns one CUDA context per session,
resident device allocations, synchronous copies, error checks and event timing.

Requires Linux, NVIDIA drivers and CUDA 12.8 toolkit headers/NVRTC. Python remains
3.14+. Build with `uv run python scripts/build_cuda.py`. The host bridge can be
built without a GPU. NVRTC specializes the shader for the selected device during
prepare; compiler/driver/architecture details and full source hashes appear in
run provenance. Runtime compilation is included in preparation timings.

Transport and atomic tallies use float32. NVRTC disables fused multiply-add,
flush-to-zero and fast division/sqrt; host shard and round reductions use float64.
Philox uses the same seed/round/photon/draw-block identity as Metal. Identical
random words do not guarantee identical floating-point histories across devices.

The supported physics and fixture limits match Metal: forward I/IQU, one aerosol
species, flat Lambertian surfaces, periodic grids and TOA source/detector/output.
Paths are enabled for clear/layered fixtures and disabled for box clouds.
Unsupported settings, unavailable devices, invalid tallies and truncated histories
fail explicitly. Failed groups invalidate the complete session.

Validation commands:

```sh
uv run python scripts/check_cuda_compile.py
uv run python scripts/verify.py --cuda
uv run python validation/benchmark.py --backend cuda --threads 1 4
```

The first command compiles all fixture specializations without executing GPU
transport. The CUDA verification gate requires a real GPU and runs the shared
helper/lifecycle tests plus CPU/reference comparisons for all three fixtures.
