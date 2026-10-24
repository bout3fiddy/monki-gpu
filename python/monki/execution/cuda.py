"""NVIDIA CUDA adapter for the shared GPU session lifecycle."""

from .gpu import ROOT, GPUBackend


class CUDA(GPUBackend):
    """Execute shared float32 transport on one NVIDIA GPU using NVRTC."""

    backend_name = "CUDA"
    native_prefix = "monki_cuda"
    required_platform = "Linux"
    default_library = ROOT / "build/cuda/libmonki_cuda.so"
    default_kernel = ROOT / "kernels/cuda/transport.cu"
