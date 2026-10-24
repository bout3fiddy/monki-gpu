"""Apple Metal adapter for the shared GPU session lifecycle."""

from .gpu import ROOT, GPUBackend


class Metal(GPUBackend):
    """Execute shared transport through the native Metal bridge."""

    backend_name = "Metal"
    native_prefix = "monki_metal"
    required_platform = "Darwin"
    default_library = ROOT / "build/metal/libmonki_metal.dylib"
    default_kernel = ROOT / "kernels/metal/transport.metal"
