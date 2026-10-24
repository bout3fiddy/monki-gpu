"""Compile production CUDA shaders with NVRTC without executing GPU transport."""

import ctypes
import json
from pathlib import Path

from monki.execution import CPU, CUDA
from monki.execution.gpu_source import specialized_source
from monki.execution.gpu_tables import pack_tables
from monki.input import ReferenceCase

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    """Check all fixture specializations; this supplies compiler evidence only."""

    native = ctypes.CDLL("libnvrtc.so.12")
    program_type = ctypes.c_void_p
    native.nvrtcCreateProgram.argtypes = [
        ctypes.POINTER(program_type),
        ctypes.c_char_p,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_void_p,
        ctypes.c_void_p,
    ]
    native.nvrtcCompileProgram.argtypes = [
        program_type,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_char_p),
    ]
    native.nvrtcGetProgramLogSize.argtypes = [program_type, ctypes.POINTER(ctypes.c_size_t)]
    native.nvrtcGetProgramLog.argtypes = [program_type, ctypes.c_char_p]
    native.nvrtcDestroyProgram.argtypes = [ctypes.POINTER(program_type)]

    for name in ("clear", "layered-cloud", "box-clouds"):
        case = ReferenceCase(**json.loads((ROOT / f"validation/cases/{name}.json").read_text()))

        with CPU(save_path=name != "box-clouds") as preparer:
            preparer.prepare(case)
            packed, (_, _, _, nz, angles) = pack_tables(preparer.transport_tables())

        source = specialized_source(CUDA.default_kernel, packed, nz, angles)
        program = program_type()
        status = native.nvrtcCreateProgram(
            ctypes.byref(program), source, b"monki.cu", 0, None, None
        )

        if status:
            raise RuntimeError(f"NVRTC create failed: {status}")

        try:
            options = (ctypes.c_char_p * 7)(
                b"--std=c++17",
                b"--fmad=false",
                b"--prec-div=true",
                b"--prec-sqrt=true",
                b"--ftz=false",
                b"--include-path=/usr/local/cuda/include",
                b"--gpu-architecture=compute_89",
            )
            status = native.nvrtcCompileProgram(program, len(options), options)
            size = ctypes.c_size_t()
            native.nvrtcGetProgramLogSize(program, ctypes.byref(size))
            log = ctypes.create_string_buffer(size.value)
            native.nvrtcGetProgramLog(program, log)

            if status:
                raise RuntimeError(log.value.decode())

            print(f"{name}: CUDA compilation passed (no GPU execution)")
        finally:
            native.nvrtcDestroyProgram(ctypes.byref(program))


if __name__ == "__main__":
    main()
