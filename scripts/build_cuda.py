"""Build the CUDA driver bridge; NVRTC specializes transport during preparation."""

import hashlib
import json
import os
import platform
import subprocess
from pathlib import Path

from build_cpu import build_cpu

ROOT = Path(__file__).resolve().parents[1]


def build_cuda() -> Path:
    """Build on Linux with the CUDA toolkit; compilation needs no attached GPU."""

    if platform.system() != "Linux":
        raise RuntimeError("CUDA requires Linux")

    build_cpu()
    build_cpu(save_path=False)
    toolkit = Path(os.environ.get("CUDA_HOME", "/usr/local/cuda"))
    source = ROOT / "kernels/cuda/bridge.cpp"
    output = ROOT / "build/cuda/libmonki_cuda.so"
    output.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "c++",
        "-std=c++17",
        "-shared",
        "-fPIC",
        "-O2",
        "-Wall",
        "-Wextra",
        "-Werror",
        "-Wno-deprecated-declarations",
        f"-I{toolkit / 'include'}",
        f'-DMONKI_CUDA_INCLUDE="{toolkit / "include"}"',
        str(source),
        f"-L{toolkit / 'lib64'}",
        f"-L{toolkit / 'lib64/stubs'}",
        "-lcuda",
        "-lnvrtc",
        f"-Wl,-rpath,{toolkit / 'lib64'}",
        "-o",
        str(output),
    ]
    # The provider credential is only used by the cloud client.
    env = {key: value for key, value in os.environ.items() if key != "RUNPOD_API_KEY"}
    subprocess.run(command, check=True, env=env)
    receipt = {
        "command": command,
        "compiler": subprocess.check_output(["c++", "--version"], text=True, env=env),
        "bridge_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "binary_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
    }
    output.with_suffix(".json").write_text(json.dumps(receipt, indent=2) + "\n")
    return output


if __name__ == "__main__":
    print(build_cuda())
