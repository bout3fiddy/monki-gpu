"""Build the small Objective-C Metal bridge; kernels compile at session preparation."""

import hashlib
import json
import platform
import subprocess
from pathlib import Path

from build_cpu import build_cpu
from monki.environment import native_environment

ROOT = Path(__file__).resolve().parents[1]


def build_metal() -> Path:
    """Use the system SDK; no full Xcode Metal compiler or Python Metal package is needed."""

    if platform.system() != "Darwin":
        raise RuntimeError("Metal requires macOS")

    build_cpu()
    build_cpu(save_path=False)
    source = ROOT / "kernels/metal/bridge.m"
    output = ROOT / "build/metal/libmonki_metal.dylib"
    output.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "xcrun",
        "clang",
        "-dynamiclib",
        "-fobjc-arc",
        "-O2",
        "-Wall",
        "-Wextra",
        "-Werror",
        "-Wno-deprecated-declarations",
        "-framework",
        "Foundation",
        "-framework",
        "Metal",
        str(source),
        "-o",
        str(output),
    ]
    subprocess.run(command, check=True, env=native_environment())
    receipt = {
        "command": command,
        "bridge_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "binary_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
    }
    output.with_suffix(".json").write_text(json.dumps(receipt, indent=2) + "\n")
    return output


if __name__ == "__main__":
    print(build_metal())
