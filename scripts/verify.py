"""Canonical development gate: integrity, lint, build, invariant tests, reference suite."""

import argparse
import os
import platform
import subprocess
import sys
from pathlib import Path

from build_cpu import build_cpu
from build_cuda import build_cuda
from build_metal import build_metal
from build_reference import build, verify_vendor

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cuda", action="store_true", help="Require CUDA build and real-device tests"
    )
    args = parser.parse_args()
    verify_vendor()

    subprocess.run(["uv", "lock", "--check"], cwd=ROOT, check=True)
    subprocess.run(
        [
            sys.executable,
            "-m",
            "ruff",
            "format",
            "--check",
            "python",
            "scripts",
            "tests",
            "validation",
        ],
        cwd=ROOT,
        check=True,
    )
    subprocess.run(
        [sys.executable, "-m", "ruff", "check", "python", "scripts", "tests", "validation"],
        cwd=ROOT,
        check=True,
    )
    subprocess.run([sys.executable, "scripts/linting/check_readability.py"], cwd=ROOT, check=True)
    subprocess.run(["uv", "run", "basedpyright"], cwd=ROOT, check=True)

    subprocess.run(
        [
            "uv",
            "run",
            "clang-format",
            "--dry-run",
            "--Werror",
            "kernels/metal/bridge.m",
            "kernels/metal/transport.metal",
            "kernels/common/transport.h",
            "kernels/cuda/platform.cuh",
            "kernels/cuda/transport.cu",
            "kernels/cuda/bridge.cpp",
        ],
        cwd=ROOT,
        check=True,
    )

    build()
    build(save_path=False)
    build_cpu()
    build_cpu(save_path=False)

    if platform.system() == "Darwin":
        build_metal()

    if args.cuda:
        build_cuda()

    environment = os.environ.copy()
    environment.pop("RUNPOD_API_KEY", None)
    environment["MONKI_TEST_CUDA"] = "1" if args.cuda else "0"
    subprocess.run([sys.executable, "-m", "pytest", "-q"], cwd=ROOT, check=True, env=environment)
    suite_command = [sys.executable, "validation/run_suite.py"]

    if platform.system() == "Darwin":
        suite_command.append("--metal")

    if args.cuda:
        suite_command.append("--cuda")

    subprocess.run(suite_command, cwd=ROOT, check=True, env=environment)


if __name__ == "__main__":
    main()
