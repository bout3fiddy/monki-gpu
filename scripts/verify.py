"""Canonical development gate: integrity, lint, build, invariant tests, reference suite."""

import platform
import subprocess
import sys
from pathlib import Path

from build_cpu import build_cpu
from build_metal import build_metal
from build_reference import build, verify_vendor

ROOT = Path(__file__).resolve().parents[1]


def main():
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

    subprocess.run([sys.executable, "-m", "pytest", "-q"], cwd=ROOT, check=True)
    suite_command = [sys.executable, "validation/run_suite.py"]

    if platform.system() == "Darwin":
        suite_command.append("--metal")

    subprocess.run(suite_command, cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
