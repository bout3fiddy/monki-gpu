"""Runpod Serverless entry point; each request executes a complete CUDA job."""

import os
import signal
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryFile

from ..environment import native_environment
from .jobs import execute


def handler(job: dict) -> dict:
    request = job["input"]

    if request == {"version": 1, "operation": "validate"}:
        # A fixed diagnostic operation allows the canonical NVIDIA gate to run on
        # Serverless workers without granting remote command execution.
        root = Path(__file__).resolve().parents[3]

        with TemporaryFile() as log:
            process = subprocess.Popen(
                [sys.executable, "scripts/verify.py", "--cuda"],
                cwd=root,
                env=native_environment(),
                stdout=log,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
            timed_out = False

            try:
                process.wait(timeout=1200)
            except subprocess.TimeoutExpired:
                timed_out = True
            finally:
                if process.poll() is None:
                    # Also stop compiler, pytest and native descendants on timeout
                    # so a reusable Serverless worker cannot retain abandoned work.
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()

            log.seek(0, 2)
            log.seek(max(0, log.tell() - 64_000))
            return {
                "version": 1,
                "validation": {
                    "passed": process.returncode == 0 and not timed_out,
                    "timed_out": timed_out,
                    "log": log.read().decode(errors="replace"),
                },
            }

    return execute(request)


def main() -> None:
    import runpod

    # The synchronous handler owns one simulation session at a time per process.
    runpod.serverless.start({"handler": handler})


if __name__ == "__main__":
    main()
