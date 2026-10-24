"""Provider-independent lifecycle for complete remote simulation jobs."""

import math
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import StrEnum


class JobState(StrEnum):
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"
    TIMED_OUT = "timed_out"


@dataclass(frozen=True)
class JobStatus:
    job_id: str
    state: JobState
    output: dict | None = None
    error: str | None = None


class Provider(ABC):
    """Submit portable JSON jobs and retrieve portable results.

    Implementations own authentication and provider state translation. Scientific
    backends execute inside the worker and never receive provider credentials.
    """

    @abstractmethod
    def submit(self, request: dict) -> str:
        """Submit once and return a job identifier; ambiguous failures must be surfaced."""

        raise NotImplementedError

    @abstractmethod
    def status(self, job_id: str) -> JobStatus:
        """Retrieve state and completed output for an existing job."""

        raise NotImplementedError

    @abstractmethod
    def cancel(self, job_id: str) -> None:
        """Request cancellation; use status to confirm the final state."""

        raise NotImplementedError

    def wait(self, job_id: str, *, timeout: float = 900, poll_interval: float = 2) -> dict:
        """Poll with bounded backoff; a local timeout leaves the job retrievable."""

        if not all(math.isfinite(value) and value > 0 for value in (timeout, poll_interval)):
            raise ValueError("timeout and poll_interval must be positive")

        deadline = time.monotonic() + timeout

        while True:
            status = self.status(job_id)

            if status.state == JobState.COMPLETED:
                if status.output is None:
                    raise RuntimeError("Completed job omitted its output")

                return status.output

            if status.state not in (JobState.QUEUED, JobState.RUNNING):
                raise RuntimeError(f"Job {job_id}: {status.state}: {status.error or 'no detail'}")

            remaining = deadline - time.monotonic()

            if remaining <= 0:
                raise TimeoutError(f"Local wait expired; job {job_id} may still be running")

            time.sleep(min(poll_interval, remaining))
            poll_interval = min(poll_interval * 1.5, 15)
