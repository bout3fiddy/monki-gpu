"""Runpod Serverless queue API adapter with environment-provided credentials."""

import json
import os
import re
from urllib.error import HTTPError, URLError
from urllib.request import HTTPRedirectHandler, Request, build_opener

from .provider import JobState, JobStatus, Provider

STATES = {
    "IN_QUEUE": JobState.QUEUED,
    "IN_PROGRESS": JobState.RUNNING,
    "COMPLETED": JobState.COMPLETED,
    "FAILED": JobState.FAILED,
    "CANCELLED": JobState.CANCELLED,
    "TIMED_OUT": JobState.TIMED_OUT,
}


class NoRedirect(HTTPRedirectHandler):
    """Keep the authorization header confined to the configured Runpod API origin."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


class RunpodServerless(Provider):
    """Use an existing queue endpoint; the API key stays in the client process."""

    def __init__(self, endpoint_id: str, *, execution_timeout: int = 600, ttl: int = 3600):
        self.endpoint_id = identifier(endpoint_id)

        for name, value, minimum in (("execution_timeout", execution_timeout, 5), ("ttl", ttl, 10)):
            if type(value) is not int or not minimum <= value <= 604800:
                raise ValueError(f"{name} must be an integer in [{minimum}, 604800] seconds")

        if ttl < execution_timeout:
            raise ValueError("ttl must cover execution_timeout and queue time")

        self.execution_timeout = execution_timeout
        self.ttl = ttl

    def submit(self, request: dict) -> str:
        response = self._request(
            "POST",
            "run",
            {
                "input": request,
                "policy": {
                    "executionTimeout": self.execution_timeout * 1000,
                    "ttl": self.ttl * 1000,
                },
            },
        )
        return identifier(response["id"])

    def status(self, job_id: str) -> JobStatus:
        job_id = identifier(job_id)
        response = self._request("GET", f"status/{job_id}")

        if response.get("id") != job_id or response.get("status") not in STATES:
            raise RuntimeError("Runpod returned an invalid job status")

        output = response.get("output")

        if output is not None and not isinstance(output, dict):
            raise RuntimeError("Runpod returned invalid MONKI output")

        # Do not echo arbitrary server messages containing request data or secrets.
        error = "Worker failed; inspect endpoint logs" if response.get("error") else None
        return JobStatus(job_id, STATES[response["status"]], output, error)

    def cancel(self, job_id: str) -> None:
        self._request("POST", f"cancel/{identifier(job_id)}", {})

    def _request(self, method: str, operation: str, payload: dict | None = None) -> dict:
        key = os.environ.get("RUNPOD_API_KEY")

        if not key:
            raise RuntimeError("Set RUNPOD_API_KEY in the launching shell")

        request = Request(
            f"https://api.runpod.ai/v2/{self.endpoint_id}/{operation}",
            data=json.dumps(payload, allow_nan=False).encode() if payload is not None else None,
            headers={
                "Authorization": f"Bearer {key}",
                "Content-Type": "application/json",
                "User-Agent": "monki-rt/0.1.0",
            },
            method=method,
        )

        try:
            with build_opener(NoRedirect()).open(request, timeout=30) as response:
                # The worker bounds inline artifacts; reject unexpectedly large replies.
                raw = response.read(12 * 1024 * 1024 + 1)

                if len(raw) > 12 * 1024 * 1024:
                    raise RuntimeError("Runpod response exceeds the inline result limit")

                result = json.loads(raw)
        except HTTPError as error:
            raise RuntimeError(f"Runpod HTTP {error.code}; request was not retried") from None
        except URLError, TimeoutError, OSError:
            raise RuntimeError(
                "Runpod connection failed; submission may have succeeded. "
                "Inspect endpoint jobs before submitting again."
            ) from None

        if not isinstance(result, dict):
            raise RuntimeError("Runpod returned a non-object response")

        return result


def identifier(value: str) -> str:
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9_-]{1,128}", value):
        raise ValueError("Invalid Runpod endpoint/job identifier")

    return value
