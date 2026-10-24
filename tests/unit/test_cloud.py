"""Cloud lifecycle, request isolation and artifact integrity without paid jobs."""

import base64
import hashlib
import json
from dataclasses import replace

import pytest
from monki.cloud import JobState, JobStatus, Provider, RunpodServerless
from monki.cloud.jobs import save_artifacts, simulation_request, validate_request
from monki.cloud.worker import handler
from monki.environment import native_environment
from monki.input import ReferenceCase

CASE = ReferenceCase("clear", "testrun/MONKI_MLS_cloudfree_340_340_test.in")


def test_provider_requires_all_lifecycle_operations():
    with pytest.raises(TypeError, match="abstract"):
        Provider()


class RecordedProvider(Provider):
    def __init__(self, statuses):
        self.statuses = iter(statuses)
        self.submissions = 0

    def submit(self, request):
        self.submissions += 1
        return "job"

    def status(self, job_id):
        return next(self.statuses)

    def cancel(self, job_id):
        pass


def test_wait_observes_existing_job_without_resubmission(monkeypatch):
    monkeypatch.setattr("monki.cloud.provider.time.sleep", lambda seconds: None)
    provider = RecordedProvider(
        [
            JobStatus("job", JobState.QUEUED),
            JobStatus("job", JobState.RUNNING),
            JobStatus("job", JobState.COMPLETED, {"version": 1}),
        ]
    )
    assert provider.wait("job") == {"version": 1}
    assert provider.submissions == 0


@pytest.mark.parametrize("state", [JobState.FAILED, JobState.CANCELLED, JobState.TIMED_OUT])
def test_terminal_failure_is_explicit(state):
    provider = RecordedProvider([JobStatus("job", state)])

    with pytest.raises(RuntimeError, match=state.value):
        provider.wait("job")


def test_local_timeout_preserves_remote_job(monkeypatch):
    ticks = iter([0, 2])
    monkeypatch.setattr("monki.cloud.provider.time.monotonic", lambda: next(ticks))
    provider = RecordedProvider([JobStatus("job", JobState.RUNNING)])

    with pytest.raises(TimeoutError, match="may still be running"):
        provider.wait("job", timeout=1)


def test_runpod_maps_requests_and_states(monkeypatch):
    provider = RunpodServerless("endpoint", execution_timeout=300, ttl=900)
    calls = []

    def request(method, operation, payload=None):
        calls.append((method, operation, payload))
        return {"id": "job", "status": "IN_PROGRESS"}

    monkeypatch.setattr(provider, "_request", request)
    payload = simulation_request(CASE)
    assert provider.submit(payload) == "job"
    assert calls[0] == (
        "POST",
        "run",
        {"input": payload, "policy": {"executionTimeout": 300000, "ttl": 900000}},
    )
    assert provider.status("job").state == JobState.RUNNING
    provider.cancel("job")
    assert calls[-1][:2] == ("POST", "cancel/job")


def test_credentials_stay_out_of_worker_and_job(monkeypatch):
    monkeypatch.setenv("RUNPOD_API_KEY", "private-test-token")
    assert "RUNPOD_API_KEY" not in native_environment()
    assert "private-test-token" not in json.dumps(simulation_request(CASE))
    assert "private-test-token" not in repr(RunpodServerless("endpoint"))


def test_missing_credential_fails_before_network(monkeypatch):
    monkeypatch.delenv("RUNPOD_API_KEY", raising=False)

    with pytest.raises(RuntimeError, match="Set RUNPOD_API_KEY"):
        RunpodServerless("endpoint").submit(simulation_request(CASE))


def test_ambiguous_submission_is_not_retried_or_leaked(monkeypatch):
    monkeypatch.setenv("RUNPOD_API_KEY", "private-test-token")
    calls = []

    class Offline:
        def open(self, request, timeout):
            calls.append(request)
            raise TimeoutError("private-test-token")

    monkeypatch.setattr("monki.cloud.runpod.build_opener", lambda *args: Offline())

    with pytest.raises(RuntimeError, match="may have succeeded") as failure:
        RunpodServerless("endpoint").submit(simulation_request(CASE))

    assert len(calls) == 1
    assert calls[0].get_header("User-agent") == "monki-rt/0.1.0"
    assert "private-test-token" not in str(failure.value)


@pytest.mark.parametrize("identifier", ["../other", "endpoint?secret", "https://host", ""])
def test_runpod_identifier_cannot_redirect_credentials(identifier):
    with pytest.raises(ValueError, match="identifier"):
        RunpodServerless(identifier)


def test_request_rejects_arbitrary_files_and_code():
    with pytest.raises(ValueError, match="supported forward"):
        simulation_request(replace(CASE, input_file="unknown.in"))

    request = simulation_request(CASE)
    request["command"] = "arbitrary command"

    with pytest.raises(ValueError, match="schema"):
        validate_request(request)


def test_worker_delegates_complete_portable_request(monkeypatch):
    request = simulation_request(CASE)
    received = []
    monkeypatch.setattr("monki.cloud.worker.execute", lambda value: received.append(value) or {})
    assert handler({"id": "job", "input": request}) == {}
    assert received == [request]


def result_bundle():
    files = {"result.nc": b"netcdf-test-bytes", "result.run.json": b'{"backend":"CUDA"}'}
    return {
        "version": 1,
        "artifacts": {
            name: {
                "base64": base64.b64encode(data).decode(),
                "sha256": hashlib.sha256(data).hexdigest(),
            }
            for name, data in files.items()
        },
    }


def test_artifacts_validate_before_writing_and_never_overwrite(tmp_path):
    result = result_bundle()
    output = tmp_path / "run"
    save_artifacts(result, output)
    assert (output / "result.nc").read_bytes() == b"netcdf-test-bytes"

    with pytest.raises(FileExistsError):
        save_artifacts(result, output)

    result["artifacts"]["result.nc"]["sha256"] = "wrong"

    with pytest.raises(ValueError, match="hash"):
        save_artifacts(result, tmp_path / "corrupt")

    assert not (tmp_path / "corrupt").exists()


def test_artifact_paths_are_fixed(tmp_path):
    result = result_bundle()
    result["artifacts"]["../outside"] = result["artifacts"].pop("result.nc")

    with pytest.raises(ValueError, match="schema"):
        save_artifacts(result, tmp_path / "run")

    assert not (tmp_path / "run").exists()


def test_remote_validation_uses_fixed_command_and_clean_environment(monkeypatch):
    calls = []
    monkeypatch.setenv("RUNPOD_API_KEY", "private-test-token")

    class Process:
        returncode = 1

        def __init__(self, command, **kwargs):
            calls.append((command, kwargs))
            kwargs["stdout"].write(b"failed invariant")

        def wait(self, timeout=None):
            return self.returncode

        def poll(self):
            return self.returncode

    monkeypatch.setattr("monki.cloud.worker.subprocess.Popen", Process)
    output = handler({"input": {"version": 1, "operation": "validate"}})
    assert output["validation"]["passed"] is False
    assert output["validation"]["log"] == "failed invariant"
    assert calls[0][0][1:] == ["scripts/verify.py", "--cuda"]
    assert "RUNPOD_API_KEY" not in calls[0][1]["env"]
    assert calls[0][1]["start_new_session"] is True


def test_remote_validation_timeout_reaps_process_group(monkeypatch):
    import subprocess

    killed = []

    class Process:
        pid = 12345
        returncode = None

        def __init__(self, *args, **kwargs):
            pass

        def wait(self, timeout=None):
            if timeout:
                raise subprocess.TimeoutExpired("verify", timeout)

            self.returncode = -9
            return self.returncode

        def poll(self):
            return self.returncode

    monkeypatch.setattr("monki.cloud.worker.subprocess.Popen", Process)
    monkeypatch.setattr("monki.cloud.worker.os.killpg", lambda pid, sig: killed.append(pid))
    result = handler({"input": {"version": 1, "operation": "validate"}})
    assert result["validation"]["timed_out"] is True
    assert result["validation"]["passed"] is False
    assert killed == [12345]
