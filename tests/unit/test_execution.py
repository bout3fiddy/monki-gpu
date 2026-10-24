"""Enforce the strategy interface and prevent pretend execution or CPU fallback."""

import pytest
from monki.execution import CPU, CUDA, Backend, Metal


def test_backend_is_abstract():
    with pytest.raises(TypeError, match="abstract"):
        Backend()


@pytest.mark.parametrize("missing_method", ["prepare", "run_batch", "close"])
def test_every_lifecycle_method_is_required(missing_method):
    methods = {
        "prepare": lambda self, scene: None,
        "run_batch": lambda self, controls: object(),
        "close": lambda self: None,
    }
    del methods[missing_method]
    incomplete_backend = type("IncompleteBackend", (Backend,), methods)

    with pytest.raises(TypeError, match=missing_method):
        incomplete_backend()


def test_context_manager_closes_on_success_and_failure():
    class RecordingBackend(Backend):
        """Exercise dispatch without pretending to implement scientific transport."""

        def __init__(self):
            self.events = []

        def prepare(self, scene):
            self.events.append(("prepare", scene))

        def run_batch(self, controls):
            self.events.append(("batch", controls))
            return controls

        def close(self):
            self.events.append(("close", None))

    scene = object()
    controls = object()
    backend = RecordingBackend()

    with backend as strategy:
        assert strategy is backend
        strategy.prepare(scene)
        assert strategy.run_batch(controls) is controls

    assert backend.events == [("prepare", scene), ("batch", controls), ("close", None)]

    with pytest.raises(RuntimeError, match="simulation failed"), backend:
        raise RuntimeError("simulation failed")

    assert backend.events[-1] == ("close", None)
    assert backend.events.count(("close", None)) == 2


@pytest.mark.parametrize("backend_type", [CPU, Metal, CUDA])
def test_cpu_and_metal_do_not_accept_cloud_settings(backend_type):
    with pytest.raises(TypeError):
        backend_type(provider="aws")


@pytest.mark.parametrize(
    "backend_type,settings",
    [
        (CPU, {"threads": 0}),
        (CPU, {"threads": True}),
        (Metal, {"device": -1}),
        (Metal, {"device": True}),
    ],
)
def test_invalid_settings_fail_before_any_execution(backend_type, settings):
    with pytest.raises(ValueError):
        backend_type(**settings)


@pytest.mark.parametrize(
    "settings",
    [
        {"seed": -1},
        {"seed": True},
        {"seed": 2**32},
        {"photons_per_dispatch": 0},
        {"photons_per_dispatch": 65537},
        {"tally_shards": 0},
        {"tally_shards": 65},
        {"tally_shards": True},
        {"profile": "yes"},
        {"max_steps": 0},
        {"max_steps": 100001},
        {"save_path": "yes"},
    ],
)
@pytest.mark.parametrize("backend_type", [Metal, CUDA])
def test_invalid_gpu_settings_fail_before_gpu_access(settings, backend_type):
    with pytest.raises(ValueError):
        backend_type(**settings)
