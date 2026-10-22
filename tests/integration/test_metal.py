"""Real-device checks: no CPU transport fallback, physics helpers, and lifecycle."""

import json
import platform
import subprocess
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from monki import ReferenceCase, simulate
from monki.execution import CPU, BatchControls, Metal
from monki.validation import check_result, compare

ROOT = Path(__file__).resolve().parents[2]
CASE = ReferenceCase("clear", "testrun/MONKI_MLS_cloudfree_340_340_test.in")
pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(platform.system() != "Darwin", reason="Metal requires macOS"),
]


def probe_kernel(tmp_path, body):
    """Dispatch production helpers in isolation, without reimplementing their equations."""

    source = (ROOT / "kernels/metal/transport.metal").read_text()
    prefix = source[: source.index("kernel void monte_carlo")]
    path = tmp_path / "probe.metal"
    path.write_text(
        prefix
        + """
kernel void monte_carlo(device const float *data [[buffer(0)]],
                       device float *output [[buffer(1)]],
                       uint lane [[thread_position_in_grid]]) {
    if (lane != 0) return;
"""
        + body
        + "\n}\n"
    )
    return path


def test_gpu_philox_matches_published_known_answer(tmp_path):
    # Random123 tests/kat_vectors, Philox4x32-10 zero counter and key.
    # https://github.com/DEShawResearch/random123/blob/main/tests/kat_vectors
    kernel = probe_kernel(
        tmp_path,
        """
    uint4 words = philox(uint4(0), uint2(0));
    for (uint i = 0; i < 4; ++i) output[i] = as_type<float>(words[i]);
    """,
    )

    with Metal(kernel=kernel, tally_shards=1) as backend:
        backend.prepare(CASE)
        batch = backend.run_batch(BatchControls(1, 1))

    actual = batch.stokes.ravel(order="F")[:4].astype(np.float32).view(np.uint32)
    np.testing.assert_array_equal(actual, [0x6627E8D5, 0xE169C58D, 0xBC57AC4C, 0x9B00DBD8])


def test_polarization_helpers_match_original_fortran(tmp_path):
    directory = ROOT / "build/reference/paths"
    receipt = json.loads((directory / "build.json").read_text())
    settings = receipt["identity"]
    executable = tmp_path / "probe"
    objects = [str(path) for path in directory.glob("*.o") if path.name != "main.o"]
    command = [
        receipt["commands"][0][0],
        *settings["flags"],
        "-I",
        str(directory),
        str(ROOT / "tests/native/metal_probe.f90"),
        *objects,
        *settings["libraries"],
        "-o",
        str(executable),
    ]
    subprocess.run(command, check=True, capture_output=True)
    expected = np.fromstring(subprocess.check_output([str(executable)], text=True), sep=" ")
    kernel = probe_kernel(
        tmp_path,
        """
    Scene s = scene_tables(data);
    s.depol = 0.03f;
    float uo[3] = {-0.7f, -1.0f, -0.4f};
    float u1[3] = {0.4f, 0.6f, 1.0f};
    float phio[3] = {0.3f, 0.2f, 0.9f};
    float phi1[3] = {1.8f, 0.8f, 0.4f};
    for (int i = 0; i < 3; ++i) {
        float cosine = uo[i]*u1[i] + sqrt(1-uo[i]*uo[i])*sqrt(1-u1[i]*u1[i])*cos(phio[i]-phi1[i]);
        float3 result = scatter_stokes(s, true, cosine, uo[i], u1[i],
            phio[i], phi1[i], float3(1, 0.2f, -0.1f));
        for (int j = 0; j < 3; ++j) output[3*i+j] = result[j];
    }
    """,
    )

    with Metal(kernel=kernel, tally_shards=1) as backend:
        backend.prepare(CASE)
        result = backend.run_batch(BatchControls(1, 1))

    np.testing.assert_allclose(result.stokes.ravel(order="F"), expected, rtol=1e-5, atol=1e-7)


def test_metal_never_uses_cpu_photon_transport(monkeypatch):
    def forbid_cpu_batch(*args):
        raise AssertionError("Metal attempted CPU photon transport")

    monkeypatch.setattr(CPU, "run_batch", forbid_cpu_batch)

    with Metal() as backend:
        backend.prepare(CASE)
        first = backend.run_batch(BatchControls(1, 200))
        snapshot = first.stokes.copy()
        second = backend.run_batch(BatchControls(2, 200))
        assert backend.provenance["gpu_dispatches"] == 2
        assert backend.provenance["device"]
        assert np.all(first.stokes[0] > 0)
        assert not np.array_equal(first.stokes, second.stokes)
        backend.prepare(CASE)
        replay = backend.run_batch(BatchControls(1, 200))
        # Atomic accumulation order may differ; random histories must not.
        np.testing.assert_allclose(first.stokes, replay.stokes, rtol=2e-5, atol=2e-5)

    np.testing.assert_array_equal(first.stokes, snapshot)
    backend.close()

    with pytest.raises(RuntimeError, match="prepared"):
        backend.run_batch(BatchControls(1, 200))


def test_chunking_preserves_photon_streams():
    with Metal(photons_per_dispatch=32) as small, Metal(photons_per_dispatch=200) as large:
        small.prepare(CASE)
        large.prepare(CASE)
        a = small.run_batch(BatchControls(1, 200))
        b = large.run_batch(BatchControls(1, 200))
        np.testing.assert_allclose(a.stokes, b.stokes, rtol=2e-5, atol=2e-5)
        np.testing.assert_allclose(a.weighted_paths, b.weighted_paths, rtol=2e-5, atol=2e-5)


def test_off_plane_polarization_matches_cpu(tmp_path):
    case = replace(
        CASE, photons_per_column=2000, solar_azimuth=37, view_azimuth=109, view_zeniths=(20, 45, 70)
    )
    cpu = simulate(case, CPU(), tmp_path / "cpu.nc")
    metal = simulate(case, Metal(seed=881), tmp_path / "metal.nc")
    result = compare(cpu, metal)
    assert result.passed, result.details
    assert check_result(metal) == []
    assert np.max(np.abs(metal.radiance[..., 2, :])) > 1e-3


def test_truncated_histories_fail_and_release_session():
    with Metal(max_steps=1) as backend:
        backend.prepare(CASE)

        with pytest.raises(RuntimeError, match="round rejected"):
            backend.run_batch(BatchControls(1, 200))

        assert backend._handle is None
        assert backend.metadata is None


def test_invalid_device_does_not_fall_back():
    with Metal(device=2**32 - 1) as backend:
        with pytest.raises(RuntimeError, match="unavailable"):
            backend.prepare(CASE)

        assert backend._handle is None


def test_batch_validation_precedes_dispatch():
    with Metal() as backend:
        backend.prepare(CASE)

        with pytest.raises(ValueError, match="sequential"):
            backend.run_batch(BatchControls(2, 200))

        assert backend.provenance["gpu_dispatches"] == 0
        assert backend.run_batch(BatchControls(1, 200)).controls.round_index == 1


def test_backward_mode_is_explicitly_unsupported():
    with Metal() as backend, pytest.raises(NotImplementedError, match="backward"):
        backend.prepare(replace(CASE, backward=True))


@pytest.mark.parametrize("field,value", [(0, 2), (4, 4), (6, 2), (7, 0), (9, 0)])
def test_unsupported_prepared_tables_are_rejected(field, value):
    from monki.execution.metal_tables import pack_tables

    with CPU() as preparer:
        preparer.prepare(CASE)
        raw = preparer.transport_tables()

    raw[field] = value

    with pytest.raises((NotImplementedError, ValueError)):
        pack_tables(raw)


def test_scalar_transport_matches_cpu(tmp_path):
    # A private input variant exercises nmat=1 without changing the vendored fixture.
    upstream = tmp_path / "upstream"
    source = ROOT / "fortran/monki/upstream"
    fixture = upstream / "INPUT" / CASE.input_file
    fixture.parent.mkdir(parents=True)
    (upstream / "DATA").symlink_to(source / "DATA")
    lines = (source / "INPUT" / CASE.input_file).read_text().splitlines()
    lines[24] = "1 :size of Stokes vector (1, 3, or 4)"
    fixture.write_text("\n".join(lines) + "\n")
    case = replace(CASE, photons_per_column=2000)
    cpu = simulate(case, CPU(upstream=upstream), tmp_path / "scalar-cpu.nc")
    gpu = simulate(case, Metal(upstream=upstream), tmp_path / "scalar-metal.nc")
    result = compare(cpu, gpu)
    assert result.passed, result.details
    assert gpu.radiance.shape[-2] == 1


@pytest.mark.parametrize(
    "source,message",
    [
        ("invalid shader syntax", "error"),
        (
            "#include <metal_stdlib>\nusing namespace metal;\nkernel void other() {}",
            "missing monte_carlo",
        ),
    ],
)
def test_invalid_kernel_fails_without_leaving_a_session(tmp_path, source, message):
    kernel = tmp_path / "invalid.metal"
    kernel.write_text(source)

    with Metal(kernel=kernel, tally_shards=1) as backend:
        with pytest.raises(RuntimeError, match=message):
            backend.prepare(CASE)

        assert backend._handle is None
        assert backend.metadata is None


@pytest.mark.parametrize(
    "setting,value",
    [
        ("seed", 1),
        ("save_path", False),
        ("photons_per_dispatch", 32),
        ("tally_shards", 1),
        ("profile", True),
    ],
)
def test_configuration_changes_require_reprepare(setting, value):
    with Metal() as backend:
        backend.prepare(CASE)
        setattr(backend, setting, value)

        with pytest.raises(ValueError, match="reprepare"):
            backend.run_batch(BatchControls(1, 200))

        assert backend.provenance["gpu_dispatches"] == 0


@pytest.mark.parametrize("save_path", [False, True])
def test_shard_reduction_preserves_all_estimators(save_path):
    case = CASE

    if not save_path:
        case = ReferenceCase(**json.loads((ROOT / "validation/cases/box-clouds.json").read_text()))

    with (
        Metal(tally_shards=1, save_path=save_path) as single,
        Metal(tally_shards=32, save_path=save_path, profile=True) as sharded,
    ):
        single.prepare(case)
        sharded.prepare(case)
        a = single.run_batch(BatchControls(1, 200))
        b = sharded.run_batch(BatchControls(1, 200))
        np.testing.assert_allclose(a.stokes, b.stokes, rtol=2e-5, atol=2e-5)

        if save_path:
            np.testing.assert_allclose(a.weighted_paths, b.weighted_paths, rtol=2e-5, atol=2e-5)
            np.testing.assert_allclose(a.weights, b.weights, rtol=2e-5, atol=2e-5)

        assert sharded.provenance["gpu_seconds"] > 0
        assert "gpu_seconds" not in single.provenance
        sharded.prepare(case)
        assert "gpu_seconds" not in sharded.provenance


def test_large_dispatch_preserves_3d_photon_streams_and_tail():
    case = ReferenceCase(**json.loads((ROOT / "validation/cases/box-clouds.json").read_text()))

    with (
        Metal(save_path=False, photons_per_dispatch=4096, tally_shards=1) as small,
        Metal(save_path=False) as large,
    ):
        small.prepare(case)
        large.prepare(case)
        a = small.run_batch(BatchControls(1, 200))
        b = large.run_batch(BatchControls(1, 200))
        # Same 120,000 histories, including a partial last dispatch. Only float32
        # tally summation order differs; uncertainty tolerances are not used here.
        np.testing.assert_allclose(a.stokes, b.stokes, rtol=2e-5, atol=2e-5)
        assert small.provenance["gpu_dispatches"] == 30
        assert large.provenance["gpu_dispatches"] == 2


def test_grouped_rounds_preserve_individual_streams_and_partial_group():
    requests = [BatchControls(index, 200) for index in range(1, 36)]

    with Metal() as grouped, Metal() as sequential:
        grouped.prepare(CASE)
        sequential.prepare(CASE)
        actual = grouped.run_batches(requests)
        expected = [sequential.run_batch(request) for request in requests]
        assert grouped.provenance["gpu_dispatches"] == 2
        assert sequential.provenance["gpu_dispatches"] == 35

        for a, b in zip(actual, expected, strict=True):
            assert a.controls == b.controls
            np.testing.assert_allclose(a.stokes, b.stokes, rtol=2e-5, atol=2e-5)
            np.testing.assert_allclose(a.weighted_paths, b.weighted_paths, rtol=2e-5, atol=2e-5)
            np.testing.assert_allclose(a.weights, b.weights, rtol=2e-5, atol=2e-5)

        saved = actual[0].stokes.copy()
        grouped.run_batch(BatchControls(36, 200))
        np.testing.assert_array_equal(actual[0].stokes, saved)


def test_grouped_rounds_validate_all_requests_before_dispatch():
    with Metal() as backend:
        backend.prepare(CASE)

        with pytest.raises(ValueError, match="sequential"):
            backend.run_batches([BatchControls(1, 200), BatchControls(3, 200)])

        assert backend.provenance["gpu_dispatches"] == 0
        assert backend.run_batches([]) == []
        requests = [BatchControls(1, 200), BatchControls(2, 100)]
        assert [batch.controls for batch in backend.run_batches(requests)] == requests
        assert backend.provenance["gpu_dispatches"] == 2


def test_grouped_failure_rejects_all_rounds_and_closes():
    with Metal(max_steps=1) as backend:
        backend.prepare(CASE)

        with pytest.raises(RuntimeError, match="round rejected"):
            backend.run_batches([BatchControls(1, 200), BatchControls(2, 200)])

        assert backend._handle is None


@pytest.mark.parametrize("name", ["clear", "box-clouds"])
def test_uniform_tail_matches_explicit_detector_traversal(tmp_path, name):
    # Compare the shortcut against the retained cell-by-cell detector on the GPU.
    # Off-plane directions exercise both periodic horizontal coordinates.
    case = ReferenceCase(**json.loads((ROOT / f"validation/cases/{name}.json").read_text()))
    case = replace(case, view_azimuth=109, view_zeniths=(20, 45, 70))
    fast = probe_kernel(
        tmp_path,
        """
    Scene s = scene_tables(data);
    Photon p;
    p.i = 0; p.j = 0; p.k = 12;
    p.direction = float3(0.3f, 0.4f, sqrt(0.75f));
    p.phi = atan2(p.direction.y, p.direction.x);
    p.remaining = float3(0.37f*s.xunit, 0.61f*s.yunit, 0.43f*s.zunit[p.k]);
    float history[MAX_LAYERS];
    for (uint k = 0; k < MAX_LAYERS; ++k) history[k] = 0.3f;
    for (int angle = 0; angle < s.angles; ++angle)
        detector(s, p, angle, float3(1, 0.2f, -0.1f), history,
                 (device atomic_float *)output, 100000);
    """,
    )
    slow = tmp_path / "explicit.metal"
    slow.write_text(fast.read_text().replace("if (ray.k >= MONKI_UNIFORM_TOP)", "if (false)"))

    with (
        Metal(kernel=fast, tally_shards=1, save_path=name != "box-clouds") as optimized,
        Metal(kernel=slow, tally_shards=1, save_path=name != "box-clouds") as explicit,
    ):
        optimized.prepare(case)
        explicit.prepare(case)
        a = optimized.run_batch(BatchControls(1, 1))
        b = explicit.run_batch(BatchControls(1, 1))
        np.testing.assert_allclose(a.stokes, b.stokes, rtol=1e-4, atol=1e-6)

        if a.weighted_paths is not None:
            np.testing.assert_allclose(a.weighted_paths, b.weighted_paths, rtol=1e-4, atol=1e-6)
            np.testing.assert_allclose(a.weights, b.weights, rtol=1e-4, atol=1e-6)


def test_uniform_tail_never_averages_a_varying_layer():
    from monki.execution.metal_tables import pack_tables

    case = ReferenceCase(**json.loads((ROOT / "validation/cases/box-clouds.json").read_text()))

    with CPU(save_path=False) as preparer:
        preparer.prepare(case)
        raw = preparer.transport_tables()

    packed, (_, nx, ny, nz, angles) = pack_tables(raw)
    assert packed[15] == 13
    beta_start = 18 + nz + angles + 3 * nx * ny
    # Move the highest nonuniform layer to the last physical layer. There must
    # now be no physical layer where the uniform-tail shortcut can run.
    raw[beta_start + (nz - 2) * nx * ny] += 0.01
    changed, _ = pack_tables(raw)
    assert changed[15] == nz - 1
    np.testing.assert_array_equal(changed[-nz:][-2:], [0, 0])


def test_stale_bridge_is_rejected_before_native_call(monkeypatch):
    import ctypes

    monkeypatch.setattr(ctypes, "CDLL", lambda path: object())

    with Metal() as backend:
        with pytest.raises(RuntimeError, match="Rebuild the Metal bridge"):
            backend.prepare(CASE)

        assert backend._handle is None
