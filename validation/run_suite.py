"""Run CPU reference cases and compare against separately supplied MONKI results."""

import argparse
import json
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path

from monki import ReferenceCase, ReferenceRunner, SimulationResult, simulate
from monki.execution import CPU, CUDA, Metal
from monki.validation import check_result, compare, select_views
from path_statistics import PathRecorder, compare_paths

ROOT = Path(__file__).resolve().parents[1]
UPSTREAM = ROOT / "fortran/monki/upstream"
GOLDENS = {
    "clear": "MONKI_MLS_cloudfree_340_340_test_paper_as_005_th0_045_ph0_000_ph_000.nc",
    "layered-cloud": "MONKI_MLS_1layercloud_340_340_test_paper_as_005_th0_045_ph0_000_ph_000.nc",
    "box-clouds": "MONKI_MLS_2boxclouds_340_340_test_paper_as_005_th0_045_ph0_000_ph_000.nc",
}


def run_suite(
    destination: Path, selected: list[str], *, metal: bool = False, cuda: bool = False
) -> bool:
    destination.mkdir(parents=True, exist_ok=False)
    report = {
        "schema_version": 1,
        "scope": "CPU vs original/archived outputs; optional GPU vs CPU",
        "metal_requested": metal,
        "cuda_requested": cuda,
        "statistical_policy": {
            "family_error_rate": 1e-6,
            "absolute_tolerance": 1e-10,
            "coverage": "archived radiance; Metal radiance and round-based layer AMF ratios",
        },
        "cases": [],
        "passed": True,
    }

    for name in selected:
        case = ReferenceCase(**json.loads((ROOT / f"validation/cases/{name}.json").read_text()))
        profile = "radiance" if name == "box-clouds" else "paths"
        runner = ReferenceRunner(ROOT / f"build/reference/{profile}/monki", UPSTREAM)
        entry = {"case": asdict(case)}

        try:
            actual = runner.run(case, destination / name)
            golden = select_views(
                SimulationResult.load(UPSTREAM / "OUTPUT/testrun" / GOLDENS[name]),
                tuple(case.view_zeniths),
            )
            cpu_rounds = PathRecorder(CPU(save_path=name != "box-clouds"))
            candidate = simulate(case, cpu_rounds, destination / f"{name}-cpu.nc")
            deterministic = compare(actual, candidate, mode="deterministic")
            comparison = compare(golden, candidate)
            entry.update(
                checks=check_result(candidate),
                original_comparison=deterministic.to_dict(),
                comparison=comparison.to_dict(),
                passed=comparison.passed and deterministic.passed,
            )

            for label, backend_type in (("metal", Metal), ("cuda", CUDA)):
                if not (metal if label == "metal" else cuda):
                    continue

                metal_rounds = PathRecorder(backend_type(save_path=name != "box-clouds"))
                gpu = simulate(case, metal_rounds, destination / f"{name}-{label}.nc")
                gpu_comparison = compare(candidate, gpu)
                path_comparison = compare_paths(candidate, gpu, cpu_rounds, metal_rounds)
                entry[f"{label}_comparison"] = gpu_comparison.to_dict()
                entry[f"{label}_path_comparison"] = path_comparison
                entry[f"{label}_device"] = metal_rounds.provenance["device"]
                entry["passed"] = (
                    entry["passed"] and gpu_comparison.passed and path_comparison["passed"]
                )
        except Exception as exc:
            entry.update(passed=False, error=str(exc))

        report["cases"].append(entry)
        report["passed"] = report["passed"] and entry["passed"]
        print(f"{name}: {'PASS' if entry['passed'] else 'FAIL'}", flush=True)
        (destination / "report.json").write_text(json.dumps(report, indent=2) + "\n")

    return report["passed"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", choices=sorted(GOLDENS))
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--metal", action="store_true", help="Also compare real Metal GPU transport"
    )
    parser.add_argument("--cuda", action="store_true", help="Require real CUDA transport")
    args = parser.parse_args()
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    destination = args.output or ROOT / "validation/runs" / stamp
    success = run_suite(destination, args.case or list(GOLDENS), metal=args.metal, cuda=args.cuda)
    print(destination / "report.json")
    raise SystemExit(0 if success else 1)
