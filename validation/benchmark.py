"""Compare fixed-work CPU and Metal performance with statistical output checks."""

import argparse
import json
import os
import platform
from dataclasses import asdict, replace
from datetime import UTC, datetime
from pathlib import Path
from statistics import median
from time import perf_counter

from monki import ReferenceCase, simulate
from monki.execution import CPU, Backend, Metal, RawTallies
from monki.validation import compare
from path_statistics import PathRecorder, compare_paths

ROOT = Path(__file__).resolve().parents[1]
CASES = ("clear", "layered-cloud", "box-clouds")


class TimedBackend(Backend):
    """Measure synchronous backend boundaries without changing transport code."""

    def __init__(self, backend: Backend):
        self.backend = backend
        self.prepare_seconds = 0.0
        self.batch_seconds: list[float] = []
        self.close_seconds = 0.0

    def prepare(self, scene: object) -> None:
        start = perf_counter()
        self.backend.prepare(scene)
        self.prepare_seconds = perf_counter() - start
        self.metadata = self.backend.metadata
        self.provenance = self.backend.provenance

    def run_batch(self, controls: object) -> RawTallies:
        start = perf_counter()
        batch = self.backend.run_batch(controls)
        self.batch_seconds.append(perf_counter() - start)

        if not isinstance(batch, RawTallies):
            raise TypeError("Benchmark requires raw photon tallies")

        return batch

    def run_batches(self, controls: list) -> list:
        # Time the complete synchronous operation, including any grouped dispatches.
        # Scientific checks remain outside this transport boundary.
        start = perf_counter()
        batches = self.backend.run_batches(controls)
        self.batch_seconds.append(perf_counter() - start)
        return batches

    def close(self) -> None:
        start = perf_counter()
        self.backend.close()
        self.close_seconds += perf_counter() - start


def run_benchmark(args: argparse.Namespace) -> Path:
    """Run fresh sessions, alternating device order, and preserve every measured output."""

    if platform.system() != "Darwin":
        raise ValueError("This comparison requires a Mac with Metal; there is no CPU fallback")

    cases = []

    for name in args.case or CASES:
        case = ReferenceCase(**json.loads((ROOT / f"validation/cases/{name}.json").read_text()))
        case = replace(
            case,
            photons_per_column=args.photons_per_column or case.photons_per_column,
            rounds=args.rounds,
        )

        if any(case.photons_per_column % threads for threads in args.threads):
            raise ValueError(
                f"{name}: photons per column must be divisible by every CPU thread count"
            )

        cases.append(case)

    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    destination = args.output or ROOT / "validation/runs" / f"benchmark-{stamp}"
    destination.mkdir(parents=True, exist_ok=False)
    report = {
        "schema_version": 1,
        "host": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "logical_cpus": os.cpu_count(),
            "python": platform.python_version(),
        },
        "policy": {
            "repeats": args.repeats,
            "warmup_sessions_per_backend": 1,
            "summary": "median of fresh-session measurements; raw samples retained",
            "transport": "synchronous run_batches including copies, dispatch and owned results",
            "end_to_end": "simulate including preparation, statistics, NetCDF, receipt and cleanup",
            "validation": "radiance and enabled layer AMFs; outside timed simulate calls",
            "precision": "CPU float64 versus Metal float32 transport; different random streams",
            "scope": "fixed photon count; compilation excluded",
        },
        "cases": [],
        "passed": True,
    }

    for case in cases:
        paths = case.name != "box-clouds"
        factories = {
            f"cpu-{n}": lambda n=n, paths=paths: CPU(threads=n, save_path=paths)
            for n in args.threads
        }
        factories["metal"] = lambda paths=paths: Metal(
            save_path=paths,
            photons_per_dispatch=args.metal_dispatch_size,
            tally_shards=args.metal_tally_shards,
            profile=args.profile_metal,
        )
        samples = {name: [] for name in factories}
        entry = {"case": asdict(case), "save_path": paths, "samples": samples, "comparisons": []}
        report["cases"].append(entry)

        # Warm the runtime, filesystem cache, and shader cache with a separate session.
        # Measured sessions still pay their own preparation and device allocation costs.
        for name, factory in factories.items():
            print(f"{case.name}: warmup {name}", flush=True)
            simulate(case, factory(), destination / f"{case.name}-{name}-warmup.nc")

        for repeat in range(args.repeats):
            order = list(factories)

            if repeat % 2:
                order.reverse()

            results = {}
            recorders = {}

            for name in order:
                timed = TimedBackend(factories[name]())
                recorder = PathRecorder(timed)
                output = destination / f"{case.name}-{name}-{repeat + 1}.nc"
                start = perf_counter()
                result = simulate(case, recorder, output)
                elapsed = perf_counter() - start
                receipt = json.loads(output.with_suffix(".run.json").read_text())
                transport = sum(timed.batch_seconds)
                sample = {
                    "repeat": repeat + 1,
                    "prepare_seconds": timed.prepare_seconds,
                    "transport_seconds": transport,
                    "close_seconds": timed.close_seconds,
                    "end_to_end_seconds": elapsed,
                    "other_seconds": elapsed
                    - timed.prepare_seconds
                    - transport
                    - timed.close_seconds,
                    "batch_seconds": timed.batch_seconds,
                    "total_photons": receipt["total_photons"],
                    "photons_per_second": receipt["total_photons"] / transport,
                    "provenance": receipt,
                    "output": str(output),
                }
                samples[name].append(sample)
                results[name] = result
                recorders[name] = recorder
                print(f"{case.name}: {name} repeat {repeat + 1}: {elapsed:.3f}s total", flush=True)

            # Validate every measured result before presenting a speed ratio as usable.
            for name in factories:
                if name == "metal":
                    continue

                radiance = compare(results[name], results["metal"])
                amf = compare_paths(
                    results[name], results["metal"], recorders[name], recorders["metal"]
                )
                passed = radiance.passed and amf["passed"]
                entry["comparisons"].append(
                    {
                        "cpu": name,
                        "repeat": repeat + 1,
                        "radiance": radiance.to_dict(),
                        "paths": amf,
                        "passed": passed,
                    }
                )
                report["passed"] = report["passed"] and passed

        summary = {
            name: {
                field: median(sample[field] for sample in runs)
                for field in (
                    "prepare_seconds",
                    "transport_seconds",
                    "end_to_end_seconds",
                    "photons_per_second",
                )
            }
            for name, runs in samples.items()
        }
        entry["medians"] = summary
        entry["speedups"] = {
            name: {
                "transport": values["transport_seconds"] / summary["metal"]["transport_seconds"],
                "end_to_end": values["end_to_end_seconds"] / summary["metal"]["end_to_end_seconds"],
            }
            for name, values in summary.items()
            if name != "metal"
        }
        (destination / "report.json").write_text(json.dumps(report, indent=2) + "\n")

    print("\nCase           CPU threads   Transport speedup   End-to-end speedup")

    for entry in report["cases"]:
        for name, speedup in entry["speedups"].items():
            print(
                f"{entry['case']['name']:14} {name[4:]:>11} "
                f"{speedup['transport']:>18.2f}x {speedup['end_to_end']:>20.2f}x"
            )

    print("Speedup = CPU time / Metal time; below 1 means Metal is slower.")
    print(f"Statistical checks: {'PASS' if report['passed'] else 'FAIL — do not accept speedups'}")
    print(destination / "report.json")

    if not report["passed"]:
        raise SystemExit(1)

    return destination


def main() -> None:
    """Parse an explicit fixed workload; require enough rounds for error estimates."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", choices=CASES)
    parser.add_argument("--threads", nargs="+", type=int, default=[1, 4])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--rounds", type=int, default=32)
    parser.add_argument("--photons-per-column", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--metal-dispatch-size", type=int, default=65536)
    parser.add_argument("--metal-tally-shards", type=int, default=32)
    parser.add_argument(
        "--profile-metal", action="store_true", help="Record completed GPU command time in receipts"
    )
    args = parser.parse_args()

    if args.repeats < 1 or args.rounds < 31 or any(n < 1 for n in args.threads):
        parser.error("repeats and threads must be positive; rounds must be at least 31")

    if args.photons_per_column is not None and args.photons_per_column < 1:
        parser.error("photons per column must be positive")

    if not 1 <= args.metal_dispatch_size <= 65536 or not 1 <= args.metal_tally_shards <= 64:
        parser.error("Metal dispatch size must be 1..65536 and tally shards 1..64")

    args.threads = list(dict.fromkeys(args.threads))
    run_benchmark(args)


if __name__ == "__main__":
    main()
