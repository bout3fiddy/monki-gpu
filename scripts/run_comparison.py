"""Run original MONKI and Metal, then plot TOA reflectance and layer AMFs."""

import argparse
import json
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from time import perf_counter

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from build_metal import build_metal
from build_reference import build
from monki import ReferenceCase, ReferenceRunner, simulate
from monki.execution import Metal
from monki.validation import compare

ROOT = Path(__file__).resolve().parents[1]


def main():
    """Build, run both solvers with identical conditions, and save comparisons."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("clear", "layered-cloud"), default="clear")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--photons-per-column", type=int, default=2_000)
    parser.add_argument("--rounds", type=int, default=32)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.rounds <= 30:
        parser.error("Use at least 31 rounds to obtain radiance standard errors")

    if args.threads < 1 or args.photons_per_column % args.threads:
        parser.error("Photon count must be divisible by a positive thread count")

    case_data = json.loads((ROOT / f"validation/cases/{args.case}.json").read_text())
    case_data.update(photons_per_column=args.photons_per_column, rounds=args.rounds)
    case = ReferenceCase(**case_data)
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    output = args.output or ROOT / "results" / f"{args.case}-{stamp}"
    output.mkdir(parents=True, exist_ok=False)

    # Compile before timing so the table measures simulation execution.
    reference_executable = build()
    build_metal()

    # Run Victor's original main program in an isolated directory.
    runner = ReferenceRunner(
        reference_executable,
        ROOT / "fortran/monki/upstream",
        threads=args.threads,
    )
    started = perf_counter()
    reference = runner.run(case, output / "reference")
    reference_seconds = perf_counter() - started

    # The same scene and photon budget go through the persistent Metal backend.
    backend = Metal()
    started = perf_counter()
    metal = simulate(case, backend=backend, output=output / "metal.nc")
    metal_seconds = perf_counter() - started

    comparison = compare(reference, metal)
    report = {
        "case": asdict(case),
        "threads": args.threads,
        "device": backend.provenance["device"],
        "timing_scope": "single run; setup, transport and output; compilation excluded",
        "reference_seconds": reference_seconds,
        "metal_seconds": metal_seconds,
        "radiance_comparison": comparison.to_dict(),
    }
    plot_comparison(reference, metal, report, output)
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Radiance comparison: {'PASS' if comparison.passed else 'FAIL'}")
    print(f"Reference: {reference_seconds:.3f} s; Metal: {metal_seconds:.3f} s")
    print(output.resolve() / "comparison.png")
    return 0 if comparison.passed else 1


def plot_comparison(reference, metal, report, output):
    """Place simulation conditions, six scientific panels and a numeric table together."""

    plt.rcParams.update({"font.family": "monospace", "font.size": 9})
    figure = plt.figure(figsize=(13, 11), layout="constrained")
    grid = figure.add_gridspec(3, 3, height_ratios=(1, 1.3, 0.8))
    a, b = reference.values, metal.values
    mu0 = np.cos(np.deg2rad(a["theta0"].item()))
    rows = []

    # MONKI radiance assumes E0=pi; dividing by mu0 gives Stokes reflectance.
    for component, label in enumerate(("I", "Q", "U")):
        axis = figure.add_subplot(grid[0, component])
        x = a["radiance"][:, 0, 0, component, 0] / mu0
        y = b["radiance"][:, 0, 0, component, 0] / mu0

        for values, name, marker in ((a, "Reference", "o"), (b, "Metal", "x")):
            axis.errorbar(
                values["theta"],
                values["radiance"][:, 0, 0, component, 0] / mu0,
                yerr=values["std_radiance"][:, 0, 0, component, 0] / mu0,
                fmt=marker,
                capsize=4,
                label=name,
            )

        se = np.hypot(a["std_radiance"], b["std_radiance"])[:, 0, 0, component, 0] / mu0
        z = np.divide(abs(y - x), se, out=np.zeros_like(se), where=se > 0)
        z[(se == 0) & (x != y)] = np.inf
        rows.append(metric_row(f"Reflectance {label}", x, y, f"{z.max():.3f}"))
        axis.set(title=f"TOA reflectance {label} (±1 SE)", xlabel="Viewing zenith (°)")
        axis.legend()
        axis.grid(alpha=0.2)

    altitude = (a["alt"][:-1] + a["alt"][1:]) / 2000

    for view, angle in enumerate(a["theta"]):
        axis = figure.add_subplot(grid[1, view])
        x = a["AMF_layer"][view, :, 0, 0, 0]
        y = b["AMF_layer"][view, :, 0, 0, 0]
        axis.plot(x, altitude, label="Reference")
        axis.plot(y, altitude, "--", label="Metal")
        axis.set(title=f"Layer AMF · view {angle:g}°", xlabel="AMF", ylabel="Altitude (km)")
        axis.legend()
        axis.grid(alpha=0.2)
        rows.append(metric_row(f"AMF {angle:g}°", x, y, "—"))

    # Original outputs contain radiance SE only; AMF differences remain descriptive.
    report["metrics"] = rows
    table_axis = figure.add_subplot(grid[2, :])
    table_axis.axis("off")
    table = table_axis.table(
        cellText=rows,
        colLabels=["Quantity", "RMSE", "Max |difference|", "Max |difference| / combined SE"],
        cellLoc="center",
        loc="center",
        colWidths=[0.18, 0.17, 0.22, 0.43],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)
    cpu, gpu = report["reference_seconds"], report["metal_seconds"]
    passed = report["radiance_comparison"]["passed"]
    table_axis.set_title(
        f"Reference ({report['threads']} threads): {cpu:.3f} s | Metal: {gpu:.3f} s | "
        f"Reference/Metal: {cpu / gpu:.2f}×\n"
        f"Single run, setup + transport + output | Radiance check: {'PASS' if passed else 'FAIL'}"
    )
    case = report["case"]
    figure.suptitle(
        f"{case['name']} · 340 nm · solar zenith {case['solar_zenith']:g}° · "
        f"albedo {case['surface_albedo']:g}\n"
        f"{case['rounds']:,} rounds × {case['photons_per_column']:,} photons · "
        f"{report['device']}"
    )

    for extension in ("png", "svg"):
        figure.savefig(output / f"comparison.{extension}", dpi=180)

    plt.close(figure)


def metric_row(name, reference, candidate, z):
    """Use absolute differences so zero-crossing Q/U values remain interpretable."""

    delta = candidate - reference
    return [name, f"{np.sqrt(np.mean(delta**2)):.3e}", f"{np.max(abs(delta)):.3e}", z]


if __name__ == "__main__":
    raise SystemExit(main())
