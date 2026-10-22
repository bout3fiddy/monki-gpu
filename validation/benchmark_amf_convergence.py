"""Measure end-to-end time to the same estimated AMF precision on CPU and Metal."""

import json
from dataclasses import asdict, replace
from datetime import UTC, datetime
from pathlib import Path
from time import perf_counter

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from monki import AMFConvergence, ReferenceCase, SimulationResult, simulate
from monki.execution import CPU, Metal
from monki.validation import check_result, compare, select_views

ROOT = Path(__file__).resolve().parents[1]


def main():
    """Use sequential, warmed runs and retain both converged and capped outcomes."""

    destination = ROOT / "validation/runs" / datetime.now(UTC).strftime("amf-timing-%Y%m%dT%H%M%SZ")
    destination.mkdir()
    policy = AMFConvergence()
    report = {
        "policy": asdict(policy),
        "cpu_threads": 4,
        "photons_per_round": 8000,
        "max_rounds": 8192,
        "timing": "end-to-end, compilation excluded; one warmup per case/backend",
        "samples": [],
        "comparisons": [],
    }

    for name in ("clear", "layered-cloud"):
        case = ReferenceCase(**json.loads((ROOT / f"validation/cases/{name}.json").read_text()))
        case = replace(case, photons_per_column=8000, rounds=8192)
        fixture = "cloudfree" if name == "clear" else "1layercloud"
        golden = select_views(
            SimulationResult.load(
                ROOT
                / "fortran/monki/upstream/OUTPUT/testrun"
                / f"MONKI_MLS_{fixture}_340_340_test_paper_as_005_th0_045_ph0_000_ph_000.nc"
            ),
            tuple(case.view_zeniths),
        )

        for label, backend in (("cpu", CPU(threads=4)), ("metal", Metal())):
            simulate(replace(case, rounds=2), backend, destination / f"{name}-{label}-warmup.nc")

        for seed in (23274, 714103, 19088743):
            results = {}
            # Alternate order to reduce consistent ordering/thermal bias.
            labels = ("metal", "cpu") if seed == 714103 else ("cpu", "metal")

            for label in labels:
                backend = CPU(threads=4, seed=seed) if label == "cpu" else Metal(seed=seed)
                path = destination / f"{name}-{seed}-{label}.nc"
                print(f"START {name} {seed} {label}", flush=True)
                start = perf_counter()
                result = simulate(case, backend, path, convergence=policy)
                elapsed = perf_counter() - start
                receipt = json.loads(path.with_suffix(".run.json").read_text())
                # The archived product predates the new AMF error variable.
                # Remove that extra field only for the legacy inventory check.
                legacy = replace(
                    result,
                    values={k: v for k, v in result.values.items() if k != "std_AMF_layer"},
                    dimensions={k: v for k, v in result.dimensions.items() if k != "std_AMF_layer"},
                    attributes={k: v for k, v in result.attributes.items() if k != "std_AMF_layer"},
                )
                archived_check = compare(golden, legacy)

                entry = {
                    "case": name,
                    "seed": seed,
                    "backend": label,
                    "seconds": elapsed,
                    "status": receipt["status"],
                    "rounds": receipt["completed_rounds"],
                    "photons": receipt["total_photons"],
                    "max_relative_se": receipt["amf_convergence"]["checks"][-1]["max_relative_se"],
                    "invariant_failures": check_result(result),
                    "archived_radiance_check": archived_check.to_dict(),
                    "checks": receipt["amf_convergence"]["checks"],
                }
                report["samples"].append(entry)
                results[label] = result
                print({key: value for key, value in entry.items() if key != "checks"}, flush=True)
                (destination / "report.json").write_text(json.dumps(report, indent=2) + "\n")

            cpu, metal = results["cpu"], results["metal"]
            delta = metal.values["AMF_layer"] - cpu.values["AMF_layer"]
            combined = np.hypot(cpu.values["std_AMF_layer"], metal.values["std_AMF_layer"])
            report["comparisons"].append(
                {
                    "case": name,
                    "seed": seed,
                    "max_absolute_difference": float(np.max(abs(delta))),
                    "max_combined_se": float(np.max(abs(delta) / combined)),
                }
            )
            (destination / "report.json").write_text(json.dumps(report, indent=2) + "\n")

    plot(report, destination)
    print(destination, flush=True)


def plot(report, destination):
    """Show precision histories and measured time, retaining seed-to-seed spread."""

    plt.rcParams.update({"font.family": "DejaVu Sans Mono", "font.size": 10})
    figure = plt.figure(figsize=(12, 8), layout="constrained")
    grid = figure.add_gridspec(2, 2, height_ratios=[1, 0.35])
    rows = []

    for column, name in enumerate(("clear", "layered-cloud")):
        axis = figure.add_subplot(grid[0, column])

        for label, colour in (("cpu", "#0072B2"), ("metal", "#D55E00")):
            selected = [s for s in report["samples"] if s["case"] == name and s["backend"] == label]

            for index, sample in enumerate(selected):
                checks = sample["checks"]
                axis.plot(
                    [c["rounds"] * report["photons_per_round"] for c in checks],
                    [100 * c["max_relative_se"] for c in checks],
                    color=colour,
                    alpha=0.65,
                    label=label.upper() if index == 0 else None,
                )

            seconds = [s["seconds"] for s in selected]
            photons = [s["photons"] for s in selected]
            rows.append(
                [
                    name,
                    "CPU (4 threads)" if label == "cpu" else "Metal",
                    f"{np.median(seconds):.3f}",
                    f"{min(seconds):.3f}–{max(seconds):.3f}",
                    f"{np.median(photons):,.0f}",
                    f"{sum(s['status'] == 'converged' for s in selected)}/3",
                ]
            )

        axis.axhline(1, color="black", linestyle=":", label="1% SE")
        axis.set(
            xscale="log",
            yscale="log",
            xlabel="Photons",
            ylabel="Maximum relative AMF SE (%)",
            title="Clear sky" if name == "clear" else "Layered cloud",
        )
        axis.grid(alpha=0.2)
        axis.legend()

    table_axis = figure.add_subplot(grid[1, :])
    table_axis.axis("off")
    table_axis.set_title("Apple M4 · end-to-end time · three seeds", loc="left")
    table = table_axis.table(
        cellText=rows,
        colLabels=["Scene", "Backend", "Median (s)", "Range (s)", "Median photons", "Converged"],
        cellLoc="center",
        bbox=[0, 0, 1, 0.95],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    figure.suptitle(
        "AMF precision stopping\n340 nm · solar zenith 45° · views 0°, 30°, 60°\n"
        "Albedo 0.05 · 33 layers · 8,000 photons/round · target SE < 1%",
        fontsize=13,
    )

    for extension in ("png", "svg"):
        figure.savefig(destination / f"amf-timing.{extension}", dpi=180)

    plt.close(figure)


if __name__ == "__main__":
    main()
