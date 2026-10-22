"""Independent-seed AMF convergence study with retained uncertainty and round tallies."""

import argparse
import json
from dataclasses import replace
from datetime import UTC, datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from monki import ReferenceCase, simulate
from monki.execution import CPU, Metal
from path_statistics import PathRecorder

ROOT = Path(__file__).resolve().parents[1]


def main():
    """Compare three independent streams at 64k and 1.024M photons per backend."""

    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    destination = ROOT / "validation/runs" / f"amf-study-{stamp}"
    destination.mkdir()
    entries = []

    for name in ("clear", "layered-cloud"):
        original = ReferenceCase(**json.loads((ROOT / f"validation/cases/{name}.json").read_text()))

        for photons, rounds in ((2000, 32), (8000, 128)):
            case = replace(original, photons_per_column=photons, rounds=rounds)

            for seed in (23274, 714103, 19088743):
                pair = {}

                for label, backend_type in (("cpu", CPU), ("metal", Metal)):
                    recorder = PathRecorder(backend_type(seed=seed))
                    prefix = f"{name}-{photons * rounds}-{seed}-{label}"
                    result = simulate(case, recorder, destination / f"{prefix}.nc")
                    error = recorder.standard_error(result)
                    amf = result.values["AMF_layer"][:, :, 0, 0, 0]
                    error = error[:, :, 0, 0, 0]
                    # Preserve raw evidence: rare high-weight paths can make a
                    # finite-round normal approximation unreliable.
                    np.savez_compressed(
                        destination / f"{prefix}-rounds.npz",
                        paths=np.stack([batch.weighted_paths for batch in recorder.rounds]),
                        weights=np.stack([batch.weights for batch in recorder.rounds]),
                        amf=amf,
                        standard_error=error,
                    )
                    pair[label] = {"amf": amf.tolist(), "se": error.tolist()}
                    altitude = result.values["alt"]

                delta = np.asarray(pair["metal"]["amf"]) - np.asarray(pair["cpu"]["amf"])
                combined = np.hypot(pair["metal"]["se"], pair["cpu"]["se"])
                entry = {
                    "case": name,
                    "photons": photons * rounds,
                    "seed": seed,
                    "pair": pair,
                    "max_absolute_difference": float(np.max(np.abs(delta))),
                    "max_combined_se": float(np.max(np.abs(delta) / combined)),
                    "altitude_km": ((altitude[:-1] + altitude[1:]) / 2000).tolist(),
                }
                entries.append(entry)
                (destination / "report.json").write_text(json.dumps(entries, indent=2) + "\n")
                print(
                    name,
                    photons * rounds,
                    seed,
                    entry["max_absolute_difference"],
                    entry["max_combined_se"],
                    flush=True,
                )

    plot(entries, destination)
    print(destination, flush=True)


def plot(entries, destination, benchmark=None):
    """Pool equal-work independent replicates and show two-standard-error bands."""

    plt.rcParams.update({"font.family": "DejaVu Sans Mono"})
    figure = plt.figure(figsize=(12, 12 if benchmark else 10), layout="constrained")
    grid = figure.add_gridspec(
        3 if benchmark else 2, 2, height_ratios=[1, 1, 0.38] if benchmark else [1, 1]
    )
    axes = np.asarray(
        [[figure.add_subplot(grid[row, column]) for column in range(2)] for row in range(2)]
    )
    summary = []

    for row, name in enumerate(("clear", "layered-cloud")):
        for column, count in enumerate((64000, 1024000)):
            selected = [
                item for item in entries if item["case"] == name and item["photons"] == count
            ]
            axis = axes[row, column]
            pooled = {}

            for label, style in (("cpu", "-"), ("metal", "--")):
                values = np.asarray([item["pair"][label]["amf"] for item in selected])
                errors = np.asarray([item["pair"][label]["se"] for item in selected])
                mean = values.mean(axis=0)
                error = np.sqrt((errors**2).sum(axis=0)) / len(selected)
                pooled[label] = (mean, error)
                altitude = selected[0]["altitude_km"]
                (line,) = axis.plot(mean[2], altitude, style, label=label.upper())
                axis.fill_betweenx(
                    altitude,
                    mean[2] - 2 * error[2],
                    mean[2] + 2 * error[2],
                    color=line.get_color(),
                    alpha=0.18,
                )

            delta = pooled["metal"][0] - pooled["cpu"][0]
            combined = np.hypot(pooled["metal"][1], pooled["cpu"][1])
            # Metrics use the 33 layers at the plotted 60-degree view.
            view_delta = delta[2]
            rmse = np.sqrt(np.mean(view_delta**2))
            maximum = np.max(np.abs(view_delta))
            relative = 100 * np.max(np.abs(view_delta) / np.abs(pooled["cpu"][0][2]))
            axis.text(
                0.04,
                0.76,
                f"RMSE  {rmse:.4f}\nMax |Δ|  {maximum:.4f}\nMax relative Δ  {relative:.2f}%",
                transform=axis.transAxes,
                va="top",
                fontsize=10,
                linespacing=1.6,
            )
            summary.append(
                {
                    "case": name,
                    "photons_per_seed": count,
                    "max_abs_difference": float(np.max(np.abs(delta))),
                    "max_combined_se": float(np.max(np.abs(delta) / combined)),
                    "target_32_5_km_cpu": float(pooled["cpu"][0][2, 27]),
                    "target_32_5_km_metal": float(pooled["metal"][0][2, 27]),
                    "target_combined_se": float(combined[2, 27]),
                }
            )
            scene_label = "Clear sky" if name == "clear" else "Layered cloud"
            axis.set(
                title=f"{scene_label} — {count:,} photons/seed",
                xlabel="Layer AMF",
                ylabel="Altitude (km)",
            )
            axis.grid(alpha=0.2)
            axis.legend()

    figure.suptitle(
        "Layer air mass factors\n"
        "340 nm · solar zenith 45° · viewing zenith 60°\n"
        "Surface albedo 0.05 · mean of 3 seeds",
        fontsize=14,
    )

    if benchmark is not None:
        timing_axis = figure.add_subplot(grid[2, :])
        timing_axis.axis("off")
        timing_axis.set_title("Execution time on Apple M4", loc="left", fontsize=13)
        rows = []

        for case in benchmark["cases"]:
            name = case["case"]["name"]
            label = {
                "clear": "Clear sky",
                "layered-cloud": "Layered cloud",
                "box-clouds": "3D box clouds",
            }[name]
            medians = case["medians"]
            count = case["samples"]["metal"][0]["total_photons"]
            rows.append(
                [
                    label,
                    f"{count:,}",
                    "Yes" if case["save_path"] else "No",
                    *[
                        f"{medians[key]['end_to_end_seconds']:.3f} s"
                        for key in ("cpu-1", "cpu-4", "metal")
                    ],
                    f"{case['speedups']['cpu-4']['end_to_end']:.2f}×",
                ]
            )

        table = timing_axis.table(
            cellText=rows,
            colLabels=[
                "Scene",
                "Photons",
                "AMFs",
                "CPU, 1 thread",
                "CPU, 4 threads",
                "Metal",
                "Speedup vs CPU-4",
            ],
            bbox=[0, 0.24, 1, 0.73],
            cellLoc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(9)

        for (row, _column), cell in table.get_celld().items():
            cell.set_edgecolor("#d7dde3")

            if row == 0:
                cell.set_facecolor("#eaf0f5")
                cell.set_text_props(weight="bold")

    for extension in ("png", "svg"):
        figure.savefig(destination / f"amf-convergence.{extension}", dpi=180)

    (destination / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    plt.close(figure)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot-run", type=Path, help="Replot an existing study without transport")
    parser.add_argument("--benchmark", type=Path, help="Saved benchmark report to include")
    args = parser.parse_args()

    if args.plot_run:
        entries = json.loads((args.plot_run / "report.json").read_text())
        benchmark = json.loads(args.benchmark.read_text()) if args.benchmark else None
        plot(entries, args.plot_run, benchmark)
    else:
        main()
