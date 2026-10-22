"""Plot saved CPU–Metal validation outputs without rerunning transport."""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import netCDF4
import numpy as np


def main():
    """Read a validation run and write reproducible PNG/SVG comparison figures."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    output = args.run / "plots"
    output.mkdir(exist_ok=True)
    cases = ("clear", "layered-cloud", "box-clouds")
    data = {}
    provenance = {}

    for case in cases:
        data[case] = {}

        for backend in ("cpu", "metal"):
            path = args.run / f"{case}-{backend}.nc"
            provenance[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()

            with netCDF4.Dataset(path) as dataset:
                data[case][backend] = {
                    name: np.asarray(variable[:]) for name, variable in dataset.variables.items()
                }

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    plot_angles(data, output)
    plot_maps(data["box-clouds"], output)
    plot_amf(data, output)
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(output.resolve())


def save(figure, output, name):
    """Keep a preview and a vector version suitable for a research report."""

    for extension in ("png", "svg"):
        figure.savefig(output / f"{name}.{extension}", dpi=180, bbox_inches="tight")

    plt.close(figure)


def plot_angles(data, output):
    """Follow Victor's radiance/cos(theta0) conversion for Stokes reflectance."""

    figure, axes = plt.subplots(2, 3, figsize=(12, 7), layout="constrained")

    for row, case in enumerate(("clear", "layered-cloud")):
        for backend, marker in (("cpu", "o"), ("metal", "x")):
            values = data[case][backend]
            mu0 = np.cos(np.deg2rad(values["theta0"].item()))
            reflectance = values["radiance"][:, 0, 0, :, 0] / mu0
            error = values["std_radiance"][:, 0, 0, :, 0] / mu0

            for component, label in enumerate(("I", "Q", "U")):
                axis = axes[row, component]
                axis.errorbar(
                    values["theta"],
                    reflectance[:, component],
                    yerr=error[:, component],
                    fmt=marker,
                    capsize=5,
                    label=backend.upper(),
                    alpha=0.85,
                )
                axis.set(
                    title=f"{case} · {label}",
                    xlabel="Viewing zenith (°)",
                    ylabel="TOA Stokes reflectance (dimensionless)",
                    xticks=[0, 30, 60],
                )
                axis.grid(alpha=0.2)

    axes[0, 0].legend()
    figure.suptitle(
        "CPU vs Metal · 340 nm · solar zenith 45°\n"
        "Markers are the three simulated views; error bars are ±1 standard error"
    )
    save(figure, output, "stokes-reflectance")


def plot_maps(data, output):
    """Show cloud structure and residuals without hiding spatial discrepancies."""

    cpu, metal = data["cpu"], data["metal"]
    mu0 = np.cos(np.deg2rad(cpu["theta0"].item()))
    reference = cpu["radiance"][:, :, :, 0, 0] / mu0
    candidate = metal["radiance"][:, :, :, 0, 0] / mu0
    uncertainty = np.hypot(cpu["std_radiance"], metal["std_radiance"])[:, :, :, 0, 0] / mu0
    residual = np.divide(
        candidate - reference,
        uncertainty,
        out=np.full_like(reference, np.nan),
        where=uncertainty > 0,
    )
    height, width = reference.shape[1:]
    extent = (0, width * cpu["x_unit"].item(), 0, height * cpu["y_unit"].item())
    figure, axes = plt.subplots(3, 3, figsize=(13, 9), layout="constrained")

    # Share physical colour limits across both backends and all viewing angles.
    lower, upper = min(reference.min(), candidate.min()), max(reference.max(), candidate.max())
    limit = max(1.0, np.nanmax(np.abs(residual)))

    for view, angle in enumerate(cpu["theta"]):
        for column, (field, label) in enumerate(
            ((reference, "CPU"), (candidate, "Metal"), (residual, "(Metal − CPU) / combined SE"))
        ):
            axis = axes[view, column]
            difference = column == 2
            artist = axis.imshow(
                field[view],
                origin="lower",
                extent=extent,
                cmap="RdBu_r" if difference else "viridis",
                vmin=-limit if difference else lower,
                vmax=limit if difference else upper,
            )
            axis.set(title=f"{label} · view {angle:g}°", xlabel="x (km)", ylabel="y (km)")
            figure.colorbar(
                artist,
                ax=axis,
                shrink=0.8,
                label="Combined SE units" if difference else "TOA reflectance I",
            )

    figure.suptitle(
        "3D box clouds · 340 nm · solar zenith 45°\n"
        "Independent CPU/Metal Monte Carlo estimates; 32 rounds per backend"
    )
    save(figure, output, "cloud-reflectance-maps")


def plot_amf(data, output):
    """Compare layer sensitivity profiles; saved files do not contain AMF errors."""

    figure, axes = plt.subplots(1, 2, figsize=(10, 6), layout="constrained")
    colours = ("#0072B2", "#D55E00", "#009E73")

    for axis, case in zip(axes, ("clear", "layered-cloud"), strict=True):
        for backend, style in (("cpu", "-"), ("metal", "--")):
            values = data[case][backend]
            altitude = (values["alt"][:-1] + values["alt"][1:]) / 2000

            for view, angle in enumerate(values["theta"]):
                axis.plot(
                    values["AMF_layer"][view, :, 0, 0, 0],
                    altitude,
                    linestyle=style,
                    color=colours[view],
                    label=f"{backend.upper()} {angle:g}°",
                )

        axis.set(
            title=case,
            xlabel="Layer air mass factor (dimensionless)",
            ylabel="Layer midpoint altitude (km)",
        )
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8)

    figure.suptitle("Layer air mass factors · 340 nm")
    save(figure, output, "layer-amf")


if __name__ == "__main__":
    main()
