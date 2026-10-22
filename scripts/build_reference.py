"""Compile an isolated CPU reference without modifying vendored sources."""

import hashlib
import json
import os
import shlex
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VENDOR = ROOT / "fortran/monki"
BUILD = ROOT / "build/reference"
SOURCES = [
    "modules/write_messages_mod.f90",
    "modules/physical_parameters_mod.f90",
    "modules/new_data_types_mod.f90",
    "modules/main_routines_general_mod.f90",
    "modules/main_routines_cartesian_mod.f90",
    "modules/ncdf_utilities_mod.f90",
    "modules/interpol_sp.f",
    "modules/get_O3_crabs.f90",
    "modules/mc_routines_general_mod.f90",
    "modules/mc_routines_cartesian_mod.f90",
    "mc.f90",
    "main.f90",
]


def verify_vendor():
    """Check both file contents and inventory against the supplied-source manifest."""

    manifest = json.loads((VENDOR / "manifest.json").read_text())
    expected = manifest["files"]
    actual = {
        str(p.relative_to(VENDOR / "upstream"))
        for p in (VENDOR / "upstream").rglob("*")
        if p.is_file()
    }

    if actual != set(expected):
        raise RuntimeError("Vendored file inventory differs from manifest")

    for relative, digest in expected.items():
        if hashlib.sha256((VENDOR / "upstream" / relative).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"Vendored file changed: {relative}")

    return manifest


def build(*, save_path=True, save_scats=False):
    """Compile a diagnostic-specific reference in build/, keeping upstream untouched."""

    build_dir = BUILD / ("paths" if save_path else "radiance")

    if save_scats:
        build_dir = build_dir.with_name(build_dir.name + "-scats")

    verify_vendor()
    compiler = shutil.which(os.environ.get("FC", "gfortran"))
    config = shutil.which("nf-config")

    if not compiler or not config:
        raise RuntimeError("Install gfortran and NetCDF Fortran (nf-config) before building")

    compiler_version = subprocess.check_output([compiler, "--version"], text=True).splitlines()[0]
    flags = [
        "-O2",
        "-g",
        "-fopenmp",
        "-finit-local-zero",
        "-ffree-line-length-none",
        "-fcheck=bounds",
        "-ffpe-trap=zero,overflow,invalid",
    ]
    includes = shlex.split(subprocess.check_output([config, "--fflags"], text=True))
    libraries = shlex.split(subprocess.check_output([config, "--flibs"], text=True))
    # Some nf-config installations omit the separate NetCDF C library directory.
    nc_config = shutil.which("nc-config")

    if nc_config:
        libraries += shlex.split(subprocess.check_output([nc_config, "--libs"], text=True))

    # Covers the retained 30x20 two-box scene, 35 layers and 3 detector angles.
    limits = {"nx": 30, "ny": 20, "nz": 35, "nmat": 4, "na": 3}
    max_include = (
        "\n".join(
            f"INTEGER, PARAMETER :: max_{name}{suffix} = {value}"
            for suffix, enabled in (("", True), ("_1", save_scats), ("_2", save_path))
            for name, maximum in limits.items()
            for value in (maximum if enabled else 1,)
        )
        + "\n"
    )

    # Reuse a binary only when its toolchain/settings, vendor manifest, and own hash
    # still match. This avoids stale validation without rebuilding on every test run.
    identity = {
        "compiler": compiler_version,
        "flags": flags,
        "includes": includes,
        "libraries": libraries,
        "max_incl": max_include,
        "manifest_sha256": hashlib.sha256((VENDOR / "manifest.json").read_bytes()).hexdigest(),
    }
    binary = build_dir / "monki"
    receipt = build_dir / "build.json"

    if binary.exists() and receipt.exists():
        previous = json.loads(receipt.read_text())

        if (
            previous.get("identity") == identity
            and previous.get("binary_sha256") == hashlib.sha256(binary.read_bytes()).hexdigest()
        ):
            return binary

    build_dir.mkdir(parents=True, exist_ok=True)
    (build_dir / "max_incl").write_text(max_include)
    commands = []
    objects = []

    # Copy source next to the generated max_incl. Fortran INCLUDE otherwise finds
    # the upstream dimensions first; numerical source files remain byte-identical.
    for source in SOURCES:
        original = VENDOR / "upstream/PROGRAM" / source
        local = build_dir / Path(source).name
        shutil.copyfile(original, local)
        command = [
            compiler,
            *flags,
            *includes,
            "-I",
            str(build_dir),
            "-J",
            str(build_dir),
            "-c",
            str(local),
        ]
        commands.append(command)
        subprocess.run(command, cwd=build_dir, check=True)
        objects.append(str(local.with_suffix(".o")))

    command = [compiler, "-fopenmp", *objects, *libraries, "-o", str(binary)]
    commands.append(command)
    subprocess.run(command, cwd=build_dir, check=True)
    receipt.write_text(
        json.dumps(
            {
                "identity": identity,
                "commands": commands,
                "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
            },
            indent=2,
        )
        + "\n"
    )
    return binary


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--no-path", action="store_true")
    parser.add_argument("--scats", action="store_true")
    args = parser.parse_args()
    print(build(save_path=not args.no_path, save_scats=args.scats))
