"""Link a persistent CPU worker against Victor's unchanged transport objects."""

import hashlib
import json
import shutil
import subprocess
from pathlib import Path

from build_reference import ROOT, SOURCES, build
from monki.environment import native_environment


def build_cpu(*, save_path: bool = True) -> Path:
    """Reuse the verified reference objects; replace only the program orchestration."""

    reference = build(save_path=save_path)
    reference_dir = reference.parent
    receipt = json.loads((reference_dir / "build.json").read_text())
    source = ROOT / "kernels/cpu/worker.f90"
    directory = ROOT / "build/cpu" / reference_dir.name
    executable = directory / "monki_cpu"
    identity = {
        "reference_build": receipt,
        "worker_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    }
    manifest = directory / "build.json"

    if executable.exists() and manifest.exists():
        previous = json.loads(manifest.read_text())

        if (
            previous.get("identity") == identity
            and previous.get("binary_sha256") == hashlib.sha256(executable.read_bytes()).hexdigest()
        ):
            return executable

    directory.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(reference_dir / "max_incl", directory / "max_incl")
    local_source = directory / "worker.f90"
    shutil.copyfile(source, local_source)
    compiler = receipt["commands"][0][0]
    settings = receipt["identity"]
    compile_command = [
        compiler,
        *settings["flags"],
        *settings["includes"],
        "-I",
        str(reference_dir),
        "-I",
        str(directory),
        "-J",
        str(directory),
        "-c",
        str(local_source),
    ]
    subprocess.run(compile_command, cwd=directory, check=True, env=native_environment())

    # Every physics/helper object comes from the reference build unchanged.
    # Only main.o is replaced with the persistent worker's orchestration.
    objects = [
        str(reference_dir / Path(name).with_suffix(".o").name)
        for name in SOURCES
        if name != "main.f90"
    ]
    link_command = [
        compiler,
        "-fopenmp",
        *objects,
        str(directory / "worker.o"),
        *settings["libraries"],
        "-o",
        str(executable),
    ]
    subprocess.run(link_command, cwd=directory, check=True, env=native_environment())
    manifest.write_text(
        json.dumps(
            {
                "identity": identity,
                "commands": [compile_command, link_command],
                "binary_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
            },
            indent=2,
        )
        + "\n"
    )
    return executable


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-path", action="store_true")
    arguments = parser.parse_args()
    print(build_cpu(save_path=not arguments.no_path))
