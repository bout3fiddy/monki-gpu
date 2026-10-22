"""Compare a future solver's NetCDF output with an aligned reference result."""

import argparse
import json
from pathlib import Path

from monki import SimulationResult
from monki.validation import compare

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--mode", choices=("statistical", "deterministic"), default="statistical")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    result = compare(
        SimulationResult.load(args.reference), SimulationResult.load(args.candidate), mode=args.mode
    )
    text = json.dumps(result.to_dict(), indent=2) + "\n"

    if args.report:
        args.report.write_text(text)

    print(text, end="")
    raise SystemExit(0 if result.passed else 1)
