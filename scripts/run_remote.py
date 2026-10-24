"""Submit, retrieve, or cancel MONKI simulations through Runpod."""

import argparse
import json
from pathlib import Path

from monki.cloud import RunpodServerless
from monki.cloud.jobs import save_artifacts, simulation_request
from monki.input import ReferenceCase


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoint", required=True)
    parser.add_argument("--execution-timeout", type=int, default=600)
    parser.add_argument("--ttl", type=int, default=3600)
    subparsers = parser.add_subparsers(dest="action", required=True)
    submit = subparsers.add_parser("submit")
    submit.add_argument("--case", type=Path, required=True)
    submit.add_argument("--seed", type=int, default=23274)
    submit.add_argument("--job-file", type=Path, required=True)
    validate = subparsers.add_parser("validate")
    validate.add_argument("--job-file", type=Path, required=True)
    fetch = subparsers.add_parser("fetch")
    fetch.add_argument("--job", required=True)
    fetch.add_argument("--output", type=Path, required=True)
    fetch.add_argument("--wait", type=float, default=900)
    cancel = subparsers.add_parser("cancel")
    cancel.add_argument("--job", required=True)
    args = parser.parse_args()
    provider = RunpodServerless(
        args.endpoint, execution_timeout=args.execution_timeout, ttl=args.ttl
    )

    if args.action in ("submit", "validate"):
        if args.action == "validate":
            request = {"version": 1, "operation": "validate"}
        else:
            case = ReferenceCase(**json.loads(args.case.read_text()))
            request = simulation_request(
                case, save_path="2boxclouds" not in case.input_file, seed=args.seed
            )

        # Reserve the tracking file before the billable submission. Never retry an
        # uncertain submission automatically; the endpoint console can reconcile it.
        with args.job_file.open("x") as record:
            json.dump(
                {"provider": "runpod", "endpoint": args.endpoint, "state": "submitting"},
                record,
            )
            record.flush()
            job = provider.submit(request)
            record.seek(0)
            record.truncate()
            json.dump(
                {
                    "provider": "runpod",
                    "endpoint": args.endpoint,
                    "job_id": job,
                    "request": request,
                },
                record,
                indent=2,
            )
            record.write("\n")

        print(job)
    elif args.action == "fetch":
        if args.output.exists():
            raise FileExistsError(args.output)

        result = provider.wait(args.job, timeout=args.wait)

        if result.get("version") == 1 and "validation" in result:
            args.output.mkdir(parents=True, exist_ok=False)
            (args.output / "validation.json").write_text(json.dumps(result, indent=2) + "\n")

            if result["validation"].get("passed") is not True:
                raise RuntimeError(f"CUDA validation failed; see {args.output / 'validation.json'}")
        else:
            save_artifacts(result, args.output)

        print(args.output)
    else:
        provider.cancel(args.job)
        print(f"Cancellation requested for {args.job}")


if __name__ == "__main__":
    main()
