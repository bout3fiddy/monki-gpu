# Runpod Serverless worker

The image executes a complete simulation per queue job, including Fortran table
preparation, CUDA transport, convergence, statistics and output. Provider clients
live in `python/monki/cloud/`; GPU backends contain no provider settings.

Build from the repository root:

```sh
docker build --platform linux/amd64 -f containers/cuda/Dockerfile -t monki-cuda:dev .
```

Publish the image to your chosen registry and create a Runpod queue-based
Serverless endpoint using its immutable image digest. Start with zero minimum
workers and one maximum worker, one GPU per worker and a bounded execution timeout.
The image includes the fixtures, native binaries and optional Runpod worker SDK.
No provider API key is embedded in the build or sent in simulation inputs.

The native host bridge is built during image creation. NVRTC compiles the shared
transport for the assigned GPU at session preparation. CUDA compilation and real
GPU validation are distinct checks; a container build alone does not establish
scientific agreement or runtime performance.

Run the full NVIDIA gate as a fixed Serverless diagnostic job after deploying:

```sh
uv run python scripts/run_remote.py --endpoint ENDPOINT_ID --execution-timeout 1200 validate \
  --job-file cuda-validation-job.json
```

Fetch that job with the same `fetch` command as simulations. It saves
`validation.json` and returns a failure if the canonical CUDA gate failed.

## Client

Expose `RUNPOD_API_KEY` in the launching shell. The temporary `.zshrc` setup stores
plaintext and exports it to shell children; use restrictive permissions and keep
it out of backups or dotfile repositories you share. MONKI strips this variable
from native worker and compiler environments. A scoped endpoint key limits access.
The client reads the environment when making a request, allowing a later secret
manager wrapper without API changes. Codex subprocess shells may need an explicit
interactive zsh invocation to load `.zshrc`; the library never sources shell files.

```sh
uv run python scripts/run_remote.py --endpoint ENDPOINT_ID submit \
  --case validation/cases/clear.json --job-file clear-job.json
uv run python scripts/run_remote.py --endpoint ENDPOINT_ID fetch \
  --job JOB_ID --output results/clear-cuda
uv run python scripts/run_remote.py --endpoint ENDPOINT_ID cancel --job JOB_ID
```

Submit saves the job identifier for later retrieval. A local wait timeout leaves
the remote job accessible; cancel requests cancellation and status confirms it.
Submission failures with an uncertain outcome are never automatically retried.
Runpod's async result retention is 30 minutes after completion; retrieve promptly.

The version-1 protocol accepts the three bundled fixtures, explicit seed and
optional AMF convergence. It bounds rounds and inline artifacts, rejects arbitrary
paths/commands, and returns checksummed NetCDF and JSON receipts. Output directories
must be new. Persistent object storage and larger streaming jobs are future work.

The `Provider` abstract class defines `submit`, `status`, `cancel`, and shared
bounded polling. Another provider implements these operations using the same job
schema and result artifacts; scientific execution remains unchanged.

References: [queue requests](https://docs.runpod.io/serverless/endpoints/send-requests)
and [worker handlers](https://docs.runpod.io/serverless/workers/handler-functions).
