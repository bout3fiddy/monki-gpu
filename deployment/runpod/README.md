# MONKI Runpod deployment

The Linux amd64 worker image has been uploaded to private AWS ECR in `eu-north-1`:

```
824648489183.dkr.ecr.eu-north-1.amazonaws.com/monki-cuda:cuda-c7f6fe9b6885
```

The repository uses immutable tags. The tag identifies local Docker image
`sha256:c7f6fe9b68851a917c4fdc27ac2c3c80847ee163e949ff297689c5423b3b9920`.
ECR confirms manifest digest
`sha256:c7f6fe9b68851a917c4fdc27ac2c3c80847ee163e949ff297689c5423b3b9920`.

`ecr-pull-policy.json` was applied on 2026-10-22 after user approval and verified
with `aws ecr get-repository-policy`. It grants image pull operations on this
repository to Runpod's two published deployment roles. It grants no push or
delete operations.

Registry registration is unresolved: both the Serverless wizard and Credentials
page returned to the AWS ECR form after Save, with no visible error or saved
credential confirmation. The prepared credential name is `monki-ecr`. Confirm
registration succeeds before deploying; the AWS policy alone does not establish
that Runpod can pull the image.

Target configuration: queue endpoint `monki-cuda`, zero active workers, one maximum
worker, one GPU per worker, bounded execution timeout, no network volume. No
endpoint has been created yet. The Runpod balance was $0 at setup time; no paid
jobs or credit purchases have been made.

See [worker build and client commands](../../containers/cuda/README.md).
Runpod's ECR guidance is at
https://docs.runpod.io/tutorials/pods/use-private-ecr-images ; its Serverless
console support must be verified during setup.
