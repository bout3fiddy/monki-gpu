# Build from the repository root for Runpod's Linux amd64 GPU workers.
# Runpod supplies the host NVIDIA driver when the worker starts.
FROM nvidia/cuda:12.8.1-base-ubuntu24.04

# Pin uv so builds use a reviewed release. 0.13.4 removed the affected
# eMudhra/emSign root; 0.13.5 also updates managed Python's pip/ensurepip bundles.
COPY --from=ghcr.io/astral-sh/uv:0.13.5 /uv /uvx /usr/local/bin/

# Fortran and NetCDF support the retained scientific code and table preparation.
# g++ builds the host bridge; NVRTC compiles transport for the assigned GPU at
# runtime. Driver development stubs allow image compilation without a GPU.
# Use Ubuntu's main/universe components and bound retries on package downloads.
RUN sed -i 's/Components: main restricted universe multiverse/Components: main universe/' /etc/apt/sources.list.d/ubuntu.sources \
    && apt-get -o Acquire::Retries=2 -o Acquire::http::Timeout=30 update \
    && apt-get -o Acquire::Retries=2 -o Acquire::http::Timeout=30 install -y --no-install-recommends \
    gfortran g++ libnetcdff-dev libnetcdf-dev ca-certificates \
    cuda-nvrtc-dev-12-8 cuda-cudart-dev-12-8 cuda-driver-dev-12-8 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/monki

# Keep managed Python at a stable image path. Copy dependencies into the venv
# so installed files remain independent of uv's cache and filesystem boundaries.
ENV UV_PYTHON_INSTALL_DIR=/opt/python UV_LINK_MODE=copy

# Cache dependency installation separately from source changes. --locked fails
# if the manifest and lock disagree; the lock includes certifi 2026.10.19.
COPY pyproject.toml uv.lock ./

# Stay on Python 3.14; uv 0.13.5 selects its updated managed interpreter build.
RUN uv python install 3.14 && uv sync --locked --extra runpod-worker --no-install-project

# Include scientific sources, reference fixtures and checks for remote validation.
# .dockerignore excludes local environments, run outputs and credential files.
COPY python python
COPY scripts scripts
COPY kernels kernels
COPY fortran fortran
COPY validation validation
COPY tests tests
COPY README.md ARCHITECTURE.md STYLEGUIDE.md AGENTS.md .clang-format ./

# Install MONKI once its sources exist, then build both Fortran variants and the
# CUDA host bridge. Real-device transport compilation happens at session setup.
RUN uv sync --locked --extra runpod-worker && uv run python scripts/build_cuda.py

# Start the queue handler with the locked Runpod SDK dependencies available.
CMD ["uv", "run", "--locked", "--extra", "runpod-worker", "python", "-m", "monki.cloud.worker"]
