# Code style

Adapted from [zdisamar](docs/style-reference/zdisamar-STYLEGUIDE.md).

- Separate logical stages with blank lines, including after function docstrings.
- Use domain names, spaced arithmetic and one argument per line in wrapped calls.
- Keep public entry points before helpers; describe their inputs and outputs.
- Comment what meaningful blocks do and why. Explain units, signs, guards,
  tolerances, source equations and ownership near the code; avoid narrating assignments.
- Keep parsing, file I/O, credentials and scheduling outside scientific routines.
- Pass explicit tables and named buffers; avoid catch-all mutable state.
- Preserve vendored files byte-for-byte. Keep adapted routines traceable to originals.

Python requires 3.14+, uv and a 100-column line limit. Ruff handles formatting and
E/F/I/UP/B/SIM rules; the AST checker enforces paragraph spacing and rejects future
annotations. Basedpyright checks types.

New Fortran uses explicit kinds, `implicit none`, declared `intent` and named modules.
GPU changes must explain precision, RNG, reductions and geometry differences.
Shared transport changes require validation on both Metal and CUDA.

```sh
uv run ruff format python scripts tests validation
uv run ruff check python scripts tests validation
uv run python scripts/linting/check_readability.py
uv run basedpyright
uv run clang-format -i kernels/metal/bridge.m kernels/metal/transport.metal
uv run python scripts/verify.py
```

The final command is the canonical gate: hashes, formatting, lint, types, builds,
tests and scientific comparisons. Pre-commit runs static checks only.
