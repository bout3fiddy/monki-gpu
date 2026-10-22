# MONKI development

## Primary design goal

Reuse the majority of Victor J.H. Trees's MONKI code so that, when the project is
handed back to him, it remains easy and familiar for him to understand and modify.

- Prefer incremental adaptation of the existing Fortran over rewriting the physics
  in another language or introducing a substantially different solver structure.
- Preserve recognizable scientific routines, terminology, equations, and comments
  wherever practical. Keep necessary changes traceable to their original routines.
- Put Python interfaces, cloud execution, and backend orchestration around the
  scientific code; avoid forcing those concerns into Victor's physics routines.
- Choose GPU tooling and abstractions with code reuse and researcher familiarity
  as primary criteria. Justify necessary departures with correctness or measured
  performance evidence.
- Apply adaptations outside the immutable vendored reference snapshot.

## Writing rules

- State facts directly in documentation, comments, docstrings, output text and replies.
- Avoid contrastive framing such as “A, not B”, “not X but Y”, and rhetorical
  “rather than” comparisons. Describe the behavior, scope or limitation directly.

## Development rules

- Read ARCHITECTURE.md for the target design; README.md describes what exists today.
- Read STYLEGUIDE.md before changing code. Use blank lines between logical stages.
- Include a healthy amount of what-and-why comments, especially around scientific assumptions.
- Use uv for all Python commands and dependency changes.
- Require Python 3.14 or newer; keep package and tooling version targets aligned.
- Preserve fortran/monki/upstream byte-for-byte. Its manifest records provenance and hashes.
- Never overwrite supplied reference outputs or silently regenerate golden results.
- Keep file parsing, cloud credentials, and job scheduling outside scientific routines.
- Pass explicit physics tables and named write buffers through native routines.
- A requested unsupported backend or physics option must fail explicitly.
- Keep legacy input support scoped to validated fixtures until a full typed Scene exists.
- Tests must assert scientific or interface invariants; statistical tolerances must be explicit.
- Run `uv run python scripts/verify.py` for the canonical local gate. It requires gfortran and nf-config.
- No GPU speedup or physical equivalence claims without executed evidence.
