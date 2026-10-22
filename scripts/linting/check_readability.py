"""Enforce zdisamar-inspired paragraph boundaries using Python's syntax tree."""

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCOPES = ("python", "scripts", "tests", "validation")
COMPOUND = (ast.If, ast.For, ast.AsyncFor, ast.While, ast.With, ast.AsyncWith, ast.Try)


def findings(source: str) -> list[tuple[int, str]]:
    """Require space between logical stages without inspecting vendor source."""

    tree = ast.parse(source)
    lines = source.splitlines()
    problems = []

    for node in ast.walk(tree):
        if (
            isinstance(node, ast.ImportFrom)
            and node.module == "__future__"
            and any(alias.name == "annotations" for alias in node.names)
        ):
            problems.append((node.lineno, "do not import future annotations"))

        for field in ("body", "orelse", "finalbody"):
            statements = getattr(node, field, None)

            if not isinstance(statements, list):
                continue

            for before, after in zip(statements, statements[1:], strict=False):
                if not isinstance(before, ast.stmt) or not isinstance(after, ast.stmt):
                    continue

                docstring = (
                    isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                    and before is statements[0]
                    and isinstance(before, ast.Expr)
                    and isinstance(before.value, ast.Constant)
                    and isinstance(before.value.value, str)
                )
                separate = docstring or isinstance(before, COMPOUND) or isinstance(after, COMPOUND)

                if separate and not any(
                    not line.strip() for line in lines[before.end_lineno : after.lineno - 1]
                ):
                    problems.append((after.lineno, "add a blank line between logical stages"))

    return sorted(set(problems))


def main() -> int:
    """Check authored Python files, excluding generated validation runs."""

    failures = 0

    for scope in SCOPES:
        for path in sorted((ROOT / scope).rglob("*.py")):
            if "runs" in path.relative_to(ROOT).parts:
                continue

            for line, message in findings(path.read_text()):
                print(f"{path.relative_to(ROOT)}:{line}: {message}")
                failures += 1

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
