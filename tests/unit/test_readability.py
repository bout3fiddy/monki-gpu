"""The style gate catches dense stages and leaves readable source alone."""

import runpy
from pathlib import Path

CHECKER = Path(__file__).resolve().parents[2] / "scripts/linting/check_readability.py"
findings = runpy.run_path(str(CHECKER))["findings"]


def test_dense_guard_and_return_are_rejected():
    source = "def run(x):\n    if x < 0:\n        raise ValueError()\n    return x\n"
    assert findings(source) == [(4, "add a blank line between logical stages")]


def test_separated_guard_and_return_are_accepted():
    source = "def run(x):\n    if x < 0:\n        raise ValueError()\n\n    return x\n"
    assert findings(source) == []


def test_future_annotations_rule_is_enforced():
    assert findings("from __future__ import annotations\n") == [
        (1, "do not import future annotations")
    ]
