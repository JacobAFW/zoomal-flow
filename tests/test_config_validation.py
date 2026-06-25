#!/usr/bin/env python3
"""
Negative-validation tests for the agnostic config schema.

Runs both included bad configs through validate_config.load_and_validate
and asserts each one fails with a useful message. Also confirms the
template config.yaml and cohort.example.yaml validate cleanly.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
AGNOSTIC = HERE.parent
sys.path.insert(0, str(AGNOSTIC / "scripts" / "py"))

import jsonschema  # noqa: E402
from validate_config import load_and_validate  # noqa: E402

SCHEMA = AGNOSTIC / "config" / "schema" / "config.schema.yaml"


def _expect_invalid(path, hint: str):
    try:
        load_and_validate(path, SCHEMA)
    except jsonschema.ValidationError as e:
        if hint not in (".".join(str(x) for x in e.absolute_path) + " " + e.message):
            raise AssertionError(
                f"{path}: error did not mention '{hint}'. Got: {e.message} (path={list(e.absolute_path)})"
            )
        return
    raise AssertionError(f"{path}: schema validated a config that should have failed")


def test_bad_input_type_rejected():
    _expect_invalid(HERE / "configs" / "bad_input_type.yaml", "input_type")


def test_missing_cohort_name_rejected():
    _expect_invalid(HERE / "configs" / "missing_cohort_name.yaml", "name")


def test_example_cohort_validates():
    load_and_validate(AGNOSTIC / "config" / "cohort.example.yaml", SCHEMA)


def test_template_config_validates():
    load_and_validate(AGNOSTIC / "config" / "config.yaml", SCHEMA)


def _main():
    failures = 0
    for fn in (
        test_bad_input_type_rejected,
        test_missing_cohort_name_rejected,
        test_example_cohort_validates,
        test_template_config_validates,
    ):
        try:
            fn()
            print(f"  PASS  {fn.__name__}")
        except AssertionError as e:
            failures += 1
            print(f"  FAIL  {fn.__name__}: {e}")
    if failures:
        print(f"\n{failures} test(s) failed", file=sys.stderr)
        sys.exit(1)
    print("\nAll tests passed.")


if __name__ == "__main__":
    _main()
