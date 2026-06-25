#!/usr/bin/env python3
"""
Validate the agnostic-pipeline config against config.schema.yaml.

Wraps `jsonschema` so the same code path runs (a) from the Snakefile at parse
time and (b) as a CLI for negative-testing the schema. On failure the message
names the offending key and the constraint it violates, so a user mis-typing
`input_type: "wgz"` gets a one-line error pointing at it.

CLI:
    python validate_config.py path/to/config.yaml path/to/config.schema.yaml

From Snakemake:
    from validate_config import load_and_validate
    config = load_and_validate("config/config.yaml",
                               "config/schema/config.schema.yaml")
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict

import yaml

try:
    import jsonschema  # type: ignore
except ImportError as e:  # pragma: no cover
    print(
        "ERROR: jsonschema is required for config validation. "
        "Install with `pip install jsonschema` inside the active env.",
        file=sys.stderr,
    )
    raise


def _load_yaml(path: str | Path) -> Dict[str, Any]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"{p} does not exist")
    with p.open() as fh:
        return yaml.safe_load(fh)


def load_and_validate(config_path: str | Path,
                      schema_path: str | Path) -> Dict[str, Any]:
    """
    Load `config_path` and validate it against `schema_path`. Returns the
    parsed config on success; raises `jsonschema.ValidationError` (after
    printing a one-line summary) on failure.
    """
    config = _load_yaml(config_path)
    schema = _load_yaml(schema_path)
    validator = jsonschema.Draft7Validator(schema)
    errors = sorted(validator.iter_errors(config), key=lambda e: list(e.absolute_path))
    if errors:
        for e in errors:
            loc = ".".join(str(x) for x in e.absolute_path) or "<root>"
            print(f"CONFIG ERROR at {loc}: {e.message}", file=sys.stderr)
        raise errors[0]
    return config


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2:
        print("Usage: validate_config.py <config.yaml> <config.schema.yaml>", file=sys.stderr)
        return 2
    try:
        load_and_validate(argv[0], argv[1])
    except jsonschema.ValidationError:
        return 1
    print(f"OK: {argv[0]} validates against {argv[1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
