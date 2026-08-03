#!/usr/bin/env python3
"""
Unit tests for scripts/py/render_walkthrough.py — the docstring parser and
the config resolver. Runs with or without pytest.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
AGNOSTIC = HERE.parent
sys.path.insert(0, str(AGNOSTIC / "scripts" / "py"))

from render_walkthrough import (  # noqa: E402
    parse_docstring, parse_rules_from_smk, resolve_dotted, load_config,
)


def test_parse_docstring_basic():
    body = """
    One-line summary.

    WHAT: bcftools view -m2 -M2
    WHY:  biallelic downstream
    TUNABLES: qc.min_mac, qc.filter_pass
    OUTPUT: {outputs}/qc/snps.qc.vcf.gz
    TRY:    bump qc.min_mac to 5
    """
    sec = parse_docstring(body)
    assert sec["summary"] == "One-line summary."
    assert sec["what"].startswith("bcftools")
    assert sec["why"].startswith("biallelic")
    assert "qc.min_mac" in sec["tunables"]
    assert "snps.qc.vcf.gz" in sec["output"]
    assert "bump qc.min_mac" in sec["try"]


def test_parse_docstring_multiline_why():
    body = """
    Summary.

    WHAT: x
    WHY:  first sentence.
          second sentence still in WHY.
    OUTPUT: out
    """
    sec = parse_docstring(body)
    assert "first sentence" in sec["why"]
    assert "second sentence" in sec["why"]


def test_parse_rules_from_actual_setup_smk():
    smk = AGNOSTIC / "workflow" / "rules" / "00_setup.smk"
    rules = parse_rules_from_smk(smk)
    names = {r.name for r in rules}
    assert "index_vcf" in names
    assert "extract_vcf_samples" in names
    assert "validate_metadata" in names
    r = next(r for r in rules if r.name == "extract_vcf_samples")
    assert "bcftools query -l" in r.what


def test_resolve_dotted():
    cfg = {"qc": {"min_mac": 2}, "cohort": {"input_type": "wgs"}}
    assert resolve_dotted(cfg, "qc.min_mac") == 2
    assert resolve_dotted(cfg, "cohort.input_type") == "wgs"
    assert resolve_dotted(cfg, "does.not.exist") is None


def test_load_config_layers():
    # config.yaml + cohort.example.yaml merge
    cfg = load_config(AGNOSTIC / "config" / "cohort.example.yaml")
    assert cfg["cohort"]["input_type"] == "wgs"
    assert cfg["structure"]["cluster_labelling"] == "reference"


def _main():
    failures = 0
    tests = [
        test_parse_docstring_basic,
        test_parse_docstring_multiline_why,
        test_parse_rules_from_actual_setup_smk,
        test_resolve_dotted,
        test_load_config_layers,
    ]
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
        except AssertionError as e:
            failures += 1
            print(f"  FAIL  {t.__name__}: {e}")
    if failures:
        print(f"\n{failures} test(s) failed", file=sys.stderr)
        sys.exit(1)
    print("\nAll tests passed.")


if __name__ == "__main__":
    _main()
