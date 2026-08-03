#!/usr/bin/env python3
"""
Unit + integration tests for scripts/py/contigs_from_fai.py.

Run from the agnostic/ directory:
    python -m pytest tests/test_contigs_from_fai.py -v
or
    python tests/test_contigs_from_fai.py        # falls back to a tiny self-runner
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
AGNOSTIC = HERE.parent
sys.path.insert(0, str(AGNOSTIC / "scripts" / "py"))

from contigs_from_fai import derive_contigs, contig_map  # noqa: E402


# --- 1. Self-contained synthetic .fai --------------------------------------

def _write_synthetic_fai(tmp_path: Path) -> Path:
    fai = tmp_path / "ref.fasta.fai"
    fai.write_text(
        "chr1\t100\t10\t60\t61\n"
        "chr2\t200\t100\t60\t61\n"
        "chrMT\t50\t300\t60\t61\n"
        "PLASTID\t75\t360\t60\t61\n"
    )
    return fai


def test_basic_exclude_drops_mt_and_plastid(tmp_path):
    fai = _write_synthetic_fai(tmp_path)
    got = derive_contigs(fai, exclude=["MT", "PLASTID"])
    assert got == ["chr1", "chr2"], got


def test_no_exclude_keeps_everything(tmp_path):
    fai = _write_synthetic_fai(tmp_path)
    got = derive_contigs(fai, exclude=[])
    assert got == ["chr1", "chr2", "chrMT", "PLASTID"], got


def test_regex_pattern(tmp_path):
    fai = _write_synthetic_fai(tmp_path)
    got = derive_contigs(fai, exclude=[r"^chr[12]$"])
    assert got == ["chrMT", "PLASTID"], got


def test_contig_map_is_1_indexed(tmp_path):
    fai = _write_synthetic_fai(tmp_path)
    contigs = derive_contigs(fai, exclude=["MT", "PLASTID"])
    cmap = contig_map(contigs)
    assert cmap == [("chr1", 1), ("chr2", 2)], cmap


# --- 2. Acceptance test against the tiny-cohort fixture --------------------
# Repointed from V1's real .fai to the self-contained tests/tiny_cohort/
# fixture so the test suite runs standalone (no V1 tree required).

TINY_FAI = AGNOSTIC / "tests" / "tiny_cohort" / "data" / "reference" / "tiny.fasta.fai"


def test_tiny_cohort_fai_two_nuclear_contigs():
    if not TINY_FAI.exists():
        import warnings
        warnings.warn(
            f"tiny cohort .fai not present at {TINY_FAI}; "
            f"regenerate via `python tests/tiny_cohort/generate.py`"
        )
        return
    got = derive_contigs(TINY_FAI, exclude=["MT"])
    assert got == ["chr1", "chr2"], got


# --- Fallback self-runner (no pytest required) ----------------------------

def _main():
    import tempfile

    failures = 0
    tests = [
        test_basic_exclude_drops_mt_and_plastid,
        test_no_exclude_keeps_everything,
        test_regex_pattern,
        test_contig_map_is_1_indexed,
    ]
    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        for t in tests:
            try:
                t(td_path)
                print(f"  PASS  {t.__name__}")
            except AssertionError as e:
                failures += 1
                print(f"  FAIL  {t.__name__}: {e}")

    # Integration test takes no arg
    try:
        test_tiny_cohort_fai_two_nuclear_contigs()
        print("  PASS  test_tiny_cohort_fai_two_nuclear_contigs")
    except AssertionError as e:
        failures += 1
        print(f"  FAIL  test_tiny_cohort_fai_two_nuclear_contigs: {e}")

    if failures:
        print(f"\n{failures} test(s) failed", file=sys.stderr)
        sys.exit(1)
    print("\nAll tests passed.")


if __name__ == "__main__":
    _main()
