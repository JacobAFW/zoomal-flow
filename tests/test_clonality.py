"""
Clonality / de-clonalization tests (docs/clonality.md).

The fixture carries a KNOWN clonal block: `tests/tiny_cohort/generate.py`
appends N_CLONAL genotype-identical replicates of one group-B sample and
records them in data/clonality_truth.tsv. hmmIBD scores them at
fract_sites_IBD = 1.0, so Stage 4 must group them and de-clonalization must
collapse them to exactly one representative.

Three things are checked, matching the three the increment claims:

  (a) the `unique` sample set drops all but one member of the clonal block;
  (b) the `full` arm still behaves exactly as before the {sampleset} axis
      existed — the pre-existing Stage-5 assertions in test_introgression.py
      cover this directly, and here we check the structural invariant that
      `full` excludes nobody;
  (c) a frequency result inflated by the clonal block visibly shifts in the
      `unique` variant, while a result that owes nothing to clonality does not.

Skips cleanly if the pipeline toolchain is not on PATH.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path

import pytest

AGNOSTIC = Path(__file__).resolve().parent.parent
TINY     = AGNOSTIC / "tests" / "tiny_cohort"
CONFIG   = TINY / "config.yaml"
OUT      = TINY / "outputs"
TRUTH    = TINY / "data" / "clonality_truth.tsv"

CLONALITY   = OUT / "clonality"
UNIQUE_LIST = CLONALITY / "unique_genotypes.txt"
ALL_LIST    = CLONALITY / "all_genotypes.txt"
EXCL_UNIQUE = CLONALITY / "exclude_unique.txt"
EXCL_FULL   = CLONALITY / "exclude_full.txt"
AUDIT       = CLONALITY / "declonalization_audit.tsv"

CMP_MOI     = OUT / "moi" / "declonalization_comparison.tsv"
CMP_STRUCT  = OUT / "structure" / "declonalization_comparison.tsv"


def _read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _lines(path: Path) -> list[str]:
    return [ln for ln in path.read_text().splitlines() if ln.strip()]


def _toolchain_ready() -> bool:
    return all(shutil.which(t) for t in ("snakemake", "Rscript", "bcftools", "plink2"))


@pytest.fixture(scope="module")
def clonality():
    """Build the clonality artifacts + the comparison tables for the fixture."""
    if not TRUTH.exists():
        pytest.skip(f"fixture truth table missing — run: python {TINY}/generate.py")
    if not _toolchain_ready():
        pytest.skip("pipeline toolchain not on PATH (source envs/activate.sh)")

    targets = [str(p.relative_to(AGNOSTIC))
               for p in (UNIQUE_LIST, ALL_LIST, EXCL_UNIQUE, EXCL_FULL, AUDIT,
                         CMP_MOI, CMP_STRUCT)]
    proc = subprocess.run(
        ["snakemake", "--configfile", str(CONFIG.relative_to(AGNOSTIC)),
         "--cores", "2", "--", *targets],
        cwd=AGNOSTIC, capture_output=True, text=True, timeout=3600,
    )
    if proc.returncode != 0:
        pytest.fail("snakemake failed building the clonality targets:\n"
                    + proc.stdout[-4000:] + proc.stderr[-4000:])

    truth = _read_tsv(TRUTH)
    return {
        "block":      sorted({r["sample_id"] for r in truth}),
        "source":     truth[0]["clonal_source"],
        "all":        _lines(ALL_LIST),
        "unique":     _lines(UNIQUE_LIST),
        "excl_full":  _lines(EXCL_FULL),
        "excl_uniq":  _lines(EXCL_UNIQUE),
        "audit":      _read_tsv(AUDIT),
        "cmp_moi":    _read_tsv(CMP_MOI),
        "cmp_struct": _read_tsv(CMP_STRUCT),
    }


# --------------------------------------------------------------------------
# (a) the unique set collapses the clonal block
# --------------------------------------------------------------------------

def test_clonal_block_collapses_to_one_representative(clonality):
    """Exactly one member of the known clonal block survives de-clonalization."""
    block = set(clonality["block"])
    survivors = block & set(clonality["unique"])
    assert len(survivors) == 1, (
        f"expected 1 representative of the clonal block, got {sorted(survivors)}")


def test_clonal_block_members_are_the_ones_dropped(clonality):
    """
    Everything dropped comes from a clonal group — de-clonalization must never
    remove a sample that is not a clonal replicate.
    """
    dropped = set(clonality["excl_uniq"])
    assert dropped, "nothing was dropped, but the fixture has a clonal block"
    grouped = set()
    for row in clonality["audit"]:
        grouped |= set(filter(None, row["dropped"].split(";")))
        grouped.add(row["representative"])
    assert dropped <= grouped, f"dropped samples outside any clonal group: {dropped - grouped}"


def test_unique_is_full_minus_the_dropped(clonality):
    """The two keep-lists and the drop-list are mutually consistent."""
    assert set(clonality["unique"]) == set(clonality["all"]) - set(clonality["excl_uniq"])
    assert len(clonality["unique"]) == len(clonality["all"]) - len(clonality["excl_uniq"])


def test_non_clonal_samples_all_survive(clonality):
    """Every sample outside a clonal group is kept, untouched."""
    block = set(clonality["block"])
    non_clonal = set(clonality["all"]) - block
    assert non_clonal <= set(clonality["unique"])


def test_audit_records_the_block(clonality):
    """The audit names the representative and the dropped members per group."""
    assert clonality["audit"], "audit table is empty"
    block = set(clonality["block"])
    hit = [r for r in clonality["audit"]
           if block & ({r["representative"]} | set(filter(None, r["dropped"].split(";"))))]
    assert hit, "the known clonal block does not appear in the audit"
    row = hit[0]
    assert int(row["n_members"]) == len(block)
    assert int(row["n_dropped"]) == len(block) - 1
    assert row["straddles_clusters"] in ("FALSE", "TRUE")


def test_representative_choice_is_deterministic(clonality):
    """
    Re-running the rule picks the same representative. The rule is lowest
    missingness then alphabetical, so this must not depend on row order or RNG.
    """
    before = set(_lines(UNIQUE_LIST))
    proc = subprocess.run(
        ["snakemake", "--configfile", str(CONFIG.relative_to(AGNOSTIC)),
         "--cores", "1", "--forcerun", "unique_genotypes", "--",
         str(UNIQUE_LIST.relative_to(AGNOSTIC))],
        cwd=AGNOSTIC, capture_output=True, text=True, timeout=1800)
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]
    assert set(_lines(UNIQUE_LIST)) == before


# --------------------------------------------------------------------------
# (b) the full arm is a genuine no-op
# --------------------------------------------------------------------------

def test_full_sample_set_excludes_nobody(clonality):
    """
    `full`'s exclusion list is EMPTY. This is the structural guarantee that the
    full arm reproduces the pre-clonality pipeline: every stage's filter step
    receives an empty drop-list and is therefore unchanged.
    """
    assert clonality["excl_full"] == []
    assert EXCL_FULL.exists() and EXCL_FULL.stat().st_size == 0


def test_full_keep_list_is_every_sample(clonality):
    """Nothing is silently missing from the full set."""
    assert set(clonality["block"]) <= set(clonality["all"])
    assert len(clonality["all"]) > len(clonality["unique"])


# --------------------------------------------------------------------------
# (c) an inflated frequency result moves; an unrelated one does not
# --------------------------------------------------------------------------

def _metric(rows, name):
    for r in rows:
        if r["metric"] == name:
            return r
    return None


def test_comparison_reports_the_collapse_size(clonality):
    """Every stage's comparison leads with the size of the intervention."""
    for rows in (clonality["cmp_moi"], clonality["cmp_struct"]):
        removed = _metric(rows, "n_samples_removed")
        assert removed is not None, "n_samples_removed missing from the comparison"
        assert float(removed["unique"]) == len(clonality["excl_uniq"])


def test_inflated_group_count_shifts(clonality):
    """
    The clonal block is all one geography, so that group's sample count MUST
    fall in the de-clonalized arm — this is the pseudo-replication the whole
    increment exists to expose.
    """
    rows = clonality["cmp_moi"]
    n_rows = [r for r in rows if r["metric"].startswith("n[geography=")]
    assert n_rows, f"no per-geography n metrics: {[r['metric'] for r in rows]}"
    moved = [r for r in n_rows if r["delta"] not in ("", "NA") and float(r["delta"]) < 0]
    assert moved, ("no geography group shrank under de-clonalization, but the "
                   f"fixture's clonal block is all one geography: {n_rows}")


def test_unaffected_group_does_not_shift(clonality):
    """
    A geography with no clonal members must be identical in both arms. This is
    the specificity half: de-clonalization must not perturb what it should not
    touch.
    """
    rows = [r for r in clonality["cmp_moi"] if r["metric"].startswith("n[geography=")]
    unchanged = [r for r in rows if r["delta"] not in ("", "NA") and float(r["delta"]) == 0]
    assert unchanged, ("every geography moved — de-clonalization should only "
                       f"touch groups containing clonal replicates: {rows}")


def test_comparison_covers_both_arms(clonality):
    """The comparison has a value for each arm, not just one."""
    for rows in (clonality["cmp_moi"], clonality["cmp_struct"]):
        assert rows, "comparison table is empty"
        assert all("full" in r and "unique" in r for r in rows)
        n = _metric(rows, "n_samples_analysed")
        assert n is not None and float(n["full"]) > float(n["unique"])
