"""
Stage 5 introgression tests.

Two layers:

1. `test_r_unit_tests` shells out to tests/R/test_introgression_units.R — the
   component tests (window binning, consensus calling, per-window distance,
   point-in-polygon, both detection rules) live in R because the code under
   test is R. This wrapper just makes them part of `pytest tests/`.

2. The integration tests are the POSITIVE CONTROL: tests/tiny_cohort/ injects
   a known introgression event (six group-A samples carry a group-B donor's
   haplotype over one 10 kb window on chr2) and records it in
   data/introgression_truth.tsv. The pipeline must recover exactly that — the
   right window, for the right samples, and nobody else.

   There is no validated reference output for this stage, so these assertions
   ARE the correctness criterion. See docs/introgression_analysis_spec.md §8.

The integration tests build the tiny cohort's Stage 5 targets with snakemake if
they are not already present (a no-op resolving in seconds when they are), and
skip if the pipeline toolchain is not on PATH.
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
TRUTH    = TINY / "data" / "introgression_truth.tsv"
INTRO    = TINY / "outputs" / "introgression"

PAIR_CALLS = INTRO / "pairs" / "RegionA1__RegionB1.tsv"
FILTERED   = INTRO / "introgressed_windows_filtered.tsv"
AUDIT      = INTRO / "filter_audit.tsv"
HEADLINE   = INTRO / "unique_windows_in_RegionA2_with_freq_and_coords.tsv"
MEMBERSHIP = INTRO / "focal_group_membership.tsv"


def _read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _toolchain_ready() -> bool:
    return all(shutil.which(t) for t in ("snakemake", "Rscript", "bcftools", "plink2"))


@pytest.fixture(scope="module")
def stage5():
    """
    Ensure the tiny cohort's Stage 5 outputs exist, then hand back the truth
    table plus the parsed outputs.
    """
    if not TRUTH.exists():
        pytest.skip(f"fixture truth table missing — run: python {TINY}/generate.py")
    if not _toolchain_ready():
        pytest.skip("pipeline toolchain not on PATH (source envs/activate.sh)")

    targets = [str(p.relative_to(AGNOSTIC)) for p in (FILTERED, AUDIT, HEADLINE)]
    proc = subprocess.run(
        ["snakemake", "--configfile", str(CONFIG.relative_to(AGNOSTIC)),
         "--cores", "2", *targets],
        cwd=AGNOSTIC, capture_output=True, text=True, timeout=1800,
    )
    if proc.returncode != 0:
        pytest.fail("snakemake failed building the Stage 5 targets:\n"
                    + proc.stdout[-4000:] + proc.stderr[-4000:])

    truth = _read_tsv(TRUTH)
    return {
        "injected_samples": sorted({r["sample_id"] for r in truth}),
        "window_bin":       truth[0]["window_bin"],
        "contig":           truth[0]["contig"],
        "start":            int(truth[0]["window_start"]),
        "end":              int(truth[0]["window_end"]),
        "pair":             _read_tsv(PAIR_CALLS),
        "filtered":         _read_tsv(FILTERED),
        "audit":            _read_tsv(AUDIT),
        "headline":         _read_tsv(HEADLINE),
        "membership":       _read_tsv(MEMBERSHIP),
    }


def _expected_window_id(stage5) -> str:
    """The pipeline's window id for the injected window: w<chrom>_<bin>."""
    cmap = TINY / "outputs" / "setup" / "contig_map.tsv"
    codes = {r[0]: r[1] for r in csv.reader(cmap.open(), delimiter="\t")}
    return f"w{codes[stage5['contig']]}_{stage5['window_bin']}"


def test_r_unit_tests():
    """Component unit tests for the distance, contour and detection cores."""
    if not shutil.which("Rscript"):
        pytest.skip("Rscript not on PATH (source envs/activate.sh)")
    script = AGNOSTIC / "tests" / "R" / "test_introgression_units.R"
    proc = subprocess.run(["Rscript", str(script)], cwd=AGNOSTIC,
                          capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stdout[-6000:] + proc.stderr[-4000:]
    assert "0 failed" in proc.stdout


def test_injected_window_is_detected(stage5):
    """Every injected sample is called at the injected window (recall = 1)."""
    window = _expected_window_id(stage5)
    called = {r["SAMPLE"] for r in stage5["pair"] if r["WINDOW"] == window}
    missing = set(stage5["injected_samples"]) - called
    assert not missing, f"injected samples not detected at {window}: {sorted(missing)}"


def test_injected_window_has_no_extra_samples(stage5):
    """No untouched sample is called at the injected window (precision = 1)."""
    window = _expected_window_id(stage5)
    called = {r["SAMPLE"] for r in stage5["pair"] if r["WINDOW"] == window}
    extra = called - set(stage5["injected_samples"])
    assert not extra, f"untouched samples called at {window}: {sorted(extra)}"


def test_distance_rule_recovers_the_injected_window(stage5, tmp_path):
    """
    The deterministic `distance` rule must clear the same recovery gate as the
    density rules: every injected sample called at the injected window, nobody
    else called there.

    The fixture cannot DISCRIMINATE between rules — every rule and threshold it
    has been swept at recovers all six injections (see
    tests/tiny_cohort/outputs/introgression/detection_rule_sweep.tsv). That is
    what the Malay benchmark and the cluster-size test are for. This test only
    asserts `distance` is not broken.
    """
    if not _toolchain_ready():
        pytest.skip("pipeline toolchain not on PATH (source envs/activate.sh)")
    out = tmp_path / "distance_pair.tsv"
    proc = subprocess.run(
        ["Rscript", str(AGNOSTIC / "scripts" / "R" / "introgression_pair.R"),
         "--genotype-table", str(TINY / "outputs" / "ibd" / "combined" / "hmmIBD_input.tsv"),
         "--clusters",       str(TINY / "outputs" / "structure" / "admix_clusters.tsv"),
         "--pair",           "RegionA1__RegionB1",
         "--window-size",    "10000",
         "--min-snps",       "5",
         "--detection-rule", "distance",
         "--contour-level-other", "5e-4",
         "--contour-level-own",   "5e-4",
         "--distance-margin",     "15",
         "--distance-adaptive",   "false",
         "--out", str(out)],
        cwd=AGNOSTIC, capture_output=True, text=True, timeout=900,
    )
    assert proc.returncode == 0, proc.stdout[-4000:] + proc.stderr[-4000:]

    window = _expected_window_id(stage5)
    rows   = _read_tsv(out)
    called = {r["SAMPLE"] for r in rows if r["WINDOW"] == window}
    assert called == set(stage5["injected_samples"]), (
        f"distance rule at {window}: got {sorted(called)}, "
        f"want {stage5['injected_samples']}")
    # The rule fits no density surface, so it reports no density levels.
    assert all(r["LEVEL_OWN"] == "NA" for r in rows), \
        "distance rule emitted density levels — it should not build a surface"


def test_call_direction_is_a_to_b(stage5):
    """The call points the right way: an A-cluster sample looking like B."""
    window = _expected_window_id(stage5)
    directions = {r["DIRECTION"] for r in stage5["pair"] if r["WINDOW"] == window}
    assert directions == {"RegionA1_like_RegionB1"}, directions


def test_detection_levels_are_sane(stage5):
    """
    The injected points must be deeper in the OTHER cluster's cloud than in
    their own — the mechanism the absolute rule tests, made explicit.
    """
    window = _expected_window_id(stage5)
    rows = [r for r in stage5["pair"] if r["WINDOW"] == window]
    assert rows
    for r in rows:
        own   = float(r["LEVEL_OWN"])   if r["LEVEL_OWN"]   != "-Inf" else float("-inf")
        other = float(r["LEVEL_OTHER"]) if r["LEVEL_OTHER"] != "-Inf" else float("-inf")
        assert other > own, f"{r['SAMPLE']}: level_other {other} !> level_own {own}"
        assert float(r["DIST_OTHER"]) < float(r["DIST_OWN"]), \
            f"{r['SAMPLE']}: closer to its own consensus than the other's"


def test_negative_controls_stay_near_zero(stage5):
    """
    Untouched samples generate at most a trickle of calls, and none of them
    survives the cross-dataset filters. The fixture has exactly one real event,
    so anything else in the filtered table is a false positive.
    """
    window = _expected_window_id(stage5)
    raw_fp = [r for r in stage5["pair"] if r["WINDOW"] != window]
    n_true = len([r for r in stage5["pair"] if r["WINDOW"] == window])
    # Pre-filter noise is allowed, but must stay well below the true signal.
    assert len(raw_fp) <= n_true, (
        f"{len(raw_fp)} raw false-positive calls vs {n_true} true calls — "
        "specificity has regressed")
    surviving_fp = [r for r in stage5["filtered"] if r["WINDOW"] != window]
    assert not surviving_fp, (
        "false positives survived the filters: "
        f"{sorted({(r['SAMPLE'], r['WINDOW']) for r in surviving_fp})}")


def test_filtered_table_is_exactly_the_injected_event(stage5):
    """After the cross-dataset filters, the only thing left is the truth."""
    window = _expected_window_id(stage5)
    got = sorted({(r["SAMPLE"], r["WINDOW"]) for r in stage5["filtered"]})
    want = sorted((s, window) for s in stage5["injected_samples"])
    assert got == want


def test_window_coordinates_match_the_injection(stage5):
    """Reported window bounds are the injected window's bounds."""
    rows = [r for r in stage5["filtered"]]
    assert rows
    for r in rows:
        assert r["CONTIG"] == stage5["contig"]
        assert int(r["START"]) == stage5["start"]
        assert int(r["END"])   == stage5["end"]


def test_filter_audit_is_monotonic(stage5):
    """Each filter can only remove calls, never add them."""
    counts = [int(r["n_calls"]) for r in stage5["audit"]]
    assert counts == sorted(counts, reverse=True), counts
    assert counts[0] > 0, "no raw calls at all — the fixture or the pair step broke"


def test_headline_recovers_the_injected_window(stage5):
    """
    The focal-group headline: RegionA2 is both the focal geography and the
    injected sample set, so the headline must be exactly the injected window,
    carried by all six injected samples.
    """
    window = _expected_window_id(stage5)
    assert len(stage5["headline"]) == 1, stage5["headline"]
    row = stage5["headline"][0]
    assert row["WINDOW"] == window
    assert row["CONTIG"] == stage5["contig"]
    assert int(row["START"]) == stage5["start"]
    assert int(row["n_samples"]) == len(stage5["injected_samples"])


def test_headline_membership_split(stage5):
    """The focal/rest split is derived from the cluster, not hardcoded."""
    sides = {}
    for r in stage5["membership"]:
        sides.setdefault(r["side"], set()).add(r["SAMPLE"])
    assert sides["focal"] == set(stage5["injected_samples"])
    assert sides["rest"], "no comparison samples — the split degenerated"
    assert not (sides["focal"] & sides["rest"])
