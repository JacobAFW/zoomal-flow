#!/usr/bin/env python3
"""
Generate the tiny-cohort test fixture.

Produces (all under tests/tiny_cohort/):
  data/reference/tiny.fasta      + tiny.fasta.fai
  data/reference/regions_to_mask.list
  data/vcf/tiny.vcf.gz           + .csi (via bcftools; created below)
  data/metadata/samples.tsv
  config.yaml

  data/introgression_truth.tsv   (Stage 5 positive-control ground truth)

The VCF injects two ancestry groups (A and B) whose allele frequencies
diverge on a set of "structure SNPs" so ADMIXTURE at K=2 recovers them.
Every sample carries a country + geography + date value; a few controls
are included so the `controls.exclude_patterns` filter has real work.

Design targets:
  - ~40 samples: 18 "groupA" + 18 "groupB" + 4 controls (ctrl_1..4)
  - 3 contigs (chr1 length 200k, chr2 100k, chrMT 5k — MT excluded via
    reference.exclude_contigs regex)
  - ~630 biallelic SNPs — chr1: 250 structure + 250 neutral; chr2: 50
    structure + 50 neutral outside the injected window, plus 15 + 15 inside
    it (MT contigs are empty since they'd be excluded anyway)
  - FILTER=PASS for all records; realistic depth/GQ tags
  - Deterministic seed for reproducibility

Stage 5 positive control (the introgression ground truth)
--------------------------------------------------------
There is no validated reference output for the introgression stage, so
correctness is anchored on a manufactured event: six group-A samples have one
10 kb window on chr2 replaced by a single group-B donor's genotypes — a
transferred haplotype. Those six samples are exactly group A's `RegionA2`
geography, which is also the config's `introgression.focal_group`, so one
fixture drives both the detection test and the headline test.

Two fixture-design points that the detection method forces, both learned the
hard way:

  1. Injection copies a DONOR SAMPLE, not the group-B consensus. A
     consensus-overwritten sample sits 0% away from B's consensus, whereas
     real B samples sit ~10-20% away — so a consensus-injected point lands off
     the edge of B's density cloud instead of inside it, and no detection rule
     can see it. Copying a donor puts the injected sample exactly where that
     donor already is: in the cloud's core. It is also what introgression
     physically is.
  2. Divergence is UNIFORM across windows (see STRUCTURE_AF_A/B). The clouds
     are estimated over all windows at once, so a window whose divergence is
     far from the genome-wide norm forms its own sparse lobe that no fixed
     density threshold reaches.

The injected samples, the donor, and the window coordinates are written to
`data/introgression_truth.tsv` — tests/test_introgression.py asserts the
pipeline recovers exactly that, and scripts/R/introgression_rule_sweep.R uses
it to score detection rules against false positives.

Deps: numpy, bcftools on PATH. Reference FASTA is written directly; the
VCF is emitted as plain text then bgzipped + indexed via bcftools.

Usage:
    python tests/tiny_cohort/generate.py
"""

from __future__ import annotations

import gzip
import random
import shutil
import subprocess
import textwrap
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

RNG_SEED = 20260803
random.seed(RNG_SEED)
np.random.seed(RNG_SEED)


CONTIGS = [
    ("chr1",  200_000),
    ("chr2",  100_000),
    ("chrMT",   5_000),   # excluded via reference.exclude_contigs
]
N_GROUP_A = 25
N_GROUP_B = 25
N_CTRL    = 4

N_STRUCTURE_SNPS_CHR1 = 250
N_NEUTRAL_SNPS_CHR1   = 250
N_STRUCTURE_SNPS_CHR2 = 50
N_NEUTRAL_SNPS_CHR2   = 50

# Structure-SNP frequencies. Fixed and strongly divergent (not the wide Beta
# draws an earlier revision used) so that EVERY window has roughly the same
# between-group divergence. Stage 5 needs that: its density clouds are built
# over all windows at once, so windows with wildly different divergence smear
# each cluster's cloud across the distance plane and the method loses both
# sensitivity and specificity. Uniform divergence gives two tight, well
# separated clouds — which is also what a real pair of diverged clusters looks
# like. 0.9/0.1 rather than 1.0/0.0 keeps realistic within-group variation.
STRUCTURE_AF_A = 0.90
STRUCTURE_AF_B = 0.10

BASES = ["A", "C", "G", "T"]

# --- Stage 5 positive control -------------------------------------------
# Group A splits into two geographies. RegionA2 is the focal group AND the set
# of samples whose INTRO_WINDOW genotypes are replaced by a group-B donor's.
# RegionA1 stays the majority, so `cluster_labelling: "auto"` still names group
# A's cluster "RegionA1".
N_REGION_A2   = 6                    # focal + injected samples (of N_GROUP_A)
INTRO_CONTIG  = "chr2"
INTRO_WINDOW_SIZE  = 10_000          # must match introgression.window_size_bp
INTRO_WINDOW_START = 40_000          # window [40000, 50000) -> midpoint 45000
N_INTRO_STRUCTURE_SNPS = 15          # a normal structure/neutral mix, just denser,
N_INTRO_NEUTRAL_SNPS   = 15          # so the window always clears min_snps_per_window


def _rand_base(exclude: str) -> str:
    b = random.choice(BASES)
    while b == exclude:
        b = random.choice(BASES)
    return b


def write_reference():
    DATA_REF = DATA / "reference"
    DATA_REF.mkdir(parents=True, exist_ok=True)
    fasta = DATA_REF / "tiny.fasta"
    fai   = DATA_REF / "tiny.fasta.fai"
    with fasta.open("w") as fh:
        for name, length in CONTIGS:
            seq = "".join(random.choices(BASES, k=length))
            fh.write(f">{name}\n")
            for i in range(0, length, 60):
                fh.write(seq[i:i + 60] + "\n")
    # Build the .fai manually — 5 columns: name, length, offset, linebases, linewidth
    lines = []
    offset = 0
    for name, length in CONTIGS:
        header_line = f">{name}\n"
        offset += len(header_line)
        lines.append(f"{name}\t{length}\t{offset}\t60\t61")
        n_lines = (length + 59) // 60
        offset += length + n_lines
    fai.write_text("\n".join(lines) + "\n")

    # A tiny regions-to-mask list (a single 500 bp span on chr1 — QC seam).
    mask = DATA_REF / "regions_to_mask.list"
    mask.write_text("chr1:100000-100500\n")


def sample_names():
    a = [f"sampA_{i:02d}" for i in range(1, N_GROUP_A + 1)]
    b = [f"sampB_{i:02d}" for i in range(1, N_GROUP_B + 1)]
    ctrl = [f"ctrl_{i}" for i in range(1, N_CTRL + 1)]
    return a, b, ctrl


def build_snp_positions(n_snps: int, contig_len: int, avoid: set = None) -> list:
    positions = set()
    while len(positions) < n_snps:
        pos = random.randint(1000, contig_len - 1000)
        if avoid and pos in avoid:
            continue
        positions.add(pos)
    return sorted(positions)


def build_snp_positions_in_range(n_snps: int, lo: int, hi: int,
                                 avoid: set = None) -> list:
    """Positions inside an explicit [lo, hi] span (the injected window)."""
    positions = set()
    while len(positions) < n_snps:
        pos = random.randint(lo, hi)
        if avoid and pos in avoid:
            continue
        positions.add(pos)
    return sorted(positions)


def _call_with_gt(gt: str) -> str:
    """Wrap a GT in a fresh GT:DP:AD field, so injected samples don't share depths."""
    if gt == "./.":
        return "./.:.:.,."
    depth = random.randint(20, 60)
    noise = min(random.randint(0, 2), depth - 1)
    if gt in ("1/1", "1|1"):
        return f"{gt}:{depth}:{noise},{depth - noise}"
    return f"{gt}:{depth}:{depth - noise},{noise}"


def write_introgression_truth(injected: list, donor: str, n_snps: int):
    """
    Ground truth for the Stage 5 positive control: which samples carry the
    transferred block, from which donor, and the window's coordinates.
    `window_bin` is the midpoint the pipeline's `window_bin()` assigns
    (floor(pos/ws)*ws + ws/2).
    """
    DATA.mkdir(parents=True, exist_ok=True)
    bin_mid = (INTRO_WINDOW_START // INTRO_WINDOW_SIZE) * INTRO_WINDOW_SIZE \
        + INTRO_WINDOW_SIZE // 2
    rows = ["sample_id\tdonor\tcontig\twindow_bin\twindow_start\twindow_end\t"
            "n_injected_snps"]
    for s in injected:
        rows.append(f"{s}\t{donor}\t{INTRO_CONTIG}\t{bin_mid}\t{INTRO_WINDOW_START}\t"
                    f"{INTRO_WINDOW_START + INTRO_WINDOW_SIZE - 1}\t{n_snps}")
    (DATA / "introgression_truth.tsv").write_text("\n".join(rows) + "\n")
    print(f"[tiny_cohort] injected introgression: {len(injected)} sample(s) "
          f"carry {donor}'s haplotype over "
          f"{INTRO_CONTIG}:{INTRO_WINDOW_START}-"
          f"{INTRO_WINDOW_START + INTRO_WINDOW_SIZE - 1} ({n_snps} SNPs)")


def gt_from_freq(af_A: float, af_B: float, group: str) -> str:
    """
    Draw a call in `GT:DP:AD` FORMAT for a P. knowlesi-style sample. Use
    hom-alt (1/1) with probability af, hom-ref (0/0) otherwise, plus a
    small missing-rate. AD is needed by moimix (bafMatrix) at Stage 2.
    """
    af = af_A if group == "A" else af_B
    r = random.random()
    if r < 0.02:
        return "./.:.:.,."
    depth = random.randint(20, 60)
    if random.random() < af:
        # Hom-alt: all reads support ALT, with a couple of REF noise reads.
        noise = min(random.randint(0, 2), depth - 1)
        return f"1/1:{depth}:{noise},{depth - noise}"
    # Hom-ref: reversed.
    noise = min(random.randint(0, 2), depth - 1)
    return f"0/0:{depth}:{depth - noise},{noise}"


def introgressed_samples(a: list) -> list:
    """The group-A samples whose INTRO_WINDOW is overwritten with B's consensus."""
    return a[-N_REGION_A2:]


def geography_of(sample: str, a: list) -> str:
    """Group A's geography: RegionA2 for the injected tail, RegionA1 otherwise."""
    return "RegionA2" if sample in introgressed_samples(a) else "RegionA1"


def write_samples_tsv(a: list, b: list, ctrl: list):
    DATA_META = DATA / "metadata"
    DATA_META.mkdir(parents=True, exist_ok=True)
    rows = ["sample_id\tcountry\tgeography\thost\tdate"]
    dates_a = ["2024-03-05", "2024-04-12", "2024-05-19", "2024-06-27",
               "2024-07-14", "2024-08-02", "2024-09-08", "2024-10-15",
               "2024-11-21"]
    dates_b = ["2025-01-11", "2025-02-04", "2025-03-15", "2025-04-06",
               "2025-05-22", "2025-06-11", "2025-07-29"]
    for i, s in enumerate(a):
        rows.append(f"{s}\tCountryA\t{geography_of(s, a)}\tHuman\t{dates_a[i % len(dates_a)]}")
    for i, s in enumerate(b):
        rows.append(f"{s}\tCountryB\tRegionB1\tHuman\t{dates_b[i % len(dates_b)]}")
    for s in ctrl:
        rows.append(f"{s}\tControl\tNA\tNA\tNA")
    (DATA_META / "samples.tsv").write_text("\n".join(rows) + "\n")


def write_vcf(a: list, b: list, ctrl: list) -> Path:
    DATA_VCF = DATA / "vcf"
    DATA_VCF.mkdir(parents=True, exist_ok=True)
    samples = a + b + ctrl
    vcf_path = DATA_VCF / "tiny.vcf"

    with vcf_path.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        for name, length in CONTIGS:
            fh.write(f"##contig=<ID={name},length={length}>\n")
        fh.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        fh.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n")
        fh.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths (REF,ALT)\">\n")
        cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                "FILTER", "INFO", "FORMAT"] + samples
        fh.write("\t".join(cols) + "\n")

        injected  = introgressed_samples(a)
        donor_idx = 0                       # b[0] donates the introgressed block
        n_intro_snps = 0

        def emit(contig, pos, kind, in_intro_window=False):
            """
            Write one record. `kind` is "structure" (divergent between groups)
            or "neutral" (same AF in both). Inside the injected window, the
            injected samples' calls are replaced by the donor B sample's call —
            a transferred haplotype, which is what introgression physically is.
            """
            nonlocal n_intro_snps
            ref = random.choice(BASES)
            alt = _rand_base(ref)
            if kind == "structure":
                af_A, af_B = STRUCTURE_AF_A, STRUCTURE_AF_B
            else:
                af = float(np.clip(np.random.beta(2, 2), 0.05, 0.95))
                af_A = af_B = af
            a_calls = [gt_from_freq(af_A, af_B, "A") for _ in a]
            b_calls = [gt_from_freq(af_A, af_B, "B") for _ in b]
            if in_intro_window:
                donor_gt = b_calls[donor_idx].split(":")[0]
                for i, s in enumerate(a):
                    if s in injected:
                        a_calls[i] = _call_with_gt(donor_gt)
                n_intro_snps += 1
            fh.write("\t".join([
                contig, str(pos), ".", ref, alt, "60", "PASS", ".", "GT:DP:AD"
            ] + a_calls + b_calls
              + [gt_from_freq(af_A, af_B, "A") for _ in ctrl]) + "\n")

        # chr1: 250 structure SNPs (divergent between groups) + 250 neutral
        chr1_len = dict(CONTIGS)["chr1"]
        structure_pos = build_snp_positions(N_STRUCTURE_SNPS_CHR1, chr1_len)
        neutral_pos   = build_snp_positions(N_NEUTRAL_SNPS_CHR1, chr1_len,
                                            avoid=set(structure_pos))
        for pos in structure_pos:
            emit("chr1", pos, "structure")
        for pos in neutral_pos:
            emit("chr1", pos, "neutral")

        # chr2: structure + neutral, placed OUTSIDE the injected window (which
        # gets its own guaranteed-dense block below, so the positive control
        # never depends on how many SNPs happened to land there).
        chr2_len = dict(CONTIGS)["chr2"]
        intro_span = set(range(INTRO_WINDOW_START,
                               INTRO_WINDOW_START + INTRO_WINDOW_SIZE))
        chr2_struct = build_snp_positions(N_STRUCTURE_SNPS_CHR2, chr2_len,
                                          avoid=intro_span)
        chr2_neutral = build_snp_positions(N_NEUTRAL_SNPS_CHR2, chr2_len,
                                           avoid=intro_span | set(chr2_struct))
        for pos in chr2_struct:
            emit("chr2", pos, "structure")
        for pos in chr2_neutral:
            emit("chr2", pos, "neutral")

        # ---- Stage 5 positive control -------------------------------------
        # The injected window itself: an ordinary mix of structure + neutral
        # SNPs (so the window's distance profile is TYPICAL of the genome, and
        # the injected samples land in the core of group B's density cloud
        # rather than off its edge), with the injected samples carrying the
        # donor's genotypes throughout.
        intro_struct = build_snp_positions_in_range(
            N_INTRO_STRUCTURE_SNPS,
            INTRO_WINDOW_START + 100,
            INTRO_WINDOW_START + INTRO_WINDOW_SIZE - 100)
        intro_neutral = build_snp_positions_in_range(
            N_INTRO_NEUTRAL_SNPS,
            INTRO_WINDOW_START + 100,
            INTRO_WINDOW_START + INTRO_WINDOW_SIZE - 100,
            avoid=set(intro_struct))
        for pos in intro_struct:
            emit(INTRO_CONTIG, pos, "structure", in_intro_window=True)
        for pos in intro_neutral:
            emit(INTRO_CONTIG, pos, "neutral", in_intro_window=True)

    write_introgression_truth(injected, b[donor_idx], n_intro_snps)

    # Sort by chrom, pos then bgzip + index via bcftools.
    bcftools = shutil.which("bcftools")
    if not bcftools:
        raise SystemExit("bcftools not found — activate the env first.")
    sorted_vcf = DATA_VCF / "tiny.sorted.vcf"
    subprocess.check_call([bcftools, "sort", "-Ov", "-o", str(sorted_vcf), str(vcf_path)])
    vcf_gz = DATA_VCF / "tiny.vcf.gz"
    subprocess.check_call([bcftools, "view", "-Oz", "-o", str(vcf_gz), str(sorted_vcf)])
    subprocess.check_call([bcftools, "index", "-c", "-f", str(vcf_gz)])
    vcf_path.unlink()
    sorted_vcf.unlink()
    return vcf_gz


def write_config():
    cfg = textwrap.dedent(f"""\
        # Tiny synthetic cohort — smoke-test config for Stages 0-4.
        cohort:
          name: "tiny_cohort"
          input_type: "wgs"
          vcf: "tests/tiny_cohort/data/vcf/tiny.vcf.gz"

        reference:
          fasta: "tests/tiny_cohort/data/reference/tiny.fasta"
          exclude_contigs: ["MT"]
          mask_regions: "tests/tiny_cohort/data/reference/regions_to_mask.list"

        metadata:
          table: "tests/tiny_cohort/data/metadata/samples.tsv"
          gis: null
          roles:
            sample_id:    "sample_id"
            group:        null
            geography:    "geography"
            country:      "country"
            host:         "host"
            date:         "date"
            case_control: null

        controls:
          exclude_patterns: ["^ctrl"]

        qc:
          min_mac: 2
          filter_pass: true

        moi:
          min_maf: 0.05
          fws_polyclonal_cutoff: 0.95
          # Cutoff dropped from V1's 0.95 default: the tiny cohort's
          # simulated AD values give an intrinsic ~0.86 Fws floor, and
          # 0.95 would exclude the whole cohort. 0.5 keeps everyone in
          # for the Stage-3+ smoke test.
          fws_exclusion_cutoff:  0.5
          seed: 2023

        structure:
          min_maf: 0.05
          max_sample_missing: 0.20
          max_variant_missing: 0.20
          admixture_k_min: 2
          admixture_k_max: 3
          admixture_cv_folds: 5
          admixture_k: null
          cluster_labelling: "auto"     # tag by majority geography role
          duplicate_id_pattern: null

        ibd:
          cluster_min_maf: 0.05
          clonal_ibd_threshold: 0.95
          min_cluster_n: 5
          focal_cluster: null

        # Stage 5: the two clusters (RegionA1 / RegionB1) give exactly one
        # pair. focal_group is RegionA2 — group A's minority geography, which
        # is also the set of samples carrying the injected window, so the
        # headline is testable against data/introgression_truth.tsv.
        introgression:
          window_size_bp: {INTRO_WINDOW_SIZE}
          min_snps_per_window: 5
          min_cluster_n: 5
          min_samples_per_window: 2
          per_cluster_min_pct: 0.05
          pairs: "all"
          pair_warn_threshold: 15
          detection_rule: "absolute"
          contour_level_other: 5.0e-4
          contour_level_own:   5.0e-4
          focal_group: "RegionA2"
          gene_family_filters: []      # no annotation in the fixture
          gff: null

        selection:
          models: []

        compute:
          threads_heavy: 2
          threads_light: 1

        paths:
          outputs: "tests/tiny_cohort/outputs"
          logs:    "tests/tiny_cohort/logs"
          reports: "tests/tiny_cohort/reports"
    """)
    (HERE / "config.yaml").write_text(cfg)


def main():
    print(f"[tiny_cohort] generating fixture under {HERE.relative_to(HERE.parent.parent)}/")
    write_reference()
    a, b, ctrl = sample_names()
    write_samples_tsv(a, b, ctrl)
    vcf = write_vcf(a, b, ctrl)
    write_config()
    print(f"[tiny_cohort] wrote {vcf.relative_to(HERE.parent.parent)}")
    print("[tiny_cohort] done.")


if __name__ == "__main__":
    main()
