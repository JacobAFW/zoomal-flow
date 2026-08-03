#!/usr/bin/env python3
"""
Generate the tiny-cohort test fixture.

Produces (all under tests/tiny_cohort/):
  data/reference/tiny.fasta      + tiny.fasta.fai
  data/reference/regions_to_mask.list
  data/vcf/tiny.vcf.gz           + .csi (via bcftools; created below)
  data/metadata/samples.tsv
  config.yaml

The VCF injects two ancestry groups (A and B) whose allele frequencies
diverge on a set of "structure SNPs" so ADMIXTURE at K=2 recovers them.
Every sample carries a country + geography + date value; a few controls
are included so the `controls.exclude_patterns` filter has real work.

Design targets:
  - ~40 samples: 18 "groupA" + 18 "groupB" + 4 controls (ctrl_1..4)
  - 3 contigs (chr1 length 200k, chr2 100k, chrMT 5k — MT excluded via
    reference.exclude_contigs regex)
  - ~600 biallelic SNPs (250 structure + 250 neutral on chr1 and 100 on
    chr2; MT contigs mostly empty since they'd be excluded anyway)
  - FILTER=PASS for all records; realistic depth/GQ tags
  - Deterministic seed for reproducibility

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
N_NEUTRAL_SNPS_CHR2   = 100

BASES = ["A", "C", "G", "T"]


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
        rows.append(f"{s}\tCountryA\tRegionA1\tHuman\t{dates_a[i % len(dates_a)]}")
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

        # chr1: 250 structure SNPs (divergent AF between groups) + 250 neutral
        chr1_len = dict(CONTIGS)["chr1"]
        structure_pos = build_snp_positions(N_STRUCTURE_SNPS_CHR1, chr1_len)
        neutral_pos   = build_snp_positions(N_NEUTRAL_SNPS_CHR1, chr1_len,
                                            avoid=set(structure_pos))
        for pos in structure_pos:
            ref = random.choice(BASES)
            alt = _rand_base(ref)
            af_A = np.clip(np.random.beta(4, 2), 0.1, 0.9)     # A skew high
            af_B = np.clip(np.random.beta(2, 4), 0.1, 0.9)     # B skew low
            fh.write("\t".join([
                "chr1", str(pos), ".", ref, alt, "60", "PASS", ".", "GT:DP:AD"
            ] + [gt_from_freq(af_A, af_B, "A") for _ in a]
              + [gt_from_freq(af_A, af_B, "B") for _ in b]
              + [gt_from_freq(0.5, 0.5, "A") for _ in ctrl]) + "\n")
        for pos in neutral_pos:
            ref = random.choice(BASES)
            alt = _rand_base(ref)
            af = np.clip(np.random.beta(2, 2), 0.05, 0.95)
            fh.write("\t".join([
                "chr1", str(pos), ".", ref, alt, "60", "PASS", ".", "GT:DP:AD"
            ] + [gt_from_freq(af, af, "A") for _ in a]
              + [gt_from_freq(af, af, "B") for _ in b]
              + [gt_from_freq(af, af, "A") for _ in ctrl]) + "\n")

        # chr2: 100 neutral SNPs
        chr2_len = dict(CONTIGS)["chr2"]
        chr2_pos = build_snp_positions(N_NEUTRAL_SNPS_CHR2, chr2_len)
        for pos in chr2_pos:
            ref = random.choice(BASES)
            alt = _rand_base(ref)
            af = np.clip(np.random.beta(2, 2), 0.05, 0.95)
            fh.write("\t".join([
                "chr2", str(pos), ".", ref, alt, "60", "PASS", ".", "GT:DP:AD"
            ] + [gt_from_freq(af, af, "A") for _ in a]
              + [gt_from_freq(af, af, "B") for _ in b]
              + [gt_from_freq(af, af, "A") for _ in ctrl]) + "\n")

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
