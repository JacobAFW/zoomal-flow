# 04_ibd.smk — Stage 4 IBD (hmmIBD per-cluster + generic clonal detection)
# --------------------------------------------------------------------------
# input_type note: hmmIBD assumes biallelic diploid-style genotypes, so a
# future microhap path would need a different relatedness estimator. The
# WGS-vs-microhap fork is contained upstream at the Stage-3 seam
# (structure_prep_wgs.smk vs structure_prep_microhap.smk) — Stage 4 reads
# the shared seam output (outputs/structure/cleaned.{bed,bim,fam}) and
# doesn't branch on input_type here.
#
# Cluster list: read from the checkpointed outputs/structure/admix_clusters.tsv
# at DAG-resolution time. No hardcoded ["Mf","Mn","Peninsular"] anywhere.
#
# Input-prep principle: V1's lab-isolates + manually-curated-drops
# machinery does NOT port. Those samples are pre-cleaned out of the input
# (see agnostic/data/INPUT_PREP.md); the Stage-3 `cleaned` bfile is already
# deduplicated. So "per-cluster exclude list" collapses to plain cluster-
# membership keep-list.
# --------------------------------------------------------------------------

import csv
import os
import sys

IBD = config["ibd"]


def _read_cluster_membership(admix_tsv):
    """
    Read `admix_clusters.tsv` (Sample, Component, Proportion, Cluster) into
    a dict {cluster_name: [sample_id, ...]} without pandas.
    """
    memb = {}
    with open(admix_tsv) as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for r in rows:
            sample  = r.get("Sample") or r.get("sample_id") or r.get("sample")
            cluster = r.get("Cluster")
            if sample is None or cluster is None:
                continue
            memb.setdefault(cluster, []).append(sample)
    return memb


def _ibd_clusters(admix_tsv):
    """
    Cluster names eligible for per-cluster hmmIBD: those with at least
    `ibd.min_cluster_n` members. Others are dropped with a logged note the
    first time this is evaluated in a run.
    """
    memb = _read_cluster_membership(admix_tsv)
    min_n = IBD["min_cluster_n"]
    kept, dropped = [], []
    for c, samples in memb.items():
        (kept if len(samples) >= min_n else dropped).append((c, len(samples)))
    if dropped:
        # Snakemake evaluates input functions repeatedly; guard the log.
        if not getattr(_ibd_clusters, "_logged", False):
            for c, n in dropped:
                sys.stderr.write(
                    f"[04_ibd] skipping cluster {c!r} (n={n} < min_cluster_n={min_n})\n"
                )
            _ibd_clusters._logged = True
    return sorted(c for c, _ in kept)


def ibd_cluster_targets(wildcards):
    """
    Snakemake input function: expand per-cluster hmmIBD outputs over the
    cluster list read from the checkpoint. `all` / FINAL_TARGETS consume
    this.
    """
    admix_tsv = checkpoints.assign_clusters.get(**wildcards).output.tsv
    clusters  = _ibd_clusters(admix_tsv)
    return expand(
        f"{PATHS['outputs']}/ibd/{{cluster}}/{{cluster}}.hmm_fract.txt",
        cluster=clusters,
    )


def ibd_fract_specs(wildcards):
    """
    Return the list of `<fract_path>:<cluster>` arg specs the clonal_clusters
    R script expects (one per included cluster).
    """
    admix_tsv = checkpoints.assign_clusters.get(**wildcards).output.tsv
    clusters  = _ibd_clusters(admix_tsv)
    return [
        f"{PATHS['outputs']}/ibd/{c}/{c}.hmm_fract.txt:{c}" for c in clusters
    ]


def ibd_fract_files(wildcards):
    """Just the file paths (for the `input:` block)."""
    admix_tsv = checkpoints.assign_clusters.get(**wildcards).output.tsv
    clusters  = _ibd_clusters(admix_tsv)
    return [f"{PATHS['outputs']}/ibd/{c}/{c}.hmm_fract.txt" for c in clusters]


# --------------------------------------------------------------------------
# Per-cluster IBD
# --------------------------------------------------------------------------

rule cluster_membership:
    """
    Write the cluster's keep-list (FID + IID per line) for PLINK --keep.

    WHAT: read admix_clusters.tsv, filter to this cluster, emit sample IDs.
    WHY:  Replaces V1's cluster_exclude_lists rule + its 4-source union
          (lab isolates + manually-curated dups). Those are gone from
          clean input; membership from admix_clusters is the only
          agnostic input to per-cluster subsetting.
    TUNABLES: (none)
    OUTPUT: {outputs}/ibd/{cluster}/keep.txt
    TRY:    peek at the file — one line per sample, `sample sample`
            format matching plink --keep expectations.
    """
    input:
        clusters = f"{PATHS['outputs']}/structure/admix_clusters.tsv",
    output:
        keep = f"{PATHS['outputs']}/ibd/{{cluster}}/keep.txt",
    log:
        f"{PATHS['logs']}/ibd/cluster_membership_{{cluster}}.log",
    message:
        "[ibd] Cluster membership: {wildcards.cluster}"
    shell:
        r"""
        mkdir -p $(dirname {output.keep})
        awk -F'\t' -v c="{wildcards.cluster}" '
            NR>1 && $4 == c {{ print $1, $1 }}
        ' {input.clusters} > {output.keep} 2> {log}
        echo "Cluster {wildcards.cluster} samples: $(wc -l < {output.keep})" >> {log}
        """


rule cluster_maf_filter:
    """
    Per-cluster bfile + VCF: keep only cluster members, apply within-
    cluster --maf.

    WHAT: plink2 --keep + --maf → make-bed + recode vcf.
    WHY:  Cluster-specific MAF is the point of per-cluster IBD; variants
          polymorphic in one cluster may be fixed in another.
    TUNABLES: ibd.cluster_min_maf
    OUTPUT: {outputs}/ibd/{cluster}/cleaned.{bed,bim,fam,vcf.gz}
    TRY:    bump ibd.cluster_min_maf to 0.05 for a tighter within-cluster
            variant set; usually shrinks hmmIBD input by ~30-40%.
    """
    input:
        bed  = f"{PATHS['outputs']}/structure/cleaned.bed",
        bim  = f"{PATHS['outputs']}/structure/cleaned.bim",
        fam  = f"{PATHS['outputs']}/structure/cleaned.fam",
        keep = rules.cluster_membership.output.keep,
    output:
        bed = f"{PATHS['outputs']}/ibd/{{cluster}}/cleaned.bed",
        bim = f"{PATHS['outputs']}/ibd/{{cluster}}/cleaned.bim",
        fam = f"{PATHS['outputs']}/ibd/{{cluster}}/cleaned.fam",
        vcf = f"{PATHS['outputs']}/ibd/{{cluster}}/cleaned.vcf.gz",
    log:
        f"{PATHS['logs']}/ibd/cluster_maf_filter_{{cluster}}.log",
    params:
        in_prefix  = f"{PATHS['outputs']}/structure/cleaned",
        out_prefix = f"{PATHS['outputs']}/ibd/{{cluster}}/cleaned",
        cluster_maf= IBD["cluster_min_maf"],
    threads: config["compute"]["threads_heavy"]
    message:
        "[ibd] Per-cluster MAF filter: {wildcards.cluster}"
    shell:
        r"""
        plink2 --bfile {params.in_prefix} \
            --allow-extra-chr \
            --keep {input.keep} \
            --maf {params.cluster_maf} \
            --make-bed \
            --threads {threads} \
            --out {params.out_prefix} > {log} 2>&1
        plink2 --bfile {params.out_prefix} \
            --allow-extra-chr \
            --recode vcf bgz \
            --threads {threads} \
            --out {params.out_prefix} >> {log} 2>&1
        echo "Cluster {wildcards.cluster} samples:  $(wc -l < {params.out_prefix}.fam)" >> {log}
        echo "Cluster {wildcards.cluster} variants: $(wc -l < {params.out_prefix}.bim)" >> {log}
        """


rule genotype_table:
    """
    Per-cluster VCF → hmmIBD genotype-table format.

    WHAT: scripts/R/genotype_table.R, using outputs/setup/contig_map.tsv
          to map CHROM strings → integers.
    WHY:  hmmIBD wants integer CHROM. The contig_map is derived once
          in Stage 0 from the reference .fai — the agnostic replacement
          for V1's hardcoded `ordered_PKNH_/_v2` gsub.
    TUNABLES: (none)
    OUTPUT: {outputs}/ibd/{cluster}/hmmIBD_input.tsv
    TRY:    `head -3` the tsv — chrom column should be small ints in
            the same order as contig_map.tsv.
    """
    input:
        vcf  = rules.cluster_maf_filter.output.vcf,
        cmap = f"{PATHS['outputs']}/setup/contig_map.tsv",
    output:
        tsv = f"{PATHS['outputs']}/ibd/{{cluster}}/hmmIBD_input.tsv",
    log:
        f"{PATHS['logs']}/ibd/genotype_table_{{cluster}}.log",
    params:
        script = str(_AGNOSTIC / "scripts" / "R" / "genotype_table.R"),
    message:
        "[ibd] Genotype table: {wildcards.cluster}"
    shell:
        r"""
        Rscript {params.script} {input.vcf} {input.cmap} {output.tsv} \
            > {log} 2>&1
        """


rule run_hmmibd:
    """
    Run hmmIBD pairwise within a cluster.

    WHAT: hmmIBD -i <genotype_table.tsv> -o <prefix>
    WHY:  Produces .hmm_fract.txt (per-pair fract_sites_IBD) which is the
          basis for both connectivity plots and clonal-pair detection.
    TUNABLES: (single-threaded; hmmIBD has no useful tuning knobs for
              cross-sample IBD estimation on this scale)
    OUTPUT: {outputs}/ibd/{cluster}/{cluster}.hmm_fract.txt + .hmm.txt
    TRY:    inspect the hmm_fract.txt — the fract_sites_IBD histogram
            should have a bulk near 0 and a small clonal spike near 1.
    """
    input:
        tsv = rules.genotype_table.output.tsv,
    output:
        fract = f"{PATHS['outputs']}/ibd/{{cluster}}/{{cluster}}.hmm_fract.txt",
        segs  = f"{PATHS['outputs']}/ibd/{{cluster}}/{{cluster}}.hmm.txt",
    log:
        f"{PATHS['logs']}/ibd/run_hmmibd_{{cluster}}.log",
    params:
        prefix = f"{PATHS['outputs']}/ibd/{{cluster}}/{{cluster}}",
    threads: 1
    message:
        "[ibd] hmmIBD: {wildcards.cluster}"
    shell:
        r"""
        hmmIBD -i {input.tsv} -o {params.prefix} > {log} 2>&1
        """


# --------------------------------------------------------------------------
# Combined all-samples genotype table (Stage 5 input)
# --------------------------------------------------------------------------

rule combined_vcf:
    """
    Export the Stage-3 `cleaned` bfile straight to VCF. No exclusion step:
    the cleaned bfile already has dups + lab isolates + manually-curated
    samples removed (see data/INPUT_PREP.md for the Indo cohort's prep).

    WHAT: plink2 --recode vcf bgz on cleaned.{bed,bim,fam}
    WHY:  Stage 5 (introgression) reads the combined genotype table for
          all Stage-3-passing samples.
    TUNABLES: (none)
    OUTPUT: {outputs}/ibd/combined/cleaned.vcf.gz
    TRY:    n = wc -l cleaned.fam should equal the sample count in the
            combined VCF header.
    """
    input:
        bed = f"{PATHS['outputs']}/structure/cleaned.bed",
        bim = f"{PATHS['outputs']}/structure/cleaned.bim",
        fam = f"{PATHS['outputs']}/structure/cleaned.fam",
    output:
        vcf = f"{PATHS['outputs']}/ibd/combined/cleaned.vcf.gz",
    log:
        f"{PATHS['logs']}/ibd/combined_vcf.log",
    params:
        in_prefix  = f"{PATHS['outputs']}/structure/cleaned",
        out_prefix = f"{PATHS['outputs']}/ibd/combined/cleaned",
    threads: config["compute"]["threads_heavy"]
    message:
        "[ibd] Combined (all-samples) VCF export"
    shell:
        r"""
        mkdir -p $(dirname {output.vcf})
        plink2 --bfile {params.in_prefix} \
            --allow-extra-chr \
            --recode vcf bgz \
            --threads {threads} \
            --out {params.out_prefix} > {log} 2>&1
        """


rule combined_genotype_table:
    """
    All-samples hmmIBD-format table (Stage 5 introgression input).

    WHAT: scripts/R/genotype_table.R on the combined VCF using contig_map.
    WHY:  Introgression needs every sample under the same genotype-encoded
          matrix, not per-cluster subsets.
    TUNABLES: (none)
    OUTPUT: {outputs}/ibd/combined/hmmIBD_input.tsv
    TRY:    Stage 5's rules read from this path; touch it to force a
            downstream re-run without redoing the plink export.
    """
    input:
        vcf  = rules.combined_vcf.output.vcf,
        cmap = f"{PATHS['outputs']}/setup/contig_map.tsv",
    output:
        tsv = f"{PATHS['outputs']}/ibd/combined/hmmIBD_input.tsv",
    log:
        f"{PATHS['logs']}/ibd/combined_genotype_table.log",
    params:
        script = str(_AGNOSTIC / "scripts" / "R" / "genotype_table.R"),
    message:
        "[ibd] Combined genotype table (Stage 5 input)"
    shell:
        r"""
        Rscript {params.script} {input.vcf} {input.cmap} {output.tsv} \
            > {log} 2>&1
        """


# --------------------------------------------------------------------------
# Generic clonal-pair detection
# --------------------------------------------------------------------------

def _focal_path(_):
    fc = IBD.get("focal_cluster")
    return f"{PATHS['outputs']}/ibd/focal_{fc}_clones.tsv" if fc else "NULL"


rule clonal_clusters:
    """
    Aggregate every cluster's hmm_fract.txt, flag pairs at or above
    ibd.clonal_ibd_threshold, derive clonal groups as connected components
    of the pair graph. Join to canonical geography/date roles for a tidy
    per-sample table.

    WHAT: scripts/R/clonal_clusters.R over every cluster's fract file.
    WHY:  V1 hardcoded a Peninsular-only pipeline plus a manual 8-sample
          tribble for downstream. The agnostic version computes clonal
          groups from the IBD graph itself; the cohort-specific tribble
          never returns.
    TUNABLES: ibd.clonal_ibd_threshold, ibd.focal_cluster
    OUTPUT: {outputs}/ibd/clonal_clusters.tsv       (always)
            {outputs}/ibd/focal_<name>_clones.tsv    (only if
                                                       ibd.focal_cluster set)
    TRY:    drop ibd.clonal_ibd_threshold to 0.5 to see near-clonal
            structure; the connected-component count will spike.
    """
    input:
        fracts   = ibd_fract_files,
        metadata = rules.validate_metadata.output.tsv,
    output:
        clusters_tsv = f"{PATHS['outputs']}/ibd/clonal_clusters.tsv",
        focal_tsv    = (f"{PATHS['outputs']}/ibd/focal_{IBD['focal_cluster']}_clones.tsv"
                        if IBD.get("focal_cluster") else []),
    log:
        f"{PATHS['logs']}/ibd/clonal_clusters.log",
    params:
        threshold     = IBD["clonal_ibd_threshold"],
        focal_cluster = IBD.get("focal_cluster") or "NULL",
        focal_path    = (f"{PATHS['outputs']}/ibd/focal_{IBD['focal_cluster']}_clones.tsv"
                         if IBD.get("focal_cluster") else "NULL"),
        specs         = ibd_fract_specs,
        script        = str(_AGNOSTIC / "scripts" / "R" / "clonal_clusters.R"),
    message:
        "[ibd] Clonal groups from connected components"
    shell:
        r"""
        mkdir -p $(dirname {output.clusters_tsv})
        Rscript {params.script} \
            {params.threshold} \
            {input.metadata} \
            {output.clusters_tsv} \
            {params.focal_path} \
            {params.focal_cluster} \
            {params.specs} \
            > {log} 2>&1
        """
