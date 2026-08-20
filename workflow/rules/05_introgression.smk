# 05_introgression.smk — Stage 5 introgression (pairwise density clouds)
# --------------------------------------------------------------------------
# Method: docs/introgression_analysis_spec.md. In one paragraph — for each
# cluster PAIR (Kx, Ky) we compute every sample-window's percent-mismatch
# distance to both clusters' consensus alleles, draw a 2D kernel-density cloud
# per cluster in (distance-to-Kx, distance-to-Ky) space, and flag a window as
# introgressed when a sample lands in the OTHER cluster's cloud rather than its
# own. V1 did this with three cluster names welded into the code; the pairwise
# form keeps the identical density method and scales to any K.
#
# DETECTION RULE — three are available (introgression.detection_rule):
# `absolute` (V1, shipped default) and `relative` both read a fitted 2D density
# surface; `distance` decides from the raw per-window distances alone and fits
# nothing, so its calls are reproducible across cohort shifts. The density
# rules' window-level calls are known to be cohort-sensitive — see the "Method
# status" block in docs/introgression_analysis_spec.md before interpreting them.
#
# VALIDATION STATUS — this stage was never validated against an HPC/reference
# output (V1's own script header records that reference outputs were
# unavailable). Correctness rests on the synthetic positive controls in
# tests/tiny_cohort/, not on matching V1. Treat real-cohort numbers as a method
# result for review, not a validated figure.
#
# GUARDRAIL — `introgression.pairs: "all"` is C(K,2) pairs: K=3→3, K=6→15,
# K=8→28, K=10→45, each a full density pass. Two costs, both real:
#   compute        — N pairwise comparisons, linear in N, quadratic in K;
#   interpretation — N contrasts carry a multiple-comparison burden, and not
#                    every cluster pair is biologically meaningful.
# A warning fires above `introgression.pair_warn_threshold`; set an explicit
# `pairs` list to focus on the comparisons you actually mean.
#
# Cluster list + pair list come from the checkpointed
# outputs/structure/admix_clusters.tsv at DAG-resolution time — no cluster
# names anywhere in this file.
#
# input_type note: like Stage 4, this stage consumes the shared Stage-3 seam
# (the combined genotype table) and does not branch on input_type.
# --------------------------------------------------------------------------

import itertools
import os

INTRO = config.get("introgression", {})

INTRO_DIR   = f"{PATHS['outputs']}/introgression"
INTRO_FOCAL = INTRO.get("focal_group")
INTRO_GFF   = INTRO.get("gff")

# The GFF is optional AND static. Declare it as a rule input only when it is
# actually present, so a configured-but-missing annotation degrades to a
# logged skip inside the R script instead of a Snakemake MissingInputException.
_INTRO_GFF_INPUT = [INTRO_GFF] if (INTRO_GFF and os.path.exists(INTRO_GFF)) else []
if INTRO_GFF and not _INTRO_GFF_INPUT:
    sys.stderr.write(
        f"[05_introgression] introgression.gff set to {INTRO_GFF!r} but the file is "
        f"absent — the gene-family filter will be skipped\n"
    )


def _intro_clusters(admix_tsv):
    """Clusters large enough to take part in a pair (introgression.min_cluster_n)."""
    memb  = _read_cluster_membership(admix_tsv)      # from 04_ibd.smk
    min_n = INTRO.get("min_cluster_n", 5)
    kept, dropped = [], []
    for c, samples in memb.items():
        (kept if len(samples) >= min_n else dropped).append((c, len(samples)))
    if dropped and not getattr(_intro_clusters, "_logged", False):
        for c, n in dropped:
            sys.stderr.write(
                f"[05_introgression] skipping cluster {c!r} (n={n} < min_cluster_n={min_n})\n"
            )
        _intro_clusters._logged = True
    return sorted(c for c, _ in kept)


def _intro_pairs(admix_tsv):
    """
    The (Kx, Ky) pairs to run. `introgression.pairs` is either "all" — every
    C(K,2) combination of eligible clusters — or an explicit [[Kx, Ky], …]
    list, whose entries are checked against the eligible clusters and skipped
    (with a warning) if either side is missing or too small.
    """
    clusters = _intro_clusters(admix_tsv)
    spec = INTRO.get("pairs", "all")

    if spec == "all" or spec is None:
        pairs = list(itertools.combinations(clusters, 2))
    else:
        pairs = []
        for entry in spec:
            if len(entry) != 2:
                sys.exit(f"ERROR: introgression.pairs entry {entry!r} is not a [Kx, Ky] pair")
            kx, ky = entry
            bad = [k for k in (kx, ky) if k not in clusters]
            if bad:
                sys.stderr.write(
                    f"[05_introgression] skipping pair {kx}__{ky}: cluster(s) "
                    f"{bad} absent or below min_cluster_n\n"
                )
                continue
            pairs.append(tuple(sorted((kx, ky))))
        pairs = sorted(set(pairs))

    warn_at = INTRO.get("pair_warn_threshold", 15)
    if len(pairs) > warn_at and not getattr(_intro_pairs, "_warned", False):
        sys.stderr.write(
            "\n"
            "########################################################################\n"
            f"[05_introgression] WARNING: {len(pairs)} pairwise comparisons queued "
            f"(> pair_warn_threshold = {warn_at}), from K = {len(clusters)} clusters.\n"
            "\n"
            f"  COMPUTE: {len(pairs)} full 2D-density passes will run, each over the\n"
            "           two clusters' whole genotype table. All-pairs cost is\n"
            "           quadratic in K, so this grows fast (K=8 -> 28, K=10 -> 45).\n"
            "\n"
            "  INTERPRETATION: these are multiple contrasts. The more pairs you run,\n"
            "           the more windows clear any fixed threshold by chance, and NOT\n"
            "           every cluster pair is a biologically meaningful comparison.\n"
            "           Consider setting an explicit `introgression.pairs` list for the\n"
            "           contrasts you actually intend to interpret.\n"
            "########################################################################\n\n"
        )
        _intro_pairs._warned = True
    return pairs


def _pair_id(kx, ky):
    """Wildcard-safe pair encoding: double underscore keeps single-underscore labels intact."""
    return f"{kx}__{ky}"


def introgression_pair_targets(wildcards):
    """FINAL_TARGETS / aggregate input: one call file per pair."""
    admix_tsv = checkpoints.assign_clusters.get(**wildcards).output.tsv
    return [f"{INTRO_DIR}/pairs/{_pair_id(kx, ky)}.tsv"
            for kx, ky in _intro_pairs(admix_tsv)]


wildcard_constraints:
    pair = r"[^/]+",


# --------------------------------------------------------------------------
# Per-pair detection (the fan-out)
# --------------------------------------------------------------------------

rule introgression_pair:
    """
    Detect introgressed windows for one cluster pair, in both directions.

    WHAT: scripts/R/introgression_pair.R — subsets the combined genotype table
          to the pair's samples, computes both clusters' consensus alleles and
          each sample-window's percent-mismatch distance to both, builds the 2D
          kernel-density cloud per cluster, and calls introgression via
          scripts/R/introgression_detect.R (Kx-looks-like-Ky AND
          Ky-looks-like-Kx).
    WHY:  V1's detection was intrinsically three-cluster — two fixed distance
          axes with Mn/Mf/Peninsular welded in by name — and could not scale as
          ADMIXTURE K climbs. One clean 2-cluster density plot per pair keeps
          the identical method, drops every cluster-name hardcode, and works at
          any K. V1's three-way "ambiguous between Mf and Mn" case comes back
          at aggregation as "flagged in both the Kx-Mf and Kx-Mn pairs".
          Snakemake fans the {pair} wildcard across cores; there is no R-level
          parallelism.
    TUNABLES: introgression.window_size_bp, introgression.min_snps_per_window,
              introgression.detection_rule, introgression.contour_level_other,
              introgression.contour_level_own, introgression.distance_margin,
              introgression.distance_adaptive,
              introgression.distance_adaptive_quantile, introgression.pairs,
              introgression.min_cluster_n
    OUTPUT: {outputs}/introgression/pairs/{{pair}}.tsv
    TRY:    switch introgression.detection_rule to "distance" and diff the call
            count for one pair. `distance` never fits a density surface — it
            asks only whether the window matches the OTHER cluster's consensus
            at least `distance_margin` percentage points better than its own —
            so unlike `absolute`/`relative` its calls cannot move when cohort
            composition changes. Set `distance_adaptive: true` to derive that
            margin from the clusters' own spread instead of a fixed number.
    """
    input:
        gt_table = rules.combined_genotype_table.output.tsv,
        clusters = f"{PATHS['outputs']}/structure/admix_clusters.tsv",
    output:
        calls = f"{INTRO_DIR}/pairs/{{pair}}.tsv",
    log:
        f"{PATHS['logs']}/introgression/pair_{{pair}}.log",
    params:
        window_size = INTRO.get("window_size_bp", 10000),
        min_snps    = INTRO.get("min_snps_per_window", 5),
        rule_name   = INTRO.get("detection_rule", "absolute"),
        lvl_other   = INTRO.get("contour_level_other", 5.0e-4),
        lvl_own     = INTRO.get("contour_level_own", 5.0e-4),
        dist_margin = INTRO.get("distance_margin", 15),
        dist_adapt  = str(bool(INTRO.get("distance_adaptive", False))).lower(),
        dist_q      = INTRO.get("distance_adaptive_quantile", 0.9),
        script      = str(_AGNOSTIC / "scripts" / "R" / "introgression_pair.R"),
    threads: 1
    message:
        "[introgression] Pair detection: {wildcards.pair}"
    shell:
        r"""
        mkdir -p $(dirname {output.calls}) $(dirname {log})
        Rscript {params.script} \
            --genotype-table      {input.gt_table} \
            --clusters            {input.clusters} \
            --pair                {wildcards.pair} \
            --window-size         {params.window_size} \
            --min-snps            {params.min_snps} \
            --detection-rule      {params.rule_name} \
            --contour-level-other {params.lvl_other} \
            --contour-level-own   {params.lvl_own} \
            --distance-margin     {params.dist_margin} \
            --distance-adaptive   {params.dist_adapt} \
            --distance-adaptive-quantile {params.dist_q} \
            --out                 {output.calls} \
            > {log} 2>&1
        """


# --------------------------------------------------------------------------
# Aggregation (the fan-in)
# --------------------------------------------------------------------------

rule introgression_aggregate:
    """
    Stitch every pair's calls and apply the four cross-dataset filters.

    WHAT: scripts/R/introgression_aggregate.R — concatenates the per-pair call
          files, then filters in V1's order: dataset low-n (> min_samples_per_
          window samples), per-cluster minimum percentage, gene-family mask
          from the GFF, hypervariable (window called in more than one cluster).
          Writes the filtered call table, a per-step audit, and the per-cluster
          / per-sample / per-chromosome / per-geography summaries.
    WHY:  the per-pair step is deliberately local. Every judgement that needs
          the whole dataset — how many samples share a window, whether it is
          hypervariable across clusters, whether it sits in a hypervariable
          gene family — has to happen once, after the fan-out. The gene-family
          filter is config keywords against the GFF rather than V1's literal
          SICA|KIR grep, and skips itself with a logged note when no GFF is set.
    TUNABLES: introgression.min_samples_per_window,
              introgression.per_cluster_min_pct,
              introgression.per_cluster_min_samples,
              introgression.gene_family_filters, introgression.gff,
              introgression.window_size_bp
    OUTPUT: {outputs}/introgression/introgressed_windows_filtered.tsv,
            {outputs}/introgression/filter_audit.tsv,
            {outputs}/introgression/window_sample_counts_raw.tsv,
            {outputs}/introgression/windows_by_cluster.tsv,
            {outputs}/introgression/windows_across_chrom.tsv,
            {outputs}/introgression/average_windows_for_clusters.tsv,
            {outputs}/introgression/intro_per_sample_summary.tsv
    TRY:    read filter_audit.tsv top-to-bottom — it shows how many calls,
            windows and samples each filter removed. If the hypervariable step
            is eating most of your signal, the clusters are probably sharing
            variable regions rather than exchanging haplotypes.
    """
    input:
        calls    = introgression_pair_targets,
        clusters = f"{PATHS['outputs']}/structure/admix_clusters.tsv",
        metadata = rules.validate_metadata.output.tsv,
        cmap     = f"{PATHS['outputs']}/setup/contig_map.tsv",
        # _FAI is the config-derived reference index (Snakefile: REF["fasta"] + ".fai").
        fai      = _FAI,
        gff      = _INTRO_GFF_INPUT,
    output:
        filtered   = f"{INTRO_DIR}/introgressed_windows_filtered.tsv",
        audit      = f"{INTRO_DIR}/filter_audit.tsv",
        raw_counts = f"{INTRO_DIR}/window_sample_counts_raw.tsv",
        by_cluster = f"{INTRO_DIR}/windows_by_cluster.tsv",
        by_chrom   = f"{INTRO_DIR}/windows_across_chrom.tsv",
        avg_clust  = f"{INTRO_DIR}/average_windows_for_clusters.tsv",
        per_sample = f"{INTRO_DIR}/intro_per_sample_summary.tsv",
        by_geo     = ([f"{INTRO_DIR}/introgression_by_geography.tsv"]
                      if ROLES.get("geography") else []),
    log:
        f"{PATHS['logs']}/introgression/aggregate.log",
    params:
        out_dir      = INTRO_DIR,
        window_size  = INTRO.get("window_size_bp", 10000),
        min_samp_win = INTRO.get("min_samples_per_window", 2),
        pct_min      = INTRO.get("per_cluster_min_pct", 0.05),
        pct_min_n    = INTRO.get("per_cluster_min_samples", 0),
        gff          = INTRO_GFF or "NULL",
        gene_kw      = ",".join(INTRO.get("gene_family_filters") or []) or "NULL",
        script       = str(_AGNOSTIC / "scripts" / "R" / "introgression_aggregate.R"),
    message:
        "[introgression] Aggregate + cross-dataset filters"
    shell:
        r"""
        mkdir -p {params.out_dir} $(dirname {log})
        Rscript {params.script} \
            --clusters               {input.clusters} \
            --metadata               {input.metadata} \
            --contig-map             {input.cmap} \
            --fai                    {input.fai} \
            --window-size            {params.window_size} \
            --min-samples-per-window {params.min_samp_win} \
            --per-cluster-min-pct    {params.pct_min} \
            --per-cluster-min-samples {params.pct_min_n} \
            --gff                    "{params.gff}" \
            --gene-family-filters    "{params.gene_kw}" \
            --out-dir                {params.out_dir} \
            {input.calls} \
            > {log} 2>&1
        """


# --------------------------------------------------------------------------
# Headline + figures
# --------------------------------------------------------------------------

rule plot_introgression_shoulder:
    """
    Windows ranked by how many samples carry them, before any filtering.

    WHAT: scripts/R/plot_introgression_shoulder.R over
          window_sample_counts_raw.tsv, with the configured threshold drawn in.
    WHY:  this is the diagnostic that justifies the dataset low-n cut. The
          curve's shoulder is where "several samples share this window" turns
          into "one sample has it"; the dashed line shows whether
          min_samples_per_window sits on that shoulder or somewhere arbitrary.
    TUNABLES: introgression.min_samples_per_window
    OUTPUT: {reports}/figures/introgression_shoulder.png/.svg
    TRY:    if the dashed line sits well left of the shoulder you are keeping
            singleton noise; move min_samples_per_window right and re-run the
            aggregate step (the per-pair calls are cached, so it is seconds).
    """
    input:
        counts = rules.introgression_aggregate.output.raw_counts,
    output:
        png = f"{PATHS['reports']}/figures/introgression_shoulder.png",
        svg = f"{PATHS['reports']}/figures/introgression_shoulder.svg",
    log:
        f"{PATHS['logs']}/introgression/plot_shoulder.log",
    params:
        threshold = INTRO.get("min_samples_per_window", 2),
        script    = str(_AGNOSTIC / "scripts" / "R" / "plot_introgression_shoulder.R"),
    message:
        "[introgression:figures] Shoulder plot"
    shell:
        r"""
        mkdir -p $(dirname {output.png}) $(dirname {log})
        Rscript {params.script} \
            {input.counts} {params.threshold} {output.png} {output.svg} \
            > {log} 2>&1
        """


if INTRO_FOCAL:

    rule introgression_headline:
        """
        Windows introgressed uniquely in the focal group (the headline result).

        WHAT: scripts/R/introgression_headline.R — splits the cluster that the
              focal group's samples belong to into "focal" vs "rest", and finds
              the windows called on exactly one side.
        WHY:  this answers the question the stage exists for ("what
              introgression does this group carry that its cluster-mates do
              not?"). V1 hardcoded `State %in% c("Aceh","Peninsular")`; here
              `focal_group` is matched against the `geography` role and the
              comparison set is derived from the cluster assignment, so any
              cohort can ask the same question of any group.
        TUNABLES: introgression.focal_group, metadata.roles.geography
        OUTPUT: {outputs}/introgression/unique_windows_in_<focal>_with_freq_and_coords.tsv,
                {outputs}/introgression/unique_windows_per_chrom_<focal>_vs_rest.tsv,
                {outputs}/introgression/focal_group_membership.tsv
        TRY:    point introgression.focal_group at a different geography value
                and re-run — the comparison cluster is re-derived, so you get
                that group's unique windows without editing any code.
        """
        input:
            calls    = rules.introgression_aggregate.output.filtered,
            clusters = f"{PATHS['outputs']}/structure/admix_clusters.tsv",
            metadata = rules.validate_metadata.output.tsv,
        output:
            unique    = f"{INTRO_DIR}/unique_windows_in_{INTRO_FOCAL}_with_freq_and_coords.tsv",
            per_chrom = f"{INTRO_DIR}/unique_windows_per_chrom_{INTRO_FOCAL}_vs_rest.tsv",
        log:
            f"{PATHS['logs']}/introgression/headline_{INTRO_FOCAL}.log",
        params:
            out_dir = INTRO_DIR,
            focal   = INTRO_FOCAL,
            script  = str(_AGNOSTIC / "scripts" / "R" / "introgression_headline.R"),
        message:
            f"[introgression] Headline: windows unique to '{INTRO_FOCAL}'"
        shell:
            r"""
            mkdir -p {params.out_dir} $(dirname {log})
            Rscript {params.script} \
                --calls         {input.calls} \
                --clusters      {input.clusters} \
                --metadata      {input.metadata} \
                --focal-group   "{params.focal}" \
                --out-dir       {params.out_dir} \
                --out-unique    {output.unique} \
                --out-per-chrom {output.per_chrom} \
                > {log} 2>&1
            """

    rule plot_introgression_focal:
        """
        Genomic map of the focal-unique windows (the headline figure).

        WHAT: scripts/R/plot_introgression_focal.R — lollipop per window,
              x = position, y = focal samples carrying it, one panel per contig
              that has a hit.
        WHY:  the headline TSV answers "how many"; this answers "where". V1's
              version hardcoded the label "Aceh" and a 14-chromosome factor;
              contigs here come from the aggregate step's contig_map join.
        TUNABLES: introgression.focal_group
        OUTPUT: {reports}/figures/introgression_focal_<focal>.png/.svg
        TRY:    clustered hits on one contig are worth a look in the GFF — a
                run of adjacent focal-unique windows is more likely one
                introgressed haplotype block than several independent events.
        """
        input:
            unique = rules.introgression_headline.output.unique,
        output:
            png = f"{PATHS['reports']}/figures/introgression_focal_{INTRO_FOCAL}.png",
            svg = f"{PATHS['reports']}/figures/introgression_focal_{INTRO_FOCAL}.svg",
        log:
            f"{PATHS['logs']}/introgression/plot_focal_{INTRO_FOCAL}.log",
        params:
            focal  = INTRO_FOCAL,
            script = str(_AGNOSTIC / "scripts" / "R" / "plot_introgression_focal.R"),
        message:
            f"[introgression:figures] Focal-unique windows: '{INTRO_FOCAL}'"
        shell:
            r"""
            mkdir -p $(dirname {output.png}) $(dirname {log})
            Rscript {params.script} \
                {input.unique} "{params.focal}" {output.png} {output.svg} \
                > {log} 2>&1
            """
