# 03_structure.smk — Stage 3 population structure, common rules + seam
# --------------------------------------------------------------------------
# SEAM CONTRACT (DESIGN §4, seam 3):
#   Both structure_prep_wgs.smk and structure_prep_microhap.smk MUST emit
#     {outputs}/structure/cleaned.{bed,bim,fam}    (final filtered set)
#     {outputs}/structure/cleaned.ld.{bed,bim,fam} (LD-pruned, for ADMIXTURE)
#   Everything below in this file (ADMIXTURE, PCA, NJ-tree distance,
#   assign_clusters, optional gis_join) reads those and does NOT branch on
#   input_type — the fork is contained in the prep seam.
#
# Common rules:
#   - admixture_run            per-K ADMIXTURE with CV
#   - admixture_cv_table       collect CV errors into one TSV
#   - select_best_k            write best_k.txt (auto from CV min, or override)
#   - admixture_cv_plot        CV-vs-K diagnostic figure
#   - assign_clusters          label-mode-driven (numbered/auto/reference)
#   - pca / pca_variance       PLINK 2 --pca + tidy variance table
#   - distance_matrix          PLINK 1.9 --distance square (for the NJ tree)
#   - gis_join (optional)      join cluster assignments to a GIS file if one
#                              is configured under metadata.gis
# --------------------------------------------------------------------------

import os

STRUCTURE     = config["structure"]
K_RANGE       = list(range(STRUCTURE["admixture_k_min"],
                           STRUCTURE["admixture_k_max"] + 1))
CV_FOLDS      = STRUCTURE["admixture_cv_folds"]
LABEL_MODE    = STRUCTURE["cluster_labelling"]
K_OVERRIDE    = STRUCTURE.get("admixture_k")
DUP_PATTERN   = STRUCTURE.get("duplicate_id_pattern")    # may be None
GIS_FILE      = META.get("gis")                          # optional; metadata.gis


# --------------------------------------------------------------------------
# input_type seam — pick the structure-prep implementation
# --------------------------------------------------------------------------
_STRUCT_INPUT_TYPE = COHORT["input_type"]
if _STRUCT_INPUT_TYPE == "wgs":
    include: "structure_prep_wgs.smk"
elif _STRUCT_INPUT_TYPE == "microhap":
    include: "structure_prep_microhap.smk"
else:
    raise ValueError(
        f"Unknown cohort.input_type={_STRUCT_INPUT_TYPE!r}; expected 'wgs' or 'microhap'."
    )


# --------------------------------------------------------------------------
# ADMIXTURE
# --------------------------------------------------------------------------

rule admixture_run:
    """
    Run ADMIXTURE for a single K with --cv cross-validation. Operates on
    the LD-pruned bfile from the prep seam. K values run in parallel via
    Snakemake's wildcard expansion.

    WHAT: admixture --cv N cleaned.bed K  (staged in a per-stage dir so .Q
          and .P land under outputs/structure/admixture/)
    WHY:  ADMIXTURE is the slide-7 ancestry-bar method. CV error vs K is
          the standard model-selection diagnostic.
    TUNABLES: structure.admixture_k_min/max, structure.admixture_cv_folds
    OUTPUT: {outputs}/structure/admixture/cleaned.{K}.{Q,P}
            + {logs}/structure/admixture_K{K}.log
    TRY:    bump admixture_k_max to 12 to scan further if the CV-vs-K curve
            hasn't visibly bottomed out by K=10 (rare; usually flatlines).
    """
    input:
        bed = f"{PATHS['outputs']}/structure/cleaned.ld.bed",
        bim = f"{PATHS['outputs']}/structure/cleaned.ld.bim",
        fam = f"{PATHS['outputs']}/structure/cleaned.ld.fam",
    output:
        Q   = f"{PATHS['outputs']}/structure/admixture/cleaned.{{K}}.Q",
        P   = f"{PATHS['outputs']}/structure/admixture/cleaned.{{K}}.P",
        log = f"{PATHS['logs']}/structure/admixture_K{{K}}.log",
    threads: 2
    params:
        admix_dir = f"{PATHS['outputs']}/structure/admixture",
        ld_prefix = f"{PATHS['outputs']}/structure/cleaned.ld",
        cv_folds  = CV_FOLDS,
    message:
        "[structure] ADMIXTURE K={wildcards.K}"
    shell:
        r"""
        mkdir -p {params.admix_dir}
        for ext in bed bim fam; do
            cp -f {params.ld_prefix}.$ext {params.admix_dir}/cleaned.$ext
        done
        ( cd {params.admix_dir} && \
          admixture --cv={params.cv_folds} -j{threads} cleaned.bed {wildcards.K} ) \
            > {output.log} 2>&1
        """


rule admixture_cv_table:
    """
    Collect CV error from every K's log into one tidy TSV. The pattern in
    ADMIXTURE's log is `CV error (K=N): <err>` — awk parses it out.

    WHAT: awk extract over the admixture_K*.log files
    WHY:  The CV curve is the model-selection input + the slide-7 panel.
    TUNABLES: (none — derived from structure.admixture_k_min/max)
    OUTPUT: {outputs}/structure/admixture/cv_error.tsv  (K, CV_error)
    TRY:    grep `CV error` admixture_K*.log to spot-check parsing.
    """
    input:
        logs = expand(f"{PATHS['logs']}/structure/admixture_K{{K}}.log", K=K_RANGE),
    output:
        tsv = f"{PATHS['outputs']}/structure/admixture/cv_error.tsv",
    log:
        f"{PATHS['logs']}/structure/admixture_cv_table.log",
    message:
        "[structure] Building ADMIXTURE CV-error table"
    shell:
        r"""
        echo -e 'K\tCV_error' > {output.tsv}
        for f in {input.logs}; do
            awk -v f="$f" '
              /CV error \(K=/ {{
                k = $0; sub(/.*K=/, "", k); sub(/\):.*/, "", k);
                e = $NF;
                printf "%s\t%s\n", k, e
              }}
            ' $f
        done | sort -k1,1n >> {output.tsv} 2> {log}
        echo "CV rows: $(($(wc -l < {output.tsv}) - 1))" >> {log}
        """


rule select_best_k:
    """
    Pick the best K — minimum CV error from cv_error.tsv unless overridden
    by structure.admixture_k. Writes a one-line file consumed by every
    downstream rule that needs to know which K to read.

    WHAT: awk min over CV_error column, OR echo of the override.
    WHY:  Per DESIGN §7 decision 1: best-K is the CV minimum, with a manual
          override knob. Centralising the choice in one file keeps every
          downstream rule consistent.
    TUNABLES: structure.admixture_k  (null = auto, int = override)
    OUTPUT: {outputs}/structure/best_k.txt
    TRY:    set structure.admixture_k to an integer one step away from the
            auto pick — re-run from assign_clusters down and inspect whether
            the extra component tracks a real geographic / host split.
    """
    input:
        tsv = rules.admixture_cv_table.output.tsv,
    output:
        txt = f"{PATHS['outputs']}/structure/best_k.txt",
    log:
        f"{PATHS['logs']}/structure/select_best_k.log",
    params:
        override = "" if K_OVERRIDE is None else str(K_OVERRIDE),
    message:
        "[structure] Selecting best K"
    shell:
        r"""
        if [ -n "{params.override}" ]; then
            echo "{params.override}" > {output.txt}
            echo "[best_k] using override K = {params.override}" > {log}
        else
            # Pick the K with the smallest CV_error
            awk 'NR>1 {{ if (min=="" || $2<min) {{ min=$2; k=$1 }} }} END {{ print k }}' \
                {input.tsv} > {output.txt} 2> {log}
            echo "[best_k] auto-picked K = $(cat {output.txt}) (min CV error)" >> {log}
        fi
        """


rule admixture_cv_plot:
    """
    CV-error-vs-K diagnostic plot (the slide-7 small panel).

    WHAT: scripts/R/plot_admixture_cv.R
    WHY:  Lets a user read off the best K visually and notice if the curve
          is flat (multiple K's nearly equivalent — adds a wrinkle to the
          interpretation that the single min-K choice doesn't show).
    TUNABLES: (none in the rule; cosmetic in the R script)
    OUTPUT: {reports}/figures/admixture_cv.png (+ .svg)
    TRY:    eyeball the plot — if K=2 and K=3 are within 0.005 of each
            other, try running assign_clusters at K=2 manually and compare.
    """
    input:
        tsv = rules.admixture_cv_table.output.tsv,
        bk  = rules.select_best_k.output.txt,
    output:
        png = f"{PATHS['reports']}/figures/admixture_cv.png",
        svg = f"{PATHS['reports']}/figures/admixture_cv.svg",
    log:
        f"{PATHS['logs']}/structure/admixture_cv_plot.log",
    params:
        script = str(_AGNOSTIC / "scripts" / "R" / "plot_admixture_cv.R"),
    message:
        "[structure] Plotting ADMIXTURE CV error"
    shell:
        r"""
        mkdir -p $(dirname {output.png})
        Rscript {params.script} {input.tsv} {input.bk} {output.png} {output.svg} > {log} 2>&1
        """


checkpoint assign_clusters:
    """
    Read the best-K .Q, assign each sample to its dominant ancestry
    component, label components per structure.cluster_labelling.

    Declared as a `checkpoint` (not a plain rule) because the labels it
    emits drive Stage 4's per-cluster IBD wildcard expansion — the set of
    cluster names isn't known until this file exists. Snakemake re-
    evaluates the DAG once the checkpoint's output is on disk.

    WHAT: scripts/R/assign_clusters.R — label-mode-driven (numbered |
          auto | reference). Reads `best_k.txt` so it never hardcodes a K.
    WHY:  The core lift of the agnostic refactor (DESIGN §5 item 2). V1
          hardcoded K=3 and literal `Mn/Mf/Peninsular` labels — neither
          generalises. The agnostic version:
            numbered   — components → Cluster_1..K (zero-metadata default).
            auto       — components labelled by majority geography role.
            reference  — components labelled by majority group role
                         (this is V1's behaviour, applied via the role map).
    TUNABLES: structure.cluster_labelling, structure.admixture_k
    OUTPUT: {outputs}/structure/admix_clusters.tsv
            (cols: sample_id, Component, Proportion, Cluster)
    TRY:    flip cluster_labelling between numbered/auto/reference on the
            same .Q — the component memberships are identical, only the
            human-readable Cluster column changes. Confirms the lift is
            label-only, never re-runs ADMIXTURE.
    """
    input:
        bk       = rules.select_best_k.output.txt,
        Q_all    = expand(f"{PATHS['outputs']}/structure/admixture/cleaned.{{K}}.Q",
                          K=K_RANGE),
        fam      = f"{PATHS['outputs']}/structure/cleaned.fam",
        metadata = rules.validate_metadata.output.tsv,
    output:
        tsv = f"{PATHS['outputs']}/structure/admix_clusters.tsv",
    log:
        f"{PATHS['logs']}/structure/assign_clusters.log",
    params:
        admix_dir = f"{PATHS['outputs']}/structure/admixture",
        mode      = LABEL_MODE,
        script    = str(_AGNOSTIC / "scripts" / "R" / "assign_clusters.R"),
    message:
        "[structure] Assigning clusters (mode={params.mode})"
    shell:
        r"""
        K=$(cat {input.bk})
        Q="{params.admix_dir}/cleaned.${{K}}.Q"
        Rscript {params.script} \
            "$Q" \
            {input.fam} \
            {input.metadata} \
            {params.mode} \
            {output.tsv} \
            > {log} 2>&1
        """


# --------------------------------------------------------------------------
# PCA + distance (numeric outputs that Stage 3b figures will consume)
# --------------------------------------------------------------------------

rule pca:
    """
    Compute PCA from the cleaned (non-LD-pruned) bfile.

    WHAT: plink2 --pca on cleaned.{bed,bim,fam}
    WHY:  PC1/PC2 separation is the slide-9 panel. Run on the cleaned set
          rather than the LD-pruned one to match V1.
    TUNABLES: (none; default 10 PCs)
    OUTPUT: {outputs}/structure/Pk.eigenvec, Pk.eigenval
    TRY:    compare PC1 variance against the by-region split — if PC1 is
            country-only and PC2 hits region, your structure is hierarchical.
    """
    input:
        bed = f"{PATHS['outputs']}/structure/cleaned.bed",
        bim = f"{PATHS['outputs']}/structure/cleaned.bim",
        fam = f"{PATHS['outputs']}/structure/cleaned.fam",
    output:
        eigenvec = f"{PATHS['outputs']}/structure/Pk.eigenvec",
        eigenval = f"{PATHS['outputs']}/structure/Pk.eigenval",
    log:
        f"{PATHS['logs']}/structure/pca.log",
    threads: config["compute"]["threads_heavy"]
    params:
        in_prefix  = f"{PATHS['outputs']}/structure/cleaned",
        out_prefix = f"{PATHS['outputs']}/structure/Pk",
    message:
        "[structure] PLINK --pca"
    shell:
        r"""
        plink2 --bfile {params.in_prefix} \
            --allow-extra-chr \
            --pca \
            --threads {threads} \
            --out {params.out_prefix} > {log} 2>&1
        """


rule pca_variance:
    """
    Convert .eigenval (one variance per line) into a tidy two-column TSV
    (PC, variance_percent) for the Quarto report and Stage-3b axis labels.

    WHAT: awk sum-and-normalise
    WHY:  Eigenvalues alone are unitless; percent-variance is what plots
          and tables display.
    TUNABLES: (none)
    OUTPUT: {outputs}/structure/pca_variance.tsv
    TRY:    sum the first three rows — typically explains the bulk of
            between-region structure on this scale of cohort.
    """
    input:
        eigenval = rules.pca.output.eigenval,
    output:
        tsv = f"{PATHS['outputs']}/structure/pca_variance.tsv",
    log:
        f"{PATHS['logs']}/structure/pca_variance.log",
    message:
        "[structure] Tidying PCA variance percentages"
    shell:
        r"""
        awk 'BEGIN{{OFS="\t"; print "PC","variance_percent"}}
             {{tot+=$1; v[NR]=$1}}
             END{{ for (i=1;i<=NR;i++) printf "%d\t%.4f\n", i, 100*v[i]/tot }}
        ' {input.eigenval} > {output.tsv} 2> {log}
        """


rule distance_matrix:
    """
    Pairwise IBS-distance matrix used for the NJ tree in Stage 3b.

    WHAT: plink (1.9) --distance square
    WHY:  plink2 alpha dropped --distance; the 1.9 implementation is the
          reference. Square output is the matrix shape ape::nj() expects.
    TUNABLES: (none)
    OUTPUT: {outputs}/structure/Pk.dist + Pk.dist.id
    TRY:    convert to a phylo and inspect a quick `ape::nj` plot to sanity-
            check the cohort splits before Stage-3b adds colour + labels.
    """
    input:
        bed = f"{PATHS['outputs']}/structure/cleaned.bed",
        bim = f"{PATHS['outputs']}/structure/cleaned.bim",
        fam = f"{PATHS['outputs']}/structure/cleaned.fam",
    output:
        dist    = f"{PATHS['outputs']}/structure/Pk.dist",
        dist_id = f"{PATHS['outputs']}/structure/Pk.dist.id",
    log:
        f"{PATHS['logs']}/structure/distance_matrix.log",
    threads: config["compute"]["threads_heavy"]
    params:
        in_prefix  = f"{PATHS['outputs']}/structure/cleaned",
        out_prefix = f"{PATHS['outputs']}/structure/Pk",
    message:
        "[structure] PLINK --distance square"
    shell:
        r"""
        plink --bfile {params.in_prefix} \
            --allow-extra-chr \
            --distance square \
            --threads {threads} \
            --out {params.out_prefix} > {log} 2>&1
        """


# --------------------------------------------------------------------------
# Optional GIS join — wired only when metadata.gis is set
# --------------------------------------------------------------------------

if GIS_FILE:

    rule gis_join:
        """
        Join cluster assignments to a per-sample GIS table (sample_id, lat,
        long). Optional — skipped entirely when metadata.gis is unset.

        WHAT: scripts/R/gis_join.R (port of V1's gis_join.R, role-driven)
        WHY:  Feeds Stage 3b's province / point map. Optional because not
              every cohort has GIS coords; absence should not fail.
        TUNABLES: metadata.gis (path to a TSV with sample_id + lat + long)
        OUTPUT: {outputs}/structure/admix_clusters_gis.tsv
        TRY:    null out metadata.gis and re-run — this rule disappears
                from the DAG and the Stage-3b map rule will skip too.
        """
        input:
            clusters = rules.assign_clusters.output.tsv,
            gis_ref  = GIS_FILE,
            metadata = rules.validate_metadata.output.tsv,
        output:
            tsv = f"{PATHS['outputs']}/structure/admix_clusters_gis.tsv",
        log:
            f"{PATHS['logs']}/structure/gis_join.log",
        params:
            script = str(_AGNOSTIC / "scripts" / "R" / "gis_join.R"),
        message:
            "[structure] Joining clusters to GIS coordinates"
        shell:
            r"""
            Rscript {params.script} \
                {input.clusters} {input.gis_ref} {input.metadata} \
                {output.tsv} > {log} 2>&1
            """
