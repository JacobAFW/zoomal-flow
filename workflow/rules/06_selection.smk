# 06_selection.smk — Stage 6 selection scan (rehh XP-EHH / Rsb / iHS)
# --------------------------------------------------------------------------
# Model-driven cross-population scan. Each entry in config.selection.models
# declares a
# `name`, `case_filter`, and `control_filter`. The filter strings are R
# expressions evaluated against the canonical role columns in the
# validated metadata table joined with admix_clusters.tsv, e.g.
#
#   case_filter:    "group == 'Peninsular' & geography == 'Aceh'"
#   control_filter: "group == 'Peninsular' & geography != 'Aceh'"
#
# Rules expand over the model names. Config default ships an empty list;
# the Indo cohort example ships the Aceh_vs_other_Peninsular model.
# --------------------------------------------------------------------------

SELECTION       = config.get("selection", {})
SELECTION_MODELS_LIST = SELECTION.get("models", [])
SELECTION_MODELS      = [m["name"] for m in SELECTION_MODELS_LIST]

def _model_cfg(name):
    for m in SELECTION_MODELS_LIST:
        if m["name"] == name:
            return m
    raise ValueError(f"Unknown selection model: {name!r}")


rule run_selection_scan:
    """
    Per-model cross-population selection scan. REPLACED the iHS-only rule on
    2026-10-21; see CONTEXT-pipeline.md for the four defects that forced it.

    WHAT: scripts/R/run_selection_scan.R — splits case and control per the
          model's role-expression filters, haploidises, scans each group
          SEPARATELY with rehh, then computes XP-EHH and Rsb for the
          two-group contrast and a within-population iHS per group.
    WHY:  iHS is a WITHIN-population statistic and cannot test a two-group
          contrast; the old rule pooled the two groups and scanned the pool,
          so the case/control split never entered the result. XP-EHH
          (ies2xpehh) and Rsb (ines2rsb) are the cross-population statistics
          that do.
    INPUT: the POST-QC SNP set, not the Stage-3 `cleaned` bfile. The cleaned
          bfile is filtered cohort-wide for ADMIXTURE and PCA and is roughly
          one marker per 6.5 kb, which is too sparse for EHH integration to
          complete — the old rule's scan was empty, not negative.
    TUNABLES: selection.{min_callrate,min_maf,max_het_rate,ihs_min_maf,
              threshold,window_size_bp,window_overlap_bp,min_extr_markers,
              min_group_n} and selection.models[*].{case,control}_filter
    OUTPUT: {outputs}/selection/{sampleset}/{model}/
              scan_summary.tsv        headline table, one row per statistic
              tail_diagnostics.tsv    kurtosis / max |z| / expected max /
                                      tail excess per statistic
              candidate_regions_{XPEHH,RSB,iHS_case,iHS_control}.tsv
              xpehh_table.tsv, rsb_table.tsv, ihs_{case,control}_table.tsv
              marker_manifest.tsv     markers in and out per contig
    READ THIS FIRST: `scan_summary.tsv` carries a `verdict` column per
          statistic. rehh standardizes to mean 0 / SD 1 and reads p-values off
          a normal; when two groups are similar that denominator collapses and
          trivial differences become extreme z-scores — a weakly-differentiated
          contrast can return hundreds of "candidate regions" that are
          artefacts of a vanishing denominator. A candidate count whose
          verdict is DEGENERATE means nothing.
    TRY:    add a second model to selection.models and re-run — the rule
            expands over the added name automatically.
    """
    input:
        vcf      = f"{PATHS['outputs']}/qc/snps.qc.vcf.gz",
        idx      = f"{PATHS['outputs']}/qc/snps.qc.vcf.gz.csi",
        metadata = rules.validate_metadata.output.tsv,
        # Cluster definitions are the FULL-set backbone in both arms — the
        # de-clonalized branch changes WHO is in the case/control sets, not
        # what the clusters are (docs/clonality.md).
        clusters = f"{PATHS['outputs']}/structure/full/admix_clusters.tsv",
        exclude  = lambda wc: sampleset_exclude(wc.sampleset),
    output:
        summary    = f"{PATHS['outputs']}/selection/{{sampleset}}/{{model}}/scan_summary.tsv",
        diagnostics = f"{PATHS['outputs']}/selection/{{sampleset}}/{{model}}/tail_diagnostics.tsv",
        candidates = f"{PATHS['outputs']}/selection/{{sampleset}}/{{model}}/candidate_regions_XPEHH.tsv",
        manifest   = f"{PATHS['outputs']}/selection/{{sampleset}}/{{model}}/marker_manifest.tsv",
        # Declared so the plot rule can depend on it by name; without a
        # declared producer Snakemake cannot resolve it and the DAG fails.
        xpehh      = f"{PATHS['outputs']}/selection/{{sampleset}}/{{model}}/xpehh_table.tsv",
    log:
        f"{PATHS['logs']}/selection/{{sampleset}}/run_selection_scan_{{model}}.log",
    params:
        out_dir      = f"{PATHS['outputs']}/selection/{{sampleset}}",
        case_filter  = lambda wc: _model_cfg(wc.model)["case_filter"],
        ctrl_filter  = lambda wc: _model_cfg(wc.model)["control_filter"],
        min_callrate = SELECTION.get("min_callrate", 0.9),
        min_maf      = SELECTION.get("min_maf", 0),
        max_het_rate = SELECTION.get("max_het_rate", 0.1),
        ihs_min_maf  = SELECTION.get("ihs_min_maf", 0.05),
        threshold    = SELECTION.get("ihs_threshold", 4),
        window_size  = SELECTION.get("window_size_bp", 10000),
        overlap      = SELECTION.get("window_overlap_bp", 1000),
        min_extr_mrk = SELECTION.get("min_extr_markers", 3),
        min_group_n  = SELECTION.get("min_group_n", 5),
        script       = str(_AGNOSTIC / "scripts" / "R" / "run_selection_scan.R"),
    threads: 1
    message:
        "[selection:{wildcards.sampleset}] XP-EHH + Rsb scan: {wildcards.model}"
    shell:
        r"""
        mkdir -p {params.out_dir}/{wildcards.model} $(dirname {log})
        Rscript {params.script} \
            --vcf            {input.vcf} \
            --metadata       {input.metadata} \
            --clusters       {input.clusters} \
            --model          "{wildcards.model}" \
            --case-filter    "{params.case_filter}" \
            --control-filter "{params.ctrl_filter}" \
            --exclude        {input.exclude} \
            --min-callrate   {params.min_callrate} \
            --min-maf        {params.min_maf} \
            --max-het-rate   {params.max_het_rate} \
            --ihs-min-maf    {params.ihs_min_maf} \
            --threshold      {params.threshold} \
            --window-size    {params.window_size} \
            --overlap        {params.overlap} \
            --min-extr-mrk   {params.min_extr_mrk} \
            --min-group-n    {params.min_group_n} \
            --out-dir        {params.out_dir} \
            > {log} 2>&1
        """


rule plot_selection_scan:
    """
    Genome-wide iHS scan + -log10(p) plot for a model.

    WHAT: scripts/R/plot_ihs_scan.R (generic: plots whatever standardized
          statistic column the table carries)
    WHY:  The scan is the standard visualisation; alternating chromosome
          fills + 0.1% / 99.9% quantile dotted lines flag outliers.
    TUNABLES: (none in the rule; cosmetic in the script)
    OUTPUT: {outputs}/selection/{model}/ihs_scan.{png,svg} +
            {outputs}/selection/{model}/ihs_pvalue.{png,svg}
    TRY:    inspect the p-value plot for peaks above the 0.001 quantile
            line; the candidate_regions_iHS.tsv should hit those peaks.
    """
    input:
        table = rules.run_selection_scan.output.summary,
        xpehh = rules.run_selection_scan.output.xpehh,
    output:
        scan_png = f"{PATHS['reports']}/figures/{{sampleset}}/selection_{{model}}_scan.png",
        pval_png = f"{PATHS['reports']}/figures/{{sampleset}}/selection_{{model}}_pvalue.png",
        scan_svg = f"{PATHS['reports']}/figures/{{sampleset}}/selection_{{model}}_scan.svg",
        pval_svg = f"{PATHS['reports']}/figures/{{sampleset}}/selection_{{model}}_pvalue.svg",
    log:
        f"{PATHS['logs']}/selection/{{sampleset}}/plot_ihs_scan_{{model}}.log",
    params:
        script = str(_AGNOSTIC / "scripts" / "R" / "plot_ihs_scan.R"),
    message:
        "[selection:{wildcards.sampleset}] scan + p-value plots: {wildcards.model}"
    shell:
        r"""
        mkdir -p $(dirname {output.scan_png}) $(dirname {log})
        Rscript {params.script} \
            {input.xpehh} \
            {output.scan_png} {output.pval_png} \
            {output.scan_svg} {output.pval_svg} \
            > {log} 2>&1
        """
