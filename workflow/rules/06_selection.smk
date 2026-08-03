# 06_selection.smk — Stage 6 selection scan (rehh iHS)
# --------------------------------------------------------------------------
# Model-driven iHS scan. Each entry in config.selection.models declares a
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


rule run_rehh_ihs:
    """
    Per-model whole-genome iHS scan via rehh. Single-threaded — rehh's
    `data2haplohh` / `scan_hh` do not parallelise usefully on this scale.

    WHAT: scripts/R/run_rehh_ihs.R — subsets the combined Stage-4 VCF to
          case + control per the model's role-expression filters, splits
          by contig, accumulates a whole-genome scan, runs
          `calc_candidate_regions`. Emits empty (header-only) TSVs +
          a one-line summary on a negative result.
    WHY:  iHS detects loci under recent positive selection. Splitting by
          case/control lets the same machinery run any comparison the
          user cares about — no cohort hardcoding beyond the filter
          expression the user writes in config.
    TUNABLES: selection.ihs_threshold, selection.ihs_p_threshold,
              selection.window_size_bp, selection.window_overlap_bp,
              selection.min_extr_markers, selection.models[*].{case,control}_filter
    OUTPUT: {outputs}/selection/{model}/candidate_regions_iHS.tsv,
            {outputs}/selection/{model}/ihs_table.tsv,
            {outputs}/selection/{model}/ihs_summary.tsv
    TRY:    add a second model to selection.models (e.g. a
            North-Kalimantan-vs-Sabah scan) and re-run — the rule
            expands over the added name automatically.
    """
    input:
        vcf      = rules.combined_vcf.output.vcf,
        metadata = rules.validate_metadata.output.tsv,
        clusters = f"{PATHS['outputs']}/structure/admix_clusters.tsv",
    output:
        candidates = f"{PATHS['outputs']}/selection/{{model}}/candidate_regions_iHS.tsv",
        ihs_table  = f"{PATHS['outputs']}/selection/{{model}}/ihs_table.tsv",
        summary    = f"{PATHS['outputs']}/selection/{{model}}/ihs_summary.tsv",
    log:
        f"{PATHS['logs']}/selection/run_rehh_ihs_{{model}}.log",
    params:
        out_dir      = f"{PATHS['outputs']}/selection",
        case_filter  = lambda wc: _model_cfg(wc.model)["case_filter"],
        ctrl_filter  = lambda wc: _model_cfg(wc.model)["control_filter"],
        threshold    = SELECTION.get("ihs_threshold", 4),
        p_threshold  = SELECTION.get("ihs_p_threshold", 0.0001),
        window_size  = SELECTION.get("window_size_bp", 10000),
        overlap      = SELECTION.get("window_overlap_bp", 1000),
        min_extr_mrk = SELECTION.get("min_extr_markers", 3),
        script       = str(_AGNOSTIC / "scripts" / "R" / "run_rehh_ihs.R"),
    threads: 1
    message:
        "[selection] iHS scan: {wildcards.model}"
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
            --threshold      {params.threshold} \
            --p-threshold    {params.p_threshold} \
            --window-size    {params.window_size} \
            --overlap        {params.overlap} \
            --min-extr-mrk   {params.min_extr_mrk} \
            --out-dir        {params.out_dir} \
            > {log} 2>&1
        """


rule plot_ihs_scan:
    """
    Genome-wide iHS scan + -log10(p) plot for a model.

    WHAT: scripts/R/plot_ihs_scan.R
    WHY:  The scan is the standard visualisation; alternating chromosome
          fills + 0.1% / 99.9% quantile dotted lines flag outliers.
    TUNABLES: (none in the rule; cosmetic in the script)
    OUTPUT: {outputs}/selection/{model}/ihs_scan.{png,svg} +
            {outputs}/selection/{model}/ihs_pvalue.{png,svg}
    TRY:    inspect the p-value plot for peaks above the 0.001 quantile
            line; the candidate_regions_iHS.tsv should hit those peaks.
    """
    input:
        ihs_table = rules.run_rehh_ihs.output.ihs_table,
    output:
        scan_png = f"{PATHS['outputs']}/selection/{{model}}/ihs_scan.png",
        pval_png = f"{PATHS['outputs']}/selection/{{model}}/ihs_pvalue.png",
        scan_svg = f"{PATHS['outputs']}/selection/{{model}}/ihs_scan.svg",
        pval_svg = f"{PATHS['outputs']}/selection/{{model}}/ihs_pvalue.svg",
    log:
        f"{PATHS['logs']}/selection/plot_ihs_scan_{{model}}.log",
    params:
        script = str(_AGNOSTIC / "scripts" / "R" / "plot_ihs_scan.R"),
    message:
        "[selection] iHS scan + p-value plots: {wildcards.model}"
    shell:
        r"""
        mkdir -p $(dirname {output.scan_png}) $(dirname {log})
        Rscript {params.script} \
            {input.ihs_table} \
            {output.scan_png} {output.pval_png} \
            {output.scan_svg} {output.pval_svg} \
            > {log} 2>&1
        """
