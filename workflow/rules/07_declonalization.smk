# 07_declonalization.smk — full vs unique comparison outputs
# --------------------------------------------------------------------------
# The point of always producing BOTH sample sets is that the difference is
# visible. These rules read the same key numbers out of each stage's `full`
# and `unique` namespaces and put them side by side with the delta, so a
# reader can judge how much of any given claim was clonal pseudo-replication
# rather than taking the pipeline on faith.
#
# They compute nothing scientific — every statistic keeps its own null, its
# own FDR and its own script. These rules count rows and read numbers the
# stages already wrote. Centralising the comparison means the comparison
# logic cannot drift between stages.
#
# Wired only when clonality.declonalize is true (see the Snakefile): with one
# arm there is nothing to compare.
# --------------------------------------------------------------------------

DECLON_STAGE_DIRS = {
    "structure":     f"{PATHS['outputs']}/structure",
    "introgression": f"{PATHS['outputs']}/introgression",
    "moi":           f"{PATHS['outputs']}/moi",
    "selection":     f"{PATHS['outputs']}/selection",
}


def _declon_inputs(wildcards):
    """
    The stage outputs the comparison reads. Declared as real inputs so the
    comparison can never run against a half-built tree, and so `snakemake
    <comparison.tsv>` pulls the whole stage — both arms — behind it.
    """
    stage = wildcards.stage
    out = []
    for ss in SAMPLE_SETS:
        if stage == "structure":
            out += [f"{PATHS['outputs']}/structure/{ss}/admix_clusters.tsv",
                    f"{PATHS['outputs']}/structure/{ss}/best_k.txt",
                    f"{PATHS['outputs']}/structure/{ss}/pca_variance.tsv"]
        elif stage == "introgression":
            out += [f"{PATHS['outputs']}/introgression/{ss}/introgressed_windows_filtered.tsv",
                    f"{PATHS['outputs']}/introgression/{ss}/windows_by_cluster.tsv"]
            if INTRO_FOCAL:
                out += [f"{PATHS['outputs']}/introgression/{ss}/focal_{INTRO_FOCAL}_enriched_windows.tsv"]
        elif stage == "moi":
            if ROLES.get("country"):
                out += [f"{PATHS['outputs']}/moi/{ss}/fws_by_country.tsv"]
            if ROLES.get("geography"):
                out += [f"{PATHS['outputs']}/moi/{ss}/fws_by_geography.tsv"]
        elif stage == "selection":
            for m in SELECTION_MODELS:
                out += [f"{PATHS['outputs']}/selection/{ss}/{m}/candidate_regions_iHS.tsv",
                        f"{PATHS['outputs']}/selection/{ss}/{m}/ihs_table.tsv"]
    return out


rule declonalization_comparison:
    """
    Full-vs-unique key numbers for one stage, side by side.

    WHAT: scripts/R/declonalization_comparison.R — per stage, one row per
          metric with the full-set value, the de-clonalized value, the delta
          and the percent change. Structure compares best K, cluster count,
          per-cluster sizes and PC1/PC2 variance; introgression compares
          window counts overall, per cluster, and the focal-enrichment result;
          moi compares per-group polyclonality rates; selection compares
          candidate regions and scanned markers per model.
    WHY:  this is the artifact that makes "both outputs" worth producing. A
          metric that barely moves was not driven by clonality; one that moves
          a lot was substantially pseudo-replication. Neither is automatically
          wrong — but the difference should never be invisible.
    READ: a metric that RISES on the de-clonalized set is worth attention. It
          usually means a threshold that scales with cluster size (Stage 5's
          per-cluster support floor is a fraction/floor over cluster n) got
          easier to clear once the clonal padding was gone.
    TUNABLES: clonality.declonalize
    OUTPUT: {outputs}/<stage>/declonalization_comparison.tsv
    TRY:    read n_samples_removed first — it is the size of the intervention.
            If it is small and a metric still moved a lot, that metric was
            resting on very few independent genotypes.
    """
    input:
        stage_outputs = _declon_inputs,
        audit         = f"{CLONALITY_DIR}/declonalization_audit.tsv",
        full_keep     = f"{CLONALITY_DIR}/all_genotypes.txt",
        unique_keep   = f"{CLONALITY_DIR}/unique_genotypes.txt",
    output:
        tsv = f"{PATHS['outputs']}/{{stage}}/declonalization_comparison.tsv",
    log:
        f"{PATHS['logs']}/clonality/comparison_{{stage}}.log",
    params:
        full_dir   = lambda wc: f"{DECLON_STAGE_DIRS[wc.stage]}/full",
        unique_dir = lambda wc: f"{DECLON_STAGE_DIRS[wc.stage]}/unique",
        script     = str(_AGNOSTIC / "scripts" / "R" / "declonalization_comparison.R"),
    wildcard_constraints:
        stage = r"structure|introgression|moi|selection",
    message:
        "[clonality] full-vs-unique comparison: {wildcards.stage}"
    shell:
        r"""
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        Rscript {params.script} \
            --stage       {wildcards.stage} \
            --full-dir    {params.full_dir} \
            --unique-dir  {params.unique_dir} \
            --audit       {input.audit} \
            --full-keep   {input.full_keep} \
            --unique-keep {input.unique_keep} \
            --out         {output.tsv} \
            > {log} 2>&1
        """


rule plot_declonalization_comparison:
    """
    One dumbbell panel per stage: how far each metric moved.

    WHAT: scripts/R/plot_declonalization_comparison.R over every stage's
          comparison TSV. Each row is a metric, the two dots are its full-set
          and de-clonalized values, the segment between them is the effect.
    WHY:  nobody scans four TSVs before reading a result. Coincident dots =
          clonality did not drive it; a long segment = it did.
    OUTPUT: {reports}/figures/declonalization_comparison.{png,svg}
    TRY:    if a per-cluster row moves far more than the cohort rows, the
            clonality is concentrated in that cluster — worth checking against
            declonalization_audit.tsv to see which groups sat there.
    """
    input:
        tsvs = ([f"{PATHS['outputs']}/{st}/declonalization_comparison.tsv"
                 for st in ("structure", "introgression", "moi")]
                + ([f"{PATHS['outputs']}/selection/declonalization_comparison.tsv"]
                   if SELECTION_MODELS else [])),
    output:
        png = f"{PATHS['reports']}/figures/declonalization_comparison.png",
        svg = f"{PATHS['reports']}/figures/declonalization_comparison.svg",
    log:
        f"{PATHS['logs']}/clonality/plot_comparison.log",
    params:
        script = str(_AGNOSTIC / "scripts" / "R" / "plot_declonalization_comparison.R"),
    message:
        "[clonality] Plotting the full-vs-unique comparison"
    shell:
        r"""
        mkdir -p $(dirname {output.png}) $(dirname {log})
        Rscript {params.script} {output.png} {output.svg} {input.tsvs} > {log} 2>&1
        """
