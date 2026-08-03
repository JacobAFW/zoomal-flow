# 04_figures.smk — Stage 4b figures (role-driven, role-presence-guarded)
# --------------------------------------------------------------------------
# Optional figure rules consuming Stage 4's per-cluster hmm_fract.txt files
# and clonal_clusters.tsv. Each rule is wired in only if its backing role
# is configured — same parse-time guard pattern used in Stages 2 + 3b.
# --------------------------------------------------------------------------

HAS_GEOGRAPHY_IBD = bool(ROLES.get("geography"))
HAS_DATE_IBD      = bool(ROLES.get("date"))
FOCAL_CLUSTER     = IBD.get("focal_cluster")


def _connectivity_clusters(wildcards):
    """
    Clusters we render connectivity plots for. If `ibd.focal_cluster` is set
    just that one; otherwise every cluster clearing `ibd.min_cluster_n` from
    the checkpoint.
    """
    if FOCAL_CLUSTER:
        return [FOCAL_CLUSTER]
    admix_tsv = checkpoints.assign_clusters.get(**wildcards).output.tsv
    return _ibd_clusters(admix_tsv)


def ibd_connectivity_targets(wildcards):
    """FINAL_TARGETS input function: expand PDF+PNG per rendered cluster."""
    if not HAS_GEOGRAPHY_IBD:
        return []
    clusters = _connectivity_clusters(wildcards)
    out = []
    for c in clusters:
        out.append(f"{PATHS['reports']}/figures/ibd_connectivity_{c}.pdf")
        out.append(f"{PATHS['reports']}/figures/ibd_connectivity_{c}.png")
    return out


if HAS_GEOGRAPHY_IBD:

    rule plot_ibd_connectivity:
        """
        Two-panel igraph connectivity network for a cluster's IBD graph
        (IBD >= 5% and >= 95%), nodes coloured by the geography role.

        WHAT: scripts/R/plot_ibd_connectivity.R
        WHY:  V1's slide-10 connectivity is the go-to picture of clonal
              vs non-clonal structure. Node colour by geography via the
              palette generator so the colour key travels with the PCA /
              map figures automatically.
        TUNABLES: ibd.focal_cluster, metadata.roles.geography
        OUTPUT: {reports}/figures/ibd_connectivity_{{cluster}}.pdf/.png
        TRY:    null metadata.roles.geography → rule disappears from DAG.
                Set ibd.focal_cluster to a different cluster to re-render
                its network.
        """
        input:
            fract    = f"{PATHS['outputs']}/ibd/{{cluster}}/{{cluster}}.hmm_fract.txt",
            metadata = rules.validate_metadata.output.tsv,
        output:
            pdf = f"{PATHS['reports']}/figures/ibd_connectivity_{{cluster}}.pdf",
            png = f"{PATHS['reports']}/figures/ibd_connectivity_{{cluster}}.png",
        log:
            f"{PATHS['logs']}/ibd/plot_ibd_connectivity_{{cluster}}.log",
        params:
            script = str(_AGNOSTIC / "scripts" / "R" / "plot_ibd_connectivity.R"),
        message:
            "[ibd:figures] Connectivity network: {wildcards.cluster}"
        shell:
            r"""
            mkdir -p $(dirname {output.pdf})
            Rscript {params.script} \
                {input.fract} {input.metadata} {wildcards.cluster} \
                {output.pdf} {output.png} \
                > {log} 2>&1
            """


if HAS_DATE_IBD:

    rule plot_clonal_temporal:
        """
        Stacked bar of clonal-group sampling dates, driven entirely by the
        computed clonal_clusters.tsv joined to the canonical date role.

        WHAT: scripts/R/plot_clonal_temporal.R
        WHY:  V1 hardcoded an 8-sample tribble + a Sabang xlsx date lookup;
              neither ports. The agnostic version reads clonal groups
              straight from Stage 4's connected-components output and
              joins to the date role on the fly. Undated samples in a
              clonal group are excluded from the plot with a logged note
              (never silently dropped).
        TUNABLES: ibd.clonal_ibd_threshold (upstream), metadata.roles.date
        OUTPUT: {reports}/figures/clonal_temporal.png/.svg
        TRY:    null metadata.roles.date → rule disappears from DAG.
                Lower ibd.clonal_ibd_threshold to widen the definition of
                a clonal group and see the temporal pattern change.
        """
        input:
            clusters_tsv = rules.clonal_clusters.output.clusters_tsv,
        output:
            png = f"{PATHS['reports']}/figures/clonal_temporal.png",
            svg = f"{PATHS['reports']}/figures/clonal_temporal.svg",
        log:
            f"{PATHS['logs']}/ibd/plot_clonal_temporal.log",
        params:
            script = str(_AGNOSTIC / "scripts" / "R" / "plot_clonal_temporal.R"),
        message:
            "[ibd:figures] Clonal temporal stacked bar"
        shell:
            r"""
            mkdir -p $(dirname {output.png})
            Rscript {params.script} \
                {input.clusters_tsv} \
                {output.png} {output.svg} \
                > {log} 2>&1
            """
