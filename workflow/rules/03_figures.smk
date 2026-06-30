# 03_figures.smk — Stage 3 figures (role-driven, role-presence-guarded)
# --------------------------------------------------------------------------
# Optional figure rules consuming the Stage 3 core outputs
# (admix_clusters.tsv, Pk.eigenvec, Pk.dist, admix_clusters_gis.tsv). Each
# rule is wired in only if its backing role(s) are configured — same
# parse-time guard pattern used in Stage 2.
#
# Generic palette generator (`scripts/R/palettes.R`) gives a deterministic
# level → colour map per role; no Malaysian/Indonesian palette hardcodes.
# --------------------------------------------------------------------------

HAS_GROUP_FIG     = bool(ROLES.get("group"))         # implicit; .Q always has Cluster
HAS_COUNTRY_FIG   = bool(ROLES.get("country"))
HAS_GEOGRAPHY_FIG = bool(ROLES.get("geography"))
HAS_GIS_FIG       = bool(META.get("gis")) and HAS_GEOGRAPHY_FIG and HAS_COUNTRY_FIG


rule plot_admixture_bars:
    """
    Single-panel stacked ADMIXTURE bar plot — samples ordered by dominant
    cluster then dominant proportion, coloured by component.

    WHAT: scripts/R/plot_admixture_bars.R
    WHY:  The slide-7 cohort-wide ancestry summary. Always renders; no
          role gating because dominant-cluster sort uses admix_clusters.tsv.
    TUNABLES: (cosmetic-only)
    OUTPUT: {reports}/figures/admixture_bars.png/.svg
    TRY:    if the bar plot looks "noisy" you may want a different K — pass
            it via `structure.admixture_k` and re-run.
    """
    input:
        bk       = rules.select_best_k.output.txt,
        fam      = f"{PATHS['outputs']}/structure/cleaned.fam",
        clusters = rules.assign_clusters.output.tsv,
    output:
        png = f"{PATHS['reports']}/figures/admixture_bars.png",
        svg = f"{PATHS['reports']}/figures/admixture_bars.svg",
    log:
        f"{PATHS['logs']}/structure/plot_admixture_bars.log",
    params:
        admix_dir = f"{PATHS['outputs']}/structure/admixture",
        script    = str(_AGNOSTIC / "scripts" / "R" / "plot_admixture_bars.R"),
    message:
        "[structure:figures] ADMIXTURE bars (single panel)"
    shell:
        r"""
        K=$(cat {input.bk})
        Q="{params.admix_dir}/cleaned.${{K}}.Q"
        mkdir -p $(dirname {output.png})
        Rscript {params.script} \
            "$Q" \
            {input.fam} \
            {input.clusters} \
            {output.png} {output.svg} \
            > {log} 2>&1
        """


if HAS_COUNTRY_FIG:

    rule plot_admixture_bars_by_country:
        """
        ADMIXTURE bars faceted by the country role.

        WHAT: scripts/R/plot_admixture_bars_by_role.R … country
        WHY:  Shows how ancestry distributions differ across cohorts'
              countries of origin. Drops out gracefully if the country
              role is null in config.
        TUNABLES: metadata.roles.country
        OUTPUT: {reports}/figures/admixture_bars_by_country.png/.svg
        TRY:    null metadata.roles.country → rule disappears from DAG.
        """
        input:
            bk       = rules.select_best_k.output.txt,
            fam      = f"{PATHS['outputs']}/structure/cleaned.fam",
            clusters = rules.assign_clusters.output.tsv,
            metadata = rules.validate_metadata.output.tsv,
        output:
            png = f"{PATHS['reports']}/figures/admixture_bars_by_country.png",
            svg = f"{PATHS['reports']}/figures/admixture_bars_by_country.svg",
        log:
            f"{PATHS['logs']}/structure/plot_admixture_bars_by_country.log",
        params:
            admix_dir = f"{PATHS['outputs']}/structure/admixture",
            script    = str(_AGNOSTIC / "scripts" / "R" / "plot_admixture_bars_by_role.R"),
        message:
            "[structure:figures] ADMIXTURE bars faceted by country"
        shell:
            r"""
            K=$(cat {input.bk})
            Q="{params.admix_dir}/cleaned.${{K}}.Q"
            mkdir -p $(dirname {output.png})
            Rscript {params.script} \
                "$Q" \
                {input.fam} \
                {input.clusters} \
                {input.metadata} \
                country \
                {output.png} {output.svg} \
                > {log} 2>&1
            """


if HAS_GEOGRAPHY_FIG:

    rule plot_admixture_bars_by_geography:
        """
        ADMIXTURE bars faceted by the geography role.

        WHAT: scripts/R/plot_admixture_bars_by_role.R … geography
        WHY:  V1's "by province" plot, generalised — facets over whatever
              geography levels exist, not a fixed 6-province assumption.
        TUNABLES: metadata.roles.geography
        OUTPUT: {reports}/figures/admixture_bars_by_geography.png/.svg
        TRY:    null metadata.roles.geography → rule disappears from DAG.
        """
        input:
            bk       = rules.select_best_k.output.txt,
            fam      = f"{PATHS['outputs']}/structure/cleaned.fam",
            clusters = rules.assign_clusters.output.tsv,
            metadata = rules.validate_metadata.output.tsv,
        output:
            png = f"{PATHS['reports']}/figures/admixture_bars_by_geography.png",
            svg = f"{PATHS['reports']}/figures/admixture_bars_by_geography.svg",
        log:
            f"{PATHS['logs']}/structure/plot_admixture_bars_by_geography.log",
        params:
            admix_dir = f"{PATHS['outputs']}/structure/admixture",
            script    = str(_AGNOSTIC / "scripts" / "R" / "plot_admixture_bars_by_role.R"),
        message:
            "[structure:figures] ADMIXTURE bars faceted by geography"
        shell:
            r"""
            K=$(cat {input.bk})
            Q="{params.admix_dir}/cleaned.${{K}}.Q"
            mkdir -p $(dirname {output.png})
            Rscript {params.script} \
                "$Q" \
                {input.fam} \
                {input.clusters} \
                {input.metadata} \
                geography \
                {output.png} {output.svg} \
                > {log} 2>&1
            """


rule plot_pca:
    """
    PC1-vs-PC2 scatter coloured by ADMIXTURE cluster.

    WHAT: scripts/R/plot_pca.R
    WHY:  Always renders — uses admix_clusters.tsv's Cluster column which
          is always produced (in numbered mode at worst).
    TUNABLES: (none — cosmetic only)
    OUTPUT: {reports}/figures/pca.png/.svg
    TRY:    swap the colour role to geography in plot_pca_by_geography
            (below) to see if PC1/PC2 maps onto region instead of cluster.
    """
    input:
        eigenvec = rules.pca.output.eigenvec,
        varpct   = rules.pca_variance.output.tsv,
        clusters = rules.assign_clusters.output.tsv,
    output:
        png = f"{PATHS['reports']}/figures/pca.png",
        svg = f"{PATHS['reports']}/figures/pca.svg",
    log:
        f"{PATHS['logs']}/structure/plot_pca.log",
    params:
        script = str(_AGNOSTIC / "scripts" / "R" / "plot_pca.R"),
    message:
        "[structure:figures] PCA scatter (by cluster)"
    shell:
        r"""
        mkdir -p $(dirname {output.png})
        Rscript {params.script} \
            {input.eigenvec} {input.varpct} {input.clusters} \
            {output.png} {output.svg} \
            > {log} 2>&1
        """


if HAS_GEOGRAPHY_FIG:

    rule plot_pca_by_geography:
        """
        PCA scatter coloured by the geography role.

        WHAT: scripts/R/plot_pca_by_geography.R
        WHY:  Lets the user see whether PC1/PC2 separates samples by
              geographic origin as well as by ADMIXTURE cluster.
        TUNABLES: metadata.roles.geography
        OUTPUT: {reports}/figures/pca_by_geography.png/.svg
        TRY:    null metadata.roles.geography → rule disappears from DAG.
        """
        input:
            eigenvec = rules.pca.output.eigenvec,
            varpct   = rules.pca_variance.output.tsv,
            metadata = rules.validate_metadata.output.tsv,
        output:
            png = f"{PATHS['reports']}/figures/pca_by_geography.png",
            svg = f"{PATHS['reports']}/figures/pca_by_geography.svg",
        log:
            f"{PATHS['logs']}/structure/plot_pca_by_geography.log",
        params:
            script = str(_AGNOSTIC / "scripts" / "R" / "plot_pca_by_geography.R"),
        message:
            "[structure:figures] PCA scatter (by geography)"
        shell:
            r"""
            mkdir -p $(dirname {output.png})
            Rscript {params.script} \
                {input.eigenvec} {input.varpct} {input.metadata} \
                {output.png} {output.svg} \
                > {log} 2>&1
            """


rule plot_nj_tree:
    """
    Unrooted Neighbour-Joining tree from PLINK's IBS distance matrix,
    coloured by ADMIXTURE cluster.

    WHAT: scripts/R/plot_nj_tree.R (ape::nj + ggtree daylight layout)
    WHY:  Visualises clade structure without assuming an outgroup.
          Coloured edges show within-cluster vs between-cluster edges.
    TUNABLES: (none — colour palette auto-generated)
    OUTPUT: {reports}/figures/njt.png/.svg
    TRY:    open the svg and inspect: tight cluster cliques + sparse
            connections between clusters typically mean clean structure.
    """
    input:
        dist    = rules.distance_matrix.output.dist,
        dist_id = rules.distance_matrix.output.dist_id,
        clusters= rules.assign_clusters.output.tsv,
    output:
        png = f"{PATHS['reports']}/figures/njt.png",
        svg = f"{PATHS['reports']}/figures/njt.svg",
    log:
        f"{PATHS['logs']}/structure/plot_nj_tree.log",
    params:
        script = str(_AGNOSTIC / "scripts" / "R" / "plot_nj_tree.R"),
    message:
        "[structure:figures] Neighbour-Joining tree"
    shell:
        r"""
        mkdir -p $(dirname {output.png})
        Rscript {params.script} \
            {input.dist} {input.dist_id} {input.clusters} \
            {output.png} {output.svg} \
            > {log} 2>&1
        """


if HAS_GIS_FIG:

    rule plot_province_map:
        """
        Admin-1 (province / state) choropleth + jittered sample points.

        WHAT: scripts/R/plot_province_map.R — rnaturalearth ne_states()
              over the countries seen in the data; map extent from sample
              bbox; region fill by geography role; point colour by cluster.
        WHY:  Geographic context for the structure result. Generalised
              from V1's Malaysia+Indonesia hardcode: any cohort with
              metadata.gis + country + geography roles renders.
        TUNABLES: metadata.gis, metadata.roles.{country,geography}
        OUTPUT: {reports}/figures/admix_map.png/.svg
        TRY:    drop metadata.gis from the config — this rule and gis_join
                both disappear from the DAG. Add lat/long to the GIS file
                for new samples to extend the map without touching code.
        """
        input:
            gis = rules.gis_join.output.tsv,
        output:
            png = f"{PATHS['reports']}/figures/admix_map.png",
            svg = f"{PATHS['reports']}/figures/admix_map.svg",
        log:
            f"{PATHS['logs']}/structure/plot_province_map.log",
        params:
            script = str(_AGNOSTIC / "scripts" / "R" / "plot_province_map.R"),
        message:
            "[structure:figures] Province choropleth + sample points"
        shell:
            r"""
            mkdir -p $(dirname {output.png})
            Rscript {params.script} \
                {input.gis} {output.png} {output.svg} \
                > {log} 2>&1
            """
