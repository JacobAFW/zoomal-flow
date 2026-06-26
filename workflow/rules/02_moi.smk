# 02_moi.smk — Stage 2 MOI / Fws, common rules + WGS/microhap include seam
# --------------------------------------------------------------------------
# SEAM CONTRACT (DESIGN §4, seam 2):
#   Both moi_wgs.smk and moi_microhap.smk MUST emit
#     {outputs}/moi/fws_MOI.tsv
#   with exactly two tab-separated columns and a header row:
#       sample      Proportion
#   `sample` matches VCF / metadata `sample_id`; `Proportion` is the
#   within-host-complexity score (Fws for WGS; the microhap-aware analog
#   for the microhap path). Lower = more polyclonal.
#
# Stage 3+ and the common rules below read exclusively from fws_MOI.tsv
# and do not branch on input_type — the fork is contained here.
#
# Common rules in this file (run regardless of input_type):
#   - fws_high_moi_list  — awk over fws_MOI.tsv → exclude_high_moi.txt
#   - fws_summary        — role-driven country/geography summaries
#   - plot_fws_density   — role-driven density panels (drops gracefully
#                          when a role is absent)
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# input_type seam — pick the MOI implementation
# --------------------------------------------------------------------------
_MOI_INPUT_TYPE = COHORT["input_type"]
if _MOI_INPUT_TYPE == "wgs":
    include: "moi_wgs.smk"
elif _MOI_INPUT_TYPE == "microhap":
    include: "moi_microhap.smk"
else:
    # Schema enforces the enum; this is belt-and-braces.
    raise ValueError(
        f"Unknown cohort.input_type={_MOI_INPUT_TYPE!r}; expected 'wgs' or 'microhap'."
    )


# --------------------------------------------------------------------------
# Common rules
# --------------------------------------------------------------------------

# Role presence drives whether the country / geography summaries + plots
# are wired in. ROLES is set in the Snakefile.
HAS_COUNTRY   = bool(ROLES.get("country"))
HAS_GEOGRAPHY = bool(ROLES.get("geography"))


rule fws_high_moi_list:
    """
    Samples with Fws < moi.fws_exclusion_cutoff are polyclonal and violate
    the single-genotype assumption of ADMIXTURE / PCA / hmmIBD in Stage 3+.
    This list is the Stage 2 → Stage 3 handoff.

    WHAT: awk 'NR>1 && $2 < <cutoff>' fws_MOI.tsv → exclude list
    WHY:  Single-genotype methods downstream silently produce wrong answers
          on polyclonal samples. Better to filter them up front than rely on
          each downstream rule to remember the threshold.
    TUNABLES: moi.fws_exclusion_cutoff
    OUTPUT: {outputs}/moi/exclude_high_moi.txt
    TRY:    raise fws_exclusion_cutoff to 0.85 (V1's HPC value) to keep
            mildly-polyclonal samples in the Stage-3 input and see how it
            shifts ADMIXTURE cluster counts.
    """
    input:
        fws = f"{PATHS['outputs']}/moi/fws_MOI.tsv",
    output:
        exclude = f"{PATHS['outputs']}/moi/exclude_high_moi.txt",
    log:
        f"{PATHS['logs']}/moi/fws_high_moi_list.log",
    params:
        cutoff = config["moi"]["fws_exclusion_cutoff"],
    message:
        "[moi] Writing high-MOI exclusion list"
    shell:
        r"""
        mkdir -p $(dirname {output.exclude})
        awk -F'\t' 'NR>1 && $2 < {params.cutoff} {{print $1}}' \
            {input.fws} > {output.exclude} 2> {log}
        echo "High-MOI samples: $(wc -l < {output.exclude})" >> {log}
        """


if HAS_COUNTRY or HAS_GEOGRAPHY:

    rule fws_summary:
        """
        Group Fws by the country role and/or the geography role and emit
        a TSV per role that's present. Falls silent for any absent role
        with a logged note — does NOT crash.

        WHAT: scripts/R/fws_summary.R reads canonical role columns
              (country / geography) from outputs/metadata/samples.tsv.
        WHY:  Per-region polyclonality rates are the core slide-6 reading.
              Hardcoding `Country` / `State` would break for any cohort
              that uses different column names.
        TUNABLES: moi.fws_polyclonal_cutoff, metadata.roles.{country,geography}
        OUTPUT: {outputs}/moi/fws_by_country.tsv (if country role present)
                {outputs}/moi/fws_by_geography.tsv (if geography role present)
        TRY:    null out metadata.roles.country in your config — the country
                summary drops out of the DAG with a logged note; the
                geography summary continues if its role is still set.
        """
        input:
            fws      = f"{PATHS['outputs']}/moi/fws_MOI.tsv",
            metadata = rules.validate_metadata.output.tsv,
        output:
            country   = (f"{PATHS['outputs']}/moi/fws_by_country.tsv"
                         if HAS_COUNTRY else []),
            geography = (f"{PATHS['outputs']}/moi/fws_by_geography.tsv"
                         if HAS_GEOGRAPHY else []),
        log:
            f"{PATHS['logs']}/moi/fws_summary.log",
        params:
            cutoff       = config["moi"]["fws_polyclonal_cutoff"],
            country_arg  = (f"{PATHS['outputs']}/moi/fws_by_country.tsv"
                            if HAS_COUNTRY else "NULL"),
            geo_arg      = (f"{PATHS['outputs']}/moi/fws_by_geography.tsv"
                            if HAS_GEOGRAPHY else "NULL"),
            script       = str(_AGNOSTIC / "scripts" / "R" / "fws_summary.R"),
        message:
            "[moi] Summarising Fws by role(s)"
        shell:
            r"""
            mkdir -p $(dirname {input.fws})
            Rscript {params.script} \
                {input.fws} \
                {input.metadata} \
                {params.country_arg} \
                {params.geo_arg} \
                {params.cutoff} \
                > {log} 2>&1
            """


if HAS_COUNTRY or HAS_GEOGRAPHY:

    rule plot_fws_density:
        """
        Density-panel Fws plots driven by the country / geography roles.
        Panel A = density coloured by country role; Panel B = density
        coloured by geography role. Each panel is rendered only if its
        role exists.

        WHAT: scripts/R/plot_fws_density.R reads canonical role columns
              and a polyclonal-cutoff line.
        WHY:  The slide-6 figure. V1's hardcoded "filter to Country ==
              Indonesia" province panel is dropped — Panel B is "density
              by geography" across the whole cohort.
        TUNABLES: moi.fws_polyclonal_cutoff, metadata.roles.{country,geography}
        OUTPUT: {reports}/figures/fws_density_by_country.{png,svg}
                {reports}/figures/fws_density_by_geography.{png,svg}
                (each pair only when its role is present)
        TRY:    null out metadata.roles.geography — the geography panel
                vanishes from the DAG; only the country panel renders.
        """
        input:
            fws      = f"{PATHS['outputs']}/moi/fws_MOI.tsv",
            metadata = rules.validate_metadata.output.tsv,
        output:
            country_png   = (f"{PATHS['reports']}/figures/fws_density_by_country.png"
                             if HAS_COUNTRY else []),
            country_svg   = (f"{PATHS['reports']}/figures/fws_density_by_country.svg"
                             if HAS_COUNTRY else []),
            geography_png = (f"{PATHS['reports']}/figures/fws_density_by_geography.png"
                             if HAS_GEOGRAPHY else []),
            geography_svg = (f"{PATHS['reports']}/figures/fws_density_by_geography.svg"
                             if HAS_GEOGRAPHY else []),
        log:
            f"{PATHS['logs']}/moi/plot_fws_density.log",
        params:
            cutoff = config["moi"]["fws_polyclonal_cutoff"],
            country_png_arg = (f"{PATHS['reports']}/figures/fws_density_by_country.png"
                               if HAS_COUNTRY else "NULL"),
            country_svg_arg = (f"{PATHS['reports']}/figures/fws_density_by_country.svg"
                               if HAS_COUNTRY else "NULL"),
            geo_png_arg     = (f"{PATHS['reports']}/figures/fws_density_by_geography.png"
                               if HAS_GEOGRAPHY else "NULL"),
            geo_svg_arg     = (f"{PATHS['reports']}/figures/fws_density_by_geography.svg"
                               if HAS_GEOGRAPHY else "NULL"),
            fig_dir = f"{PATHS['reports']}/figures",
            script  = str(_AGNOSTIC / "scripts" / "R" / "plot_fws_density.R"),
        message:
            "[moi] Plotting Fws density panels"
        shell:
            r"""
            mkdir -p {params.fig_dir}
            Rscript {params.script} \
                {input.fws} \
                {input.metadata} \
                {params.country_png_arg} \
                {params.country_svg_arg} \
                {params.geo_png_arg} \
                {params.geo_svg_arg} \
                {params.cutoff} \
                > {log} 2>&1
            """
