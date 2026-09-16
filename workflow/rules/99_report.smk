# =============================================================================
# Stage 7 — the run report.
#
# One HTML document that assembles whatever the pipeline actually produced for
# this cohort. It is NOT a narrative for any particular dataset: every heading,
# group name, cluster label and column name in it is read from the active
# config and from the output files on disk at render time. A stage that did not
# run is reported as such, with the reason, rather than omitted silently.
#
# It is a NAMED target rather than part of `rule all`:
#
#     snakemake report --configfile <cfg>
#
# (target first: --configfile takes one or more files, so a target written
# after it is swallowed as a second config file.)
#
# Two reasons. The report is the last thing you want, so it depends on
# FINAL_TARGETS and running it pulls the whole pipeline anyway — but keeping it
# out of `all` means an existing cohort's run does not suddenly acquire a
# hard dependency on quarto, and re-rendering the report after an edit does not
# look like the pipeline is out of date.
#
# FINAL_TARGETS is itself config-conditional (roles, stages, sample sets), so
# depending on it is what makes "assemble whatever ran" correct rather than
# best-effort: by the time the document renders, everything this config asked
# for exists.
# =============================================================================

_REPORT_QMD = str(_AGNOSTIC / "workflow" / "report" / "REPORT.qmd")


rule render_report:
    """
    Render the cohort-agnostic HTML report.

    WHAT: quarto render workflow/report/REPORT.qmd -P config:<cfg> -P root:<repo>
    WHY:  A run produces ~50 tables and ~25 figures across six stages and two
          sample-set arms. The report is the single artefact you can hand to a
          collaborator: it puts each result next to what it means, shows the
          full and de-clonalized arms side by side wherever the {sampleset}
          axis applies, and says plainly which stages did not run and why.
    TUNABLES: everything — the document is driven entirely by the active
          config. cohort.name titles it; metadata.roles decide which
          role-driven sections exist; clonality.declonalize decides whether
          there are two arms to compare; introgression.focal_group and
          selection.models switch their own sections on.
    OUTPUT: {reports}/report.html  (self-contained — figures are embedded, so
          the single file can be emailed or archived on its own)
    TRY:    unset a role in your config (say metadata.roles.date) and
            re-render — the clonal timeline section turns into a note
            explaining what would enable it, rather than disappearing.
    """
    input:
        # Unnamed and splatted, exactly as `rule all` takes them: FINAL_TARGETS
        # holds input *functions* alongside plain paths (the per-cluster and
        # per-pair fan-outs), and snakemake only resolves those when they are
        # unnamed positional inputs.
        _REPORT_QMD,
        *FINAL_TARGETS,
    output:
        html = f"{PATHS['reports']}/report.html",
    params:
        qmd  = _REPORT_QMD,
        cfg  = ACTIVE_CONFIG,
        root = str(_AGNOSTIC),
    log:
        f"{PATHS['logs']}/report/render_report.log",
    message:
        "[report] rendering {output.html}"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.html})" "$(dirname {log})"

        # Render in a scratch directory. Quarto writes its intermediates beside
        # the .qmd, and the .qmd lives in the tracked workflow tree — rendering
        # in place would litter it, and two cohorts rendering at once would
        # collide over the same intermediate names.
        work="$(mktemp -d)"
        trap 'rm -rf "$work"' EXIT
        cp {params.qmd} "$work/report.qmd"

        quarto render "$work/report.qmd" \
            --to html \
            -P config:"{params.cfg}" \
            -P root:"{params.root}" \
            > {log} 2>&1

        mv "$work/report.html" {output.html}
        """
