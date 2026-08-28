# 04b_clonality.smk — de-clonalization: the {sampleset} keep-lists
# --------------------------------------------------------------------------
# Clonal samples are near-identical genotypes. Counting each as independent
# evidence is pseudo-replication, and it biases every FREQUENCY- or COUNT-based
# analysis: allele frequencies, ADMIXTURE ancestry proportions, per-group
# polyclonality rates, introgression window support, iHS case/control
# frequencies. Plasmodium cohorts are frequently clonal, so an agnostic
# pop-gen pipeline that ignores this is scientifically incomplete.
#
# THE CIRCULARITY, AND HOW IT IS HANDLED. Clonality is DETECTED from the data
# (Stage 4 hmmIBD), and detection needs every sample — you cannot find clonal
# groups in a set you have already thinned. So:
#
#   1. the FULL set runs through Stage 4 exactly as before, and the full-set
#      clusters remain the operational backbone (Stage 4's own clonal
#      detection and the Stage-5 introgression cluster definitions key off
#      them — nothing in this file changes that);
#   2. AFTER Stage 4, this file derives the unique-genotypes keep-list;
#   3. the frequency stages then run a SECOND, ONE-PASS branch on it.
#
# One pass, deliberately. Clonality is NOT re-derived on the de-clonalized
# clusters — that would be circular, and clonal groups are robust to it
# (field norm). See docs/clonality.md.
#
# THIS FILE PRODUCES ONLY KEEP-LISTS. Everything else is the `{sampleset}`
# wildcard on the frequency-stage rules; the code path is identical bar which
# keep-list is handed to PLINK / the R scripts.
# --------------------------------------------------------------------------

CLONALITY_DIR = f"{PATHS['outputs']}/clonality"


rule all_genotypes:
    """
    The FULL sample set, as a keep-list.

    WHAT: the sample column of the shared PLINK .fam, sorted and de-duplicated.
    WHY:  so `full` and `unique` are the same KIND of object and the frequency
          stages are one code path with a swapped file, not two branches.
          Crucially this list is derived from the SEAM, not from Stage 4 —
          `full` structure feeds Stage 4, Stage 4 feeds the unique keep-list,
          so a full keep-list that depended on Stage 4 would close a cycle in
          the DAG. That asymmetry in the inputs IS the architecture: `full`
          comes before clonality, `unique` comes after it.
    OUTPUT: {outputs}/clonality/all_genotypes.txt (readable list of the set),
            {outputs}/clonality/exclude_full.txt (empty — the wired artifact)
    """
    input:
        fam = f"{PATHS['outputs']}/structure/Pk.fam",
    output:
        txt     = f"{CLONALITY_DIR}/all_genotypes.txt",
        exclude = f"{CLONALITY_DIR}/exclude_full.txt",
    log:
        f"{PATHS['logs']}/clonality/all_genotypes.log",
    message:
        "[clonality] Full sample keep-list"
    shell:
        r"""
        mkdir -p $(dirname {output.txt}) $(dirname {log})
        awk '{{print $1}}' {input.fam} | sort -u > {output.txt} 2> {log}
        # The `full` arm excludes nobody. An EMPTY drop-list is what makes
        # `full` a guaranteed no-op in every stage, whatever that stage's own
        # sample universe happens to be.
        : > {output.exclude}
        echo "full sample set: $(wc -l < {output.txt}) samples, 0 excluded" >> {log}
        """


rule unique_genotypes:
    """
    Collapse each clonal group to one representative; write both keep-lists.

    WHAT: scripts/R/unique_genotypes.R over Stage 4's clonal_clusters.tsv.
          Writes unique_genotypes.txt (one
          representative per clonal group + every sample in no group) and an
          audit table naming the representative and the dropped members of
          every group. The full-set list is the sibling `all_genotypes` rule.
    WHY:  a clonal block of three in a 35-sample cluster is one genotype
          speaking three times. Any frequency- or count-based result computed
          over it is inflated by construction. Making the de-clonalized set a
          first-class, always-built artifact means every such claim can be
          checked rather than trusted.
    REPRESENTATIVE RULE: lowest per-sample missingness (PLINK .smiss F_MISS),
          ties broken alphabetically — the same question, and the same answer,
          as the Stage-3 duplicate pick in find_duplicates.R. Deterministic:
          re-running picks the same sample.
    STRADDLING GROUPS: a clonal group spanning more than one ADMIXTURE cluster
          takes its representative from the majority cluster; cluster ties go
          to the alphabetically-first name. Counted and logged, never silent.
    NOTE: if there is no clonal signal at all, unique_genotypes.txt equals
          all_genotypes.txt and the `unique` branch simply reproduces `full` —
          a legitimate result, not an error.
    TUNABLES: clonality.declonalize (gates the whole branch),
              ibd.clonal_ibd_threshold (upstream — what counts as clonal)
    OUTPUT: {outputs}/clonality/unique_genotypes.txt (readable list of the set),
            {outputs}/clonality/exclude_unique.txt (the dropped members — the
                    artifact the frequency stages actually consume),
            {outputs}/clonality/declonalization_audit.tsv
    TRY:    `wc -l` the two keep-lists — the difference is how much
            pseudo-replication the full-set results carry. Then read the
            per-stage declonalization_comparison.tsv to see what it moved.
    """
    input:
        clonal = f"{PATHS['outputs']}/ibd/clonal_clusters.tsv",
        # Sample-set-independent seam artifacts: the .fam defines the universe,
        # the .smiss ranks candidates inside a clonal group.
        fam    = f"{PATHS['outputs']}/structure/Pk.fam",
        smiss  = f"{PATHS['outputs']}/structure/Pk.smiss",
    output:
        unique  = f"{CLONALITY_DIR}/unique_genotypes.txt",
        exclude = f"{CLONALITY_DIR}/exclude_unique.txt",
        audit   = f"{CLONALITY_DIR}/declonalization_audit.tsv",
    log:
        f"{PATHS['logs']}/clonality/unique_genotypes.log",
    params:
        script = str(_AGNOSTIC / "scripts" / "R" / "unique_genotypes.R"),
    message:
        "[clonality] Deriving the unique-genotypes keep-list"
    shell:
        r"""
        mkdir -p $(dirname {output.unique}) $(dirname {log})
        Rscript {params.script} \
            --clonal-clusters {input.clonal} \
            --all-samples     {input.fam} \
            --smiss           {input.smiss} \
            --out-unique      {output.unique} \
            --out-exclude     {output.exclude} \
            --out-audit       {output.audit} \
            > {log} 2>&1
        """
