# 00_setup.smk — preflight rules (agnostic)
# --------------------------------------------------------------------------
# Three preflight steps every downstream rule depends on:
#   1. Extract the VCF sample list (bcftools query -l).
#   2. Persist the nuclear-contig list + contig→integer map to outputs/setup/
#      so they're inspectable and reusable by later stages.
#   3. Validate the user's samples.tsv against the configured role map,
#      check VCF concordance, and emit a canonical role-renamed copy.
#
# V1's six-file `build_metadata.R` is intentionally NOT ported (DESIGN.md §3b).
# The agnostic contract is one tidy samples.tsv supplied by the user.
# --------------------------------------------------------------------------

rule index_vcf:
    """
    Build the .csi index for the cohort VCF if it isn't already there.

    WHAT: bcftools index -c <vcf>
    WHY:  bcftools view -r <region> (used by subset_nuclear) requires the
          input VCF to be indexed. Without this rule, a cohort supplied
          without a sitting .csi fails cryptically at subset_nuclear. With
          it, the index is a tracked dependency that Snakemake will rebuild
          if the VCF mtime moves past the index's.
    TUNABLES: cohort.vcf
    OUTPUT: <cohort.vcf>.csi
    TRY:    `touch -d "1 hour ago" <cohort.vcf>` then re-run — the rule
            re-builds the index because its declared input is now newer
            than its declared output.
    """
    input:
        vcf = COHORT["vcf"],
    output:
        idx = COHORT["vcf"] + ".csi",
    log:
        f"{PATHS['logs']}/setup/index_vcf.log",
    threads: config["compute"]["threads_light"]
    message:
        "[setup] Indexing cohort VCF (.csi)"
    shell:
        r"""
        mkdir -p $(dirname {log})
        bcftools index --threads {threads} -c -f {input.vcf} 2> {log}
        """


rule extract_vcf_samples:
    """
    Extract one VCF sample ID per line; downstream rules read this rather
    than re-opening the (multi-GB) VCF.

    WHAT: bcftools query -l <vcf>
    WHY:  Every QC/metadata step needs the sample roster; doing it once,
          here, avoids redundant VCF reads.
    TUNABLES: cohort.vcf
    OUTPUT: {outputs}/setup/vcf_samples.txt
    TRY:    swap cohort.vcf to a subsetted VCF (bcftools view -s s1,s2,...)
            to dry-run the pipeline on a sample subset without re-deriving anything.
    """
    input:
        vcf = COHORT["vcf"],
    output:
        samples = f"{PATHS['outputs']}/setup/vcf_samples.txt",
    log:
        f"{PATHS['logs']}/setup/extract_vcf_samples.log",
    message:
        "[setup] Extracting VCF sample list"
    shell:
        r"""
        mkdir -p $(dirname {output.samples})
        bcftools query -l {input.vcf} > {output.samples} 2> {log}
        echo "Sample count: $(wc -l < {output.samples})" >> {log}
        """


rule build_contig_map:
    """
    Persist the nuclear-contig list and contig→integer map derived from the
    reference .fai (minus reference.exclude_contigs).

    WHAT: scripts/py/contigs_from_fai.py <fasta>.fai --exclude <patterns>
    WHY:  Replaces V1's hand-listed nuclear_contigs array. The integer map
          is the input later stages need for PLINK chromosome codes.
    TUNABLES: reference.fasta, reference.exclude_contigs
    OUTPUT: {outputs}/setup/nuclear_contigs.txt, {outputs}/setup/contig_map.tsv
    TRY:    add a contig regex to reference.exclude_contigs (e.g. "_14_") and
            re-run from here to see Stage 1 drop that contig from the VCF.
    """
    input:
        fai = REF["fasta"] + ".fai",
    output:
        contigs = f"{PATHS['outputs']}/setup/nuclear_contigs.txt",
        cmap    = f"{PATHS['outputs']}/setup/contig_map.tsv",
    log:
        f"{PATHS['logs']}/setup/build_contig_map.log",
    params:
        exclude_args = " ".join(f"--exclude {e}" for e in REF.get("exclude_contigs", [])),
        script       = str(_AGNOSTIC / "scripts" / "py" / "contigs_from_fai.py"),
    message:
        "[setup] Deriving nuclear contigs from .fai"
    shell:
        r"""
        mkdir -p $(dirname {output.contigs})
        python {params.script} {input.fai} \
            {params.exclude_args} \
            --out-list {output.contigs} \
            --out-map  {output.cmap} \
            > {log} 2>&1
        """


rule validate_metadata:
    """
    Validate samples.tsv against the configured role map + VCF concordance.
    Emit a canonical, role-renamed copy at {outputs}/metadata/samples.tsv
    so downstream rules read by role (sample_id, group, geography, ...) and
    never reach back into the cohort-specific column names.

    WHAT: scripts/R/validate_metadata.R <samples.tsv> <vcf_samples.txt>
          <roles.json> <out_tsv> <out_log>
    WHY:  Hard-errors only on missing sample_id; missing optional roles log
          a note and let dependent analyses skip gracefully (DESIGN §3b).
    TUNABLES: metadata.table, metadata.roles
    OUTPUT: {outputs}/metadata/samples.tsv
    TRY:    null out metadata.roles.geography in your config — the rule logs
            "geography role skipped" and downstream geography facets vanish
            without a crash.
    """
    input:
        meta    = META["table"],
        samples = rules.extract_vcf_samples.output.samples,
        roles   = str(ROLES_JSON),
    output:
        tsv     = f"{PATHS['outputs']}/metadata/samples.tsv",
        rpt_log = f"{PATHS['logs']}/setup/validate_metadata.report.log",
    log:
        f"{PATHS['logs']}/setup/validate_metadata.log",
    params:
        script = str(_AGNOSTIC / "scripts" / "R" / "validate_metadata.R"),
    message:
        "[setup] Validating metadata + VCF concordance"
    shell:
        r"""
        mkdir -p $(dirname {output.tsv})
        Rscript {params.script} \
            {input.meta} \
            {input.samples} \
            {input.roles} \
            {output.tsv} \
            {output.rpt_log} \
            > {log} 2>&1
        """
