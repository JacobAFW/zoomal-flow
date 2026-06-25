# 01_qc.smk — Stage 1 QC, common rules + WGS/microhap include seam
# --------------------------------------------------------------------------
# SEAM CONTRACT (DESIGN §4):
#   Both qc_wgs.smk and qc_microhap.smk MUST produce
#     {outputs}/qc/snps.qc.vcf.gz  (+ .csi)
#   as their terminal output. Stage 2+ reads exclusively from that path and
#   does not branch on input_type — the fork is contained here.
#
# Common rules below run regardless of input_type:
#   - filter_samples_list   drop controls (regex on VCF sample IDs)
#   - subset_nuclear        keep nuclear contigs (auto-derived from .fai)
#   - prepare_mask_bed      colon-dash region file → BED  (only if mask set)
#   - mask_regions          bcftools view -T ^bed         (only if mask set)
#
# After the mask step we include qc_wgs.smk OR qc_microhap.smk based on
# config.cohort.input_type. The WGS path is implemented; the microhap path
# is stubbed and fails fast.
# --------------------------------------------------------------------------

# Did the user set a mask file? Drives both prepare_mask_bed and the input
# of the input_type-specific seam: the WGS rule's input is the mask output
# if a mask was configured, otherwise the subset_nuclear output.
MASK_REGIONS = REF.get("mask_regions")
HAS_MASK     = bool(MASK_REGIONS)

# The exact file the input_type-specific seam reads from. Computed once here
# so both qc_wgs.smk and qc_microhap.smk can refer to it as QC_INPUT_VCF.
QC_INPUT_VCF = (
    f"{PATHS['outputs']}/qc/snps.masked.vcf.gz"
    if HAS_MASK
    else f"{PATHS['outputs']}/qc/snps.nuclear.vcf.gz"
)


rule filter_samples_list:
    """
    Build a keep-list by dropping anything in the VCF roster matching the
    controls.exclude_patterns regex (ctrl|cpos|cneg by default).

    WHAT: grep -viE '<patterns>' vcf_samples.txt
    WHY:  Controls in a pop-gen VCF break allele-frequency stats; dropping
          them up front means every downstream count is over the analysis
          set.
    TUNABLES: controls.exclude_patterns
    OUTPUT: {outputs}/qc/keep_samples.txt
    TRY:    add a cohort-specific pattern (e.g. "blank") to the list and
            re-run from here to see the sample count drop accordingly.
    """
    input:
        samples = rules.extract_vcf_samples.output.samples,
    output:
        keep = f"{PATHS['outputs']}/qc/keep_samples.txt",
    log:
        f"{PATHS['logs']}/qc/filter_samples_list.log",
    params:
        exclude_regex = "|".join(CONTROLS["exclude_patterns"]),
    message:
        "[qc] Filtering sample list (drop controls)"
    shell:
        r"""
        mkdir -p $(dirname {output.keep})
        if [ -z "{params.exclude_regex}" ]; then
            cp {input.samples} {output.keep}
        else
            grep -viE '{params.exclude_regex}' {input.samples} > {output.keep} 2> {log} || true
        fi
        echo "  Kept:     $(wc -l < {output.keep})"  >> {log}
        echo "  Excluded: $(($(wc -l < {input.samples}) - $(wc -l < {output.keep})))" >> {log}
        """


rule subset_nuclear:
    """
    Restrict the VCF to nuclear contigs + non-control samples.

    WHAT: bcftools view -S keep -r <contigs> ; contigs derived from
          <fasta>.fai minus reference.exclude_contigs at parse time.
    WHY:  Mitochondrial / apicoplast contigs are haploid and break
          diploid-SNP assumptions in MOI / structure / IBD downstream.
          bcftools `-r` requires a .csi/tbi beside the VCF; that index is
          declared as an explicit input below (built by rule index_vcf)
          so Snakemake tracks it and rebuilds it when stale, rather than
          relying on one being silently sitting on disk.
    TUNABLES: reference.exclude_contigs, controls.exclude_patterns
    OUTPUT: {outputs}/qc/snps.nuclear.vcf.gz
    TRY:    widen reference.exclude_contigs to drop a noisy contig and
            re-run from here to see the effect on the Stage-1 variant count.
    """
    input:
        vcf  = COHORT["vcf"],
        idx  = rules.index_vcf.output.idx,
        keep = rules.filter_samples_list.output.keep,
    output:
        vcf = f"{PATHS['outputs']}/qc/snps.nuclear.vcf.gz",
        idx = f"{PATHS['outputs']}/qc/snps.nuclear.vcf.gz.csi",
    log:
        f"{PATHS['logs']}/qc/subset_nuclear.log",
    params:
        contigs = ",".join(NUCLEAR_CONTIGS),
    threads: config["compute"]["threads_heavy"]
    message:
        "[qc] Subsetting to nuclear contigs + non-control samples"
    shell:
        r"""
        bcftools view \
            --threads {threads} \
            -S {input.keep} --force-samples \
            -r {params.contigs} \
            -Oz -o {output.vcf} \
            {input.vcf} \
            2> {log}
        bcftools index --threads {threads} -c {output.vcf}
        """


if HAS_MASK:

    rule prepare_mask_bed:
        """
        Convert a bcftools region-string list (chrom:start-end, 1-based
        inclusive) to BED (0-based half-open) so `bcftools view -T ^bed`
        can consume it. Sorted by chrom then start.

        WHAT: awk colon-dash parse → sort -k1,1 -k2,2n
        WHY:  V1's regions_to_mask.list uses bcftools `-r` syntax; `-T`
              wants BED. Convert once and reuse.
        TUNABLES: reference.mask_regions
        OUTPUT: {outputs}/qc/regions_to_mask.bed
        TRY:    drop a region from reference.mask_regions and rerun from
                here to see Stage 1's variant count climb a little.
        """
        input:
            mask = MASK_REGIONS,
        output:
            bed = f"{PATHS['outputs']}/qc/regions_to_mask.bed",
        log:
            f"{PATHS['logs']}/qc/prepare_mask_bed.log",
        message:
            "[qc] Converting mask list to BED"
        shell:
            r"""
            mkdir -p $(dirname {output.bed})
            awk -F'[:-]' '
                /^[^#]/ && NF >= 3 {{
                    printf "%s\t%d\t%s\n", $1, $2 - 1, $3
                }}
            ' {input.mask} \
              | sort -k1,1 -k2,2n > {output.bed} 2> {log}
            echo "BED regions written: $(wc -l < {output.bed})" >> {log}
            """

    rule mask_regions:
        """
        Exclude variants in masked regions (telomeric repeats, subtelomeric
        gene families, etc.).

        WHAT: bcftools view -T ^<bed>
        WHY:  These regions accumulate paralog reads and false variants.
              Masking before any allele-frequency or MOI calc avoids
              biasing downstream estimates.
        TUNABLES: reference.mask_regions (set to null in config to skip
                  this rule entirely).
        OUTPUT: {outputs}/qc/snps.masked.vcf.gz
        TRY:    set reference.mask_regions: null in the config and dry-run
                — the prepare_mask_bed + mask_regions rules disappear from
                the DAG and the input_type-specific QC step pulls directly
                from snps.nuclear.vcf.gz.
        """
        input:
            vcf = rules.subset_nuclear.output.vcf,
            bed = rules.prepare_mask_bed.output.bed,
        output:
            vcf = f"{PATHS['outputs']}/qc/snps.masked.vcf.gz",
            idx = f"{PATHS['outputs']}/qc/snps.masked.vcf.gz.csi",
        log:
            f"{PATHS['logs']}/qc/mask_regions.log",
        threads: config["compute"]["threads_heavy"]
        message:
            "[qc] Excluding masked regions"
        shell:
            r"""
            bcftools view \
                --threads {threads} \
                -T ^{input.bed} \
                -Oz -o {output.vcf} \
                {input.vcf} \
                2> {log}
            bcftools index --threads {threads} -c {output.vcf}
            """


# --------------------------------------------------------------------------
# input_type seam — pick the variant-filter implementation
# --------------------------------------------------------------------------
_INPUT_TYPE = COHORT["input_type"]
if _INPUT_TYPE == "wgs":
    include: "qc_wgs.smk"
elif _INPUT_TYPE == "microhap":
    include: "qc_microhap.smk"
else:
    # Should be unreachable because the config schema enforces the enum.
    raise ValueError(
        f"Unknown cohort.input_type={_INPUT_TYPE!r}; expected 'wgs' or 'microhap'."
    )
