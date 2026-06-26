# moi_wgs.smk — WGS-specific MOI/Fws path (biallelic moimix)
# --------------------------------------------------------------------------
# Included by 02_moi.smk when cohort.input_type == "wgs". Produces the seam
# output {outputs}/moi/fws_MOI.tsv consumed by the common Stage 2 rules and
# by Stage 3+.
#
# Chain:  snps.qc.vcf.gz  --MAF-->  snps.moi.vcf.gz  --SeqArray-->
#         snps.moi.gds    --moimix-->  fws_MOI.tsv  (+ BAF_dataframe.tsv)
# --------------------------------------------------------------------------


rule moi_filter_vcf:
    """
    Stricter MAF filter on the Stage-1 seam VCF before moimix.

    WHAT: bcftools +fill-tags -- -t AF | bcftools view -e 'MAF<min_maf'
    WHY:  moimix estimates BAF from read-depth ratios; low-MAF sites have
          unstable B-allele frequencies and add noise. WGS-specific (the
          microhap path doesn't use biallelic MAF), hence behind this seam.
    TUNABLES: moi.min_maf
    OUTPUT: {outputs}/moi/snps.moi.vcf.gz
    TRY:    raise moi.min_maf from 0.05 to 0.10 to keep only well-sampled
            variants; usually trims a small % of the Fws tail.
    """
    input:
        vcf = f"{PATHS['outputs']}/qc/snps.qc.vcf.gz",
    output:
        vcf = f"{PATHS['outputs']}/moi/snps.moi.vcf.gz",
        idx = f"{PATHS['outputs']}/moi/snps.moi.vcf.gz.csi",
    log:
        f"{PATHS['logs']}/moi/moi_filter_vcf.log",
    params:
        min_maf = config["moi"]["min_maf"],
    threads: config["compute"]["threads_heavy"]
    message:
        "[moi:wgs] MAF filter on QC'd VCF"
    shell:
        r"""
        mkdir -p $(dirname {output.vcf})
        bcftools +fill-tags {input.vcf} -- -t AF \
          | bcftools view --threads {threads} -e 'MAF<{params.min_maf}' \
            -Oz -o {output.vcf} 2> {log}
        bcftools index --threads {threads} -c {output.vcf}
        echo "Variants after MAF filter: $(bcftools view -H {output.vcf} | wc -l)" >> {log}
        """


rule vcf_to_gds:
    """
    Convert the moimix-filtered VCF to SeqArray GDS.

    WHAT: SeqArray::seqVCF2GDS(storage.option = "LZ4_RA")
    WHY:  GDS = random-access by variant + sample. moimix's BAF/Fws calls
          do per-sample passes over the variant matrix; tabix VCF is
          orders of magnitude slower for that pattern.
    TUNABLES: (none — storage option is fixed; threads from compute block)
    OUTPUT: {outputs}/moi/snps.moi.gds
    TRY:    run `Rscript scripts/R/vcf_to_gds.R <vcf> <out>` standalone on
            a subset VCF to see GDS conversion timing per million variants.
    """
    input:
        vcf = rules.moi_filter_vcf.output.vcf,
    output:
        gds = f"{PATHS['outputs']}/moi/snps.moi.gds",
    log:
        f"{PATHS['logs']}/moi/vcf_to_gds.log",
    params:
        script = str(_AGNOSTIC / "scripts" / "R" / "vcf_to_gds.R"),
    threads: config["compute"]["threads_heavy"]
    message:
        "[moi:wgs] VCF -> GDS"
    shell:
        r"""
        mkdir -p $(dirname {output.gds})
        Rscript {params.script} {input.vcf} {output.gds} > {log} 2>&1
        """


rule run_moimix:
    """
    Compute per-sample BAF matrix + Fws via moimix. Writes the seam output.

    WHAT: bafMatrix(gds) + getFws(gds), with set.seed(moi.seed) for
          reproducibility.
    WHY:  Fws is the within-host complexity score — lower = more polyclonal.
          The seam contract requires this output to be (sample, Proportion)
          so Stage 3+ doesn't branch on input_type.
    TUNABLES: moi.seed
    OUTPUT: {outputs}/moi/fws_MOI.tsv         (the seam)
            {outputs}/moi/BAF_dataframe.tsv   (long-format BAF, for the
                                               Quarto narrative / QC eyes)
    TRY:    change moi.seed and compare fws_MOI.tsv pre/post — should be
            bit-identical, since moimix uses a deterministic estimator;
            divergence here would suggest a non-deterministic code path
            (none expected with the current moimix version).
    """
    input:
        gds = rules.vcf_to_gds.output.gds,
    output:
        fws = f"{PATHS['outputs']}/moi/fws_MOI.tsv",
        baf = f"{PATHS['outputs']}/moi/BAF_dataframe.tsv",
    log:
        f"{PATHS['logs']}/moi/run_moimix.log",
    params:
        seed   = config["moi"]["seed"],
        script = str(_AGNOSTIC / "scripts" / "R" / "run_moimix.R"),
    message:
        "[moi:wgs] moimix: BAF + Fws"
    shell:
        r"""
        Rscript {params.script} \
            {input.gds} \
            {output.fws} \
            {output.baf} \
            {params.seed} \
            > {log} 2>&1
        """
