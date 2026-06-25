# qc_wgs.smk — WGS-specific final QC + MAF/missingness panels
# --------------------------------------------------------------------------
# Included by 01_qc.smk when cohort.input_type == "wgs". Produces the seam
# output (snps.qc.vcf.gz), variant_count.txt, and the four QC density panels.
# Logic is V1-faithful: biallelic SNPs + PASS + MAC, then PLINK --freq /
# --missing for the diagnostic plots.
# --------------------------------------------------------------------------


rule biallelic_pass_mac:
    """
    Keep biallelic SNPs with FILTER=PASS and MAC ≥ qc.min_mac. Writes the
    Stage-1 seam output that Stage 2+ reads (snps.qc.vcf.gz).

    WHAT: bcftools view -m2 -M2 -v snps -f PASS | bcftools view -e 'MAC<N'
    WHY:  Downstream stages (MOI/Fws, ADMIXTURE, PCA) assume biallelic
          diploid SNPs. MAC ≥ 2 removes singletons (mostly sequencing error).
    TUNABLES: qc.min_mac, qc.filter_pass
    OUTPUT: {outputs}/qc/snps.qc.vcf.gz   (the seam — see top of 01_qc.smk)
    TRY:    bump qc.min_mac to 5 to look at only common variants; or set
            qc.filter_pass: false to keep non-PASS calls and re-run from
            here to see how much sequencing noise the PASS filter removes.
    """
    input:
        vcf = QC_INPUT_VCF,
    output:
        vcf = f"{PATHS['outputs']}/qc/snps.qc.vcf.gz",
        idx = f"{PATHS['outputs']}/qc/snps.qc.vcf.gz.csi",
    log:
        f"{PATHS['logs']}/qc/biallelic_pass_mac.log",
    params:
        min_mac     = QC["min_mac"],
        pass_filter = "-f PASS" if QC.get("filter_pass", True) else "",
    threads: config["compute"]["threads_heavy"]
    message:
        "[qc:wgs] Biallelic + PASS + MAC filter"
    shell:
        r"""
        bcftools view \
            --threads {threads} \
            -m2 -M2 -v snps {params.pass_filter} \
            {input.vcf} \
          | bcftools view -e 'MAC<{params.min_mac}' -Oz -o {output.vcf} \
            2> {log}
        bcftools index --threads {threads} -c {output.vcf}
        echo "Final variant count: $(bcftools view -H {output.vcf} | wc -l)" >> {log}
        """


rule count_biallelic_variants:
    """
    Single-line variant count for the Stage-1 seam output.

    WHAT: bcftools view -H | wc -l
    WHY:  The report renders the post-Stage-1 count without re-running
          bcftools at knit time.
    TUNABLES: (none)
    OUTPUT: {outputs}/qc/variant_count.txt
    TRY:    after changing qc.min_mac, diff this number against an earlier
            run to see the impact of the MAC cutoff in one integer.
    """
    input:
        vcf = rules.biallelic_pass_mac.output.vcf,
    output:
        txt = f"{PATHS['outputs']}/qc/variant_count.txt",
    log:
        f"{PATHS['logs']}/qc/count_biallelic_variants.log",
    message:
        "[qc:wgs] Counting variants"
    shell:
        r"""
        bcftools view -H {input.vcf} | wc -l | awk '{{print $1}}' > {output.txt} 2> {log}
        """


rule compute_maf_missingness:
    """
    Compute MAF + per-sample/per-variant missingness with PLINK 1.9.

    WHAT: plink --vcf <qc.vcf.gz> --double-id --allow-extra-chr --freq --missing
    WHY:  Cheap diagnostic numbers + density panels that surface problem
          contigs or low-coverage samples before later stages run.
    TUNABLES: (none — these are diagnostic, not filtering, rules)
    OUTPUT: {outputs}/qc/plink/Pk.{{frq,imiss,lmiss}}
    TRY:    after a config change, eyeball maf_density_overall.png — a kink
            on the low end usually means the MAC cutoff should be higher.
    """
    input:
        vcf = rules.biallelic_pass_mac.output.vcf,
    output:
        frq    = f"{PATHS['outputs']}/qc/plink/Pk.frq",
        imiss  = f"{PATHS['outputs']}/qc/plink/Pk.imiss",
        lmiss  = f"{PATHS['outputs']}/qc/plink/Pk.lmiss",
    log:
        f"{PATHS['logs']}/qc/compute_maf_missingness.log",
    params:
        out_prefix = f"{PATHS['outputs']}/qc/plink/Pk",
    threads: config["compute"]["threads_heavy"]
    message:
        "[qc:wgs] PLINK --freq + --missing"
    shell:
        r"""
        mkdir -p $(dirname {params.out_prefix})
        plink \
            --vcf {input.vcf} \
            --double-id \
            --allow-extra-chr \
            --freq \
            --missing \
            --out {params.out_prefix} \
            --threads {threads} \
            > {log} 2>&1
        """


rule plot_qc_distributions:
    """
    Render the four QC density panels (MAF overall + by chrom; missingness
    per sample + per variant) used by the Quarto chapter.

    WHAT: scripts/R/plot_qc_distributions.R <frq> <imiss> <lmiss> <4 png paths>
    WHY:  Lets a reader eyeball the QC distributions before trusting any
          downstream number.
    TUNABLES: (none in the rule; cosmetic tweaks live in the R script)
    OUTPUT: {reports}/figures/maf_density_overall.png + 3 others
    TRY:    open maf_density_by_chrom.png — if one contig sits noticeably
            higher than the rest, that's a candidate for reference.exclude_contigs.
    """
    input:
        frq   = rules.compute_maf_missingness.output.frq,
        imiss = rules.compute_maf_missingness.output.imiss,
        lmiss = rules.compute_maf_missingness.output.lmiss,
    output:
        maf_overall      = f"{PATHS['reports']}/figures/maf_density_overall.png",
        maf_by_chrom     = f"{PATHS['reports']}/figures/maf_density_by_chrom.png",
        miss_per_sample  = f"{PATHS['reports']}/figures/miss_density_per_sample.png",
        miss_per_variant = f"{PATHS['reports']}/figures/miss_density_per_variant.png",
    log:
        f"{PATHS['logs']}/qc/plot_qc_distributions.log",
    params:
        script = str(_AGNOSTIC / "scripts" / "R" / "plot_qc_distributions.R"),
    message:
        "[qc:wgs] Plotting MAF + missingness density"
    shell:
        r"""
        mkdir -p $(dirname {output.maf_overall})
        Rscript {params.script} \
            {input.frq} {input.imiss} {input.lmiss} \
            {output.maf_overall} {output.maf_by_chrom} \
            {output.miss_per_sample} {output.miss_per_variant} \
            > {log} 2>&1
        """
