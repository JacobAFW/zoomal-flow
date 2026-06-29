# structure_prep_wgs.smk — WGS structure-input prep
# --------------------------------------------------------------------------
# Included by 03_structure.smk when cohort.input_type == "wgs". Produces
# the seam outputs (cleaned.{bed,bim,fam} + cleaned.ld.{bed,bim,fam}).
#
# Chain:
#   snps.qc.vcf.gz
#     → exclude_high_moi          (drop Stage 2's exclude list)
#     → normalise_vcf             (split-multi + left-align + set ID)
#     → build_variant_chrom_map   (variant ID → integer code via contig_map)
#     → vcf_to_plink              (plink2 --update-chr --sort-vars --make-bed)
#     → plink_freq_missing        (afreq / smiss)
#     → find_duplicates           (R, role/pattern-driven)
#     → final_filters             (V→S→V→M chain, legacy order)
#     → ld_prune                  (plink2 --indep-pairwise 50 5 0.5)
#
# The variant-ID → integer map is built by joining the bcftools-emitted
# CHROM column to outputs/setup/contig_map.tsv (Stage 0). This replaces
# V1's hardcoded `gsub("ordered_PKNH_","")/gsub("_v2","")` — that worked
# only for the PKA1H1 reference; the agnostic version reuses the Stage-0
# derivation.
# --------------------------------------------------------------------------


rule exclude_high_moi:
    """
    Drop high-MOI samples (from Stage 2) from the Stage-1 seam VCF.

    WHAT: bcftools view -S ^exclude --force-samples
    WHY:  Polyclonal samples break single-genotype assumptions in
          ADMIXTURE / PCA / NJ. Filter once, here, at the entry of Stage 3.
    TUNABLES: moi.fws_exclusion_cutoff (chooses the exclusion-list contents)
    OUTPUT: {outputs}/structure/snps.moi_excluded.vcf.gz
    TRY:    after raising fws_exclusion_cutoff, count the remaining samples
            and watch Stage-3's final cohort drop.
    """
    input:
        vcf     = f"{PATHS['outputs']}/qc/snps.qc.vcf.gz",
        exclude = f"{PATHS['outputs']}/moi/exclude_high_moi.txt",
    output:
        vcf = f"{PATHS['outputs']}/structure/snps.moi_excluded.vcf.gz",
        idx = f"{PATHS['outputs']}/structure/snps.moi_excluded.vcf.gz.csi",
    log:
        f"{PATHS['logs']}/structure/exclude_high_moi.log",
    threads: config["compute"]["threads_heavy"]
    message:
        "[structure:wgs] Excluding high-MOI samples"
    shell:
        r"""
        mkdir -p $(dirname {output.vcf})
        bcftools view --threads {threads} \
            -S ^{input.exclude} --force-samples \
            -Oz -o {output.vcf} {input.vcf} 2> {log}
        bcftools index --threads {threads} -c {output.vcf}
        echo "Samples after MOI exclusion: $(bcftools query -l {output.vcf} | wc -l)" >> {log}
        """


rule normalise_vcf:
    """
    Split multiallelic records, left-align, set ID = CHROM:POS:REF:ALT.

    WHAT: bcftools norm -m-any | bcftools norm --check-ref w -f <ref>
          | bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT'
    WHY:  PLINK needs unique variant IDs and biallelic records; ADMIXTURE
          is strict about both. Doing it here, once, keeps every downstream
          tool happy without re-normalisation.
    TUNABLES: reference.fasta (used for left-alignment ref check)
    OUTPUT: {outputs}/structure/snps.normalised.vcf.gz
    TRY:    head the .vcf.gz — IDs should look like
            `ordered_PKNH_01_v2:12345:A:C` (contig:pos:ref:alt).
    """
    input:
        vcf = rules.exclude_high_moi.output.vcf,
        ref = REF["fasta"],
    output:
        vcf = f"{PATHS['outputs']}/structure/snps.normalised.vcf.gz",
        idx = f"{PATHS['outputs']}/structure/snps.normalised.vcf.gz.csi",
    log:
        f"{PATHS['logs']}/structure/normalise_vcf.log",
    threads: config["compute"]["threads_heavy"]
    message:
        "[structure:wgs] Normalising VCF (split + left-align + ID)"
    shell:
        r"""
        bcftools norm --threads {threads} -m-any {input.vcf} 2>> {log} \
          | bcftools norm --threads {threads} --check-ref w -f {input.ref} 2>> {log} \
          | bcftools annotate --threads {threads} \
                -x ID -I +'%CHROM:%POS:%REF:%ALT' \
                -Oz -o {output.vcf} 2>> {log}
        bcftools index --threads {threads} -c {output.vcf}
        """


rule build_variant_chrom_map:
    """
    Build the variant-ID → integer-chrom-code map for PLINK by joining the
    normalised VCF's CHROM column to {outputs}/setup/contig_map.tsv.

    WHAT: bcftools query -f '%ID\t%CHROM\n' | awk join against contig_map
    WHY:  PLINK / ADMIXTURE want integer chromosome codes. Stage 0 already
          derived (contig_name → integer) from the reference .fai (in fai
          order); we reuse that single source of truth instead of V1's
          PKA1H1-specific `gsub("ordered_PKNH_",""); gsub("_v2","")`.
    TUNABLES: (none — driven by reference.fasta + reference.exclude_contigs
              upstream)
    OUTPUT: {outputs}/structure/chrom_update.txt   (variant_ID  int_code)
    TRY:    diff this against V1's outputs/structure/chrom_update.txt (in
            agnostic dev) — the integer codes should match the .fai order.
    """
    input:
        vcf   = rules.normalise_vcf.output.vcf,
        cmap  = f"{PATHS['outputs']}/setup/contig_map.tsv",
    output:
        tsv = f"{PATHS['outputs']}/structure/chrom_update.txt",
    log:
        f"{PATHS['logs']}/structure/build_variant_chrom_map.log",
    message:
        "[structure:wgs] Building variant→int-chrom map (from contig_map.tsv)"
    shell:
        r"""
        # awk join: contig_map.tsv is (contig\tint); we stream variant
        # records (ID\tCHROM) and substitute CHROM→int via the map.
        bcftools query -f '%ID\t%CHROM\n' {input.vcf} \
          | awk -v cmap={input.cmap} '
              BEGIN{{
                while ((getline line < cmap) > 0) {{
                    n = split(line, a, "\t")
                    if (n >= 2) map[a[1]] = a[2]
                }}
                close(cmap)
                OFS = "\t"
              }}
              {{
                if (!($2 in map)) {{
                    print "WARN: contig not in map: " $2 > "/dev/stderr"
                    next
                }}
                print $1, map[$2]
              }}
            ' > {output.tsv} 2> {log}
        echo "Chrom update rows: $(wc -l < {output.tsv})" >> {log}
        """


rule vcf_to_plink:
    """
    Convert the normalised VCF to PLINK (.bed/.bim/.fam) with integer
    chromosome codes and sample IDs preserved verbatim.

    WHAT: plink2 --vcf … --double-id --update-chr <map> --sort-vars
          --make-pgen → plink2 --pfile … --make-bed
          (two-step because plink2's --make-bed with --update-chr +
          --sort-vars requires a sorted .pgen intermediate)
    WHY:  ADMIXTURE / downstream PLINK rules want a clean integer-coded
          bfile with unique sample IDs.
    TUNABLES: (none in the rule)
    OUTPUT: {outputs}/structure/Pk.{bed,bim,fam}
    TRY:    head Pk.bim — chrom column should be 1..N integers, no contig
            strings; sample IDs in Pk.fam should match VCF sample IDs verbatim.
    """
    input:
        vcf   = rules.normalise_vcf.output.vcf,
        chrom = rules.build_variant_chrom_map.output.tsv,
    output:
        bed = f"{PATHS['outputs']}/structure/Pk.bed",
        bim = f"{PATHS['outputs']}/structure/Pk.bim",
        fam = f"{PATHS['outputs']}/structure/Pk.fam",
    log:
        f"{PATHS['logs']}/structure/vcf_to_plink.log",
    threads: config["compute"]["threads_heavy"]
    params:
        prefix     = f"{PATHS['outputs']}/structure/Pk",
        tmp_prefix = f"{PATHS['outputs']}/structure/Pk.tmp",
    message:
        "[structure:wgs] VCF → PLINK (integer chrom codes)"
    shell:
        r"""
        plink2 --vcf {input.vcf} \
            --double-id \
            --allow-extra-chr 0 \
            --update-chr {input.chrom} 2 1 \
            --sort-vars \
            --make-pgen \
            --threads {threads} \
            --out {params.tmp_prefix} > {log} 2>&1
        plink2 --pfile {params.tmp_prefix} \
            --allow-extra-chr 0 \
            --make-bed \
            --threads {threads} \
            --out {params.prefix} >> {log} 2>&1
        rm -f {params.tmp_prefix}.pgen {params.tmp_prefix}.pvar \
              {params.tmp_prefix}.psam {params.tmp_prefix}.log
        """


rule plink_freq_missing:
    """
    Initial --freq + --missing on the full PLINK set — drives duplicate
    detection (lowest-missingness replicate) + the variant-filter chain.

    WHAT: plink2 --bfile Pk --freq --missing
    WHY:  Per-sample F_MISS is what dedup needs; allele frequencies feed
          the final --maf filter further down.
    TUNABLES: (none — diagnostic, not filtering)
    OUTPUT: Pk.afreq, Pk.smiss, Pk.vmiss
    TRY:    eyeball Pk.smiss F_MISS distribution; any extreme outliers
            usually flag a sequencing failure.
    """
    input:
        bed = rules.vcf_to_plink.output.bed,
    output:
        smiss = f"{PATHS['outputs']}/structure/Pk.smiss",
        afreq = f"{PATHS['outputs']}/structure/Pk.afreq",
    log:
        f"{PATHS['logs']}/structure/plink_freq_missing.log",
    params:
        prefix = f"{PATHS['outputs']}/structure/Pk",
    threads: config["compute"]["threads_heavy"]
    message:
        "[structure:wgs] PLINK --freq --missing"
    shell:
        r"""
        plink2 --bfile {params.prefix} \
            --allow-extra-chr \
            --freq --missing \
            --threads {threads} \
            --out {params.prefix} > {log} 2>&1
        """


rule find_duplicates:
    """
    Identify replicate samples (same biological isolate sequenced twice)
    and write the IDs of the higher-missingness replicates to a PLINK
    --remove file.

    WHAT: scripts/R/find_duplicates.R reads .fam + .smiss, strips the
          configured duplicate_id_pattern (default: don't strip), then
          normalises out underscores/hyphens to derive a dedup base ID.
          Keeps the lowest-missing replicate per base ID.
    WHY:  V1 stripped Illumina lane suffix `_DK.*` which is PKA1H1-cohort-
          specific. The agnostic version takes the regex from config; null
          = no strip (still dedups on the underscore/hyphen normalisation).
    TUNABLES: structure.duplicate_id_pattern (regex; null = no strip)
    OUTPUT: {outputs}/structure/Pk.dups
    TRY:    set duplicate_id_pattern to ".*" — every sample becomes a dup
            of every other; rule will keep one. Useful sanity check that
            the strip is happening at all.
    """
    input:
        fam   = rules.vcf_to_plink.output.fam,
        smiss = rules.plink_freq_missing.output.smiss,
    output:
        dups = f"{PATHS['outputs']}/structure/Pk.dups",
    log:
        f"{PATHS['logs']}/structure/find_duplicates.log",
    params:
        script  = str(_AGNOSTIC / "scripts" / "R" / "find_duplicates.R"),
        pattern = "" if DUP_PATTERN is None else DUP_PATTERN,
    message:
        "[structure:wgs] Finding duplicate replicates"
    shell:
        r"""
        Rscript {params.script} \
            {input.fam} \
            {input.smiss} \
            {output.dups} \
            "{params.pattern}" \
            > {log} 2>&1
        """


rule final_filters:
    """
    Sequential 4-step filter chain (legacy order V→S→V→M):
      1. Remove duplicate replicates + lenient variant filter (--geno 0.20).
      2. Sample-missingness filter (--mind) on the cleaned variants.
      3. Stricter variant filter (--geno) on the post-sample set.
      4. MAF filter (--maf).

    WHAT: four plink2 passes, intermediates cleaned at the end.
    WHY:  V1 found that a strict sample-first order dropped 73% of the
          Indonesia cohort because the 1.4M-variant unfiltered set carries
          many low-coverage sites that drag down per-sample missingness.
          The V→S→V→M order matches HPC behaviour and recovers the cohort.
    TUNABLES: structure.max_sample_missing (step 2 --mind),
              structure.max_variant_missing (step 3 --geno),
              structure.min_maf (step 4 --maf)
    OUTPUT: {outputs}/structure/cleaned.{bed,bim,fam}
    TRY:    bump max_sample_missing to 0.05 and re-run — every step shows
            the cohort shrinking by tens.
    """
    input:
        bed  = rules.vcf_to_plink.output.bed,
        dups = rules.find_duplicates.output.dups,
    output:
        bed = f"{PATHS['outputs']}/structure/cleaned.bed",
        bim = f"{PATHS['outputs']}/structure/cleaned.bim",
        fam = f"{PATHS['outputs']}/structure/cleaned.fam",
    log:
        f"{PATHS['logs']}/structure/final_filters.log",
    params:
        prefix       = f"{PATHS['outputs']}/structure/Pk",
        out_prefix   = f"{PATHS['outputs']}/structure/cleaned",
        step1_prefix = f"{PATHS['outputs']}/structure/Pk.v1",
        step2_prefix = f"{PATHS['outputs']}/structure/Pk.v1s",
        step3_prefix = f"{PATHS['outputs']}/structure/Pk.v1sv2",
        mind         = STRUCTURE["max_sample_missing"],
        geno         = STRUCTURE["max_variant_missing"],
        maf          = STRUCTURE["min_maf"],
    threads: config["compute"]["threads_heavy"]
    message:
        "[structure:wgs] Final V→S→V→M filter chain"
    shell:
        r"""
        plink2 --bfile {params.prefix} \
            --allow-extra-chr \
            --remove {input.dups} \
            --geno 0.20 \
            --make-bed --threads {threads} \
            --out {params.step1_prefix} > {log} 2>&1
        plink2 --bfile {params.step1_prefix} \
            --allow-extra-chr \
            --mind {params.mind} \
            --make-bed --threads {threads} \
            --out {params.step2_prefix} >> {log} 2>&1
        plink2 --bfile {params.step2_prefix} \
            --allow-extra-chr \
            --geno {params.geno} \
            --make-bed --threads {threads} \
            --out {params.step3_prefix} >> {log} 2>&1
        plink2 --bfile {params.step3_prefix} \
            --allow-extra-chr \
            --maf {params.maf} \
            --make-bed --threads {threads} \
            --out {params.out_prefix} >> {log} 2>&1
        rm -f {params.step1_prefix}.bed {params.step1_prefix}.bim {params.step1_prefix}.fam {params.step1_prefix}.log \
              {params.step2_prefix}.bed {params.step2_prefix}.bim {params.step2_prefix}.fam {params.step2_prefix}.log {params.step2_prefix}.mindrem.id \
              {params.step3_prefix}.bed {params.step3_prefix}.bim {params.step3_prefix}.fam {params.step3_prefix}.log
        echo "Final cohort:" >> {log}
        echo "  samples:  $(wc -l < {params.out_prefix}.fam)" >> {log}
        echo "  variants: $(wc -l < {params.out_prefix}.bim)" >> {log}
        """


rule ld_prune:
    """
    LD-prune the cleaned bfile before ADMIXTURE.

    WHAT: plink2 --indep-pairwise 50 5 0.5 → plink --extract --make-bed
          (plink 1.9 for the extract step; plink2 alpha segfaults on this
          combination, both produce identical bfiles).
    WHY:  ADMIXTURE docs recommend unlinked SNPs. Without pruning, K=5
          alone took >3 h on a single core in V1; pruning takes it to
          ~minutes total. The standard 50/5/0.5 window typically retains
          5-10% of variants with full population-structure signal.
    TUNABLES: (none in the rule; window params hardcoded to PLINK defaults)
    OUTPUT: {outputs}/structure/cleaned.ld.{bed,bim,fam}
            + {outputs}/structure/cleaned.prune.in (kept-variant list)
    TRY:    after pruning, check `wc -l cleaned.prune.in` vs cleaned.bim —
            ~5-10% kept is typical for an outbred eukaryote at this scale.
    """
    input:
        bed = rules.final_filters.output.bed,
        bim = rules.final_filters.output.bim,
        fam = rules.final_filters.output.fam,
    output:
        prune_in = f"{PATHS['outputs']}/structure/cleaned.prune.in",
        bed      = f"{PATHS['outputs']}/structure/cleaned.ld.bed",
        bim      = f"{PATHS['outputs']}/structure/cleaned.ld.bim",
        fam      = f"{PATHS['outputs']}/structure/cleaned.ld.fam",
    log:
        f"{PATHS['logs']}/structure/ld_prune.log",
    params:
        in_prefix  = f"{PATHS['outputs']}/structure/cleaned",
        out_prefix = f"{PATHS['outputs']}/structure/cleaned.ld",
    threads: config["compute"]["threads_heavy"]
    message:
        "[structure:wgs] LD-pruning before ADMIXTURE"
    shell:
        r"""
        plink2 --bfile {params.in_prefix} \
            --allow-extra-chr \
            --indep-pairwise 50 5 0.5 \
            --threads {threads} \
            --out {params.in_prefix} > {log} 2>&1
        plink --bfile {params.in_prefix} \
            --allow-extra-chr \
            --extract {output.prune_in} \
            --make-bed --threads {threads} \
            --out {params.out_prefix} >> {log} 2>&1
        echo "Pruned variants kept: $(wc -l < {output.prune_in})" >> {log}
        """
