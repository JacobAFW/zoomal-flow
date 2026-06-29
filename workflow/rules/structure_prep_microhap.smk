# structure_prep_microhap.smk — microhap structure-prep seam (STUB)
# --------------------------------------------------------------------------
# Contract (DESIGN §4, seam 3):
#   When implemented, this file must emit the same seam outputs as the WGS
#   path:
#     {outputs}/structure/cleaned.{bed,bim,fam}     (filtered set)
#     {outputs}/structure/cleaned.ld.{bed,bim,fam}  (LD-pruned for ADMIXTURE)
#   ADMIXTURE / PLINK / PCA all assume biallelic SNPs. For microhap data
#   that means either:
#     (a) decompose multiallelic microhaps into biallelic SNPs (one record
#         per alternate allele), then run the WGS prep verbatim from
#         normalise_vcf onward; OR
#     (b) compute a pairwise distance matrix directly from multiallelic
#         genotypes and replace the bfile-driven downstream with PCA on the
#         distance matrix (NJ on the same matrix is unchanged).
#   The seam output names are fixed so 03_structure.smk's ADMIXTURE / PCA /
#   distance / assign_clusters rules don't fork on input_type.
# --------------------------------------------------------------------------


rule microhap_structure_prep_stub:
    """
    STUB. Emits the seam paths so the downstream DAG wires up identically;
    actually running it fails fast with the contract message.

    WHAT: echo MicrohapNotImplemented >&2 ; exit 1
    WHY:  Refuses to let a half-built microhap path produce a biallelic-
          shaped bfile from a non-biallelic dataset.
    TUNABLES: (none until implemented)
    OUTPUT: {outputs}/structure/cleaned.{bed,bim,fam} + cleaned.ld.{bed,bim,fam}
            (never actually written)
    TRY:    swap cohort.input_type to "wgs" — this rule disappears and the
            real WGS prep wires in instead.
    """
    input:
        vcf = f"{PATHS['outputs']}/qc/snps.qc.vcf.gz",
    output:
        bed    = f"{PATHS['outputs']}/structure/cleaned.bed",
        bim    = f"{PATHS['outputs']}/structure/cleaned.bim",
        fam    = f"{PATHS['outputs']}/structure/cleaned.fam",
        ld_bed = f"{PATHS['outputs']}/structure/cleaned.ld.bed",
        ld_bim = f"{PATHS['outputs']}/structure/cleaned.ld.bim",
        ld_fam = f"{PATHS['outputs']}/structure/cleaned.ld.fam",
    log:
        f"{PATHS['logs']}/structure/microhap_structure_prep_stub.log",
    message:
        "[structure:microhap] STUB — see structure_prep_microhap.smk for contract"
    shell:
        r"""
        echo "MicrohapNotImplemented: must emit a biallelic-decomposed PLINK set OR a distance matrix the structure analysis can consume; ADMIXTURE/PLINK assume biallelic SNPs — see DESIGN §4 seam 3." >&2
        exit 1
        """
