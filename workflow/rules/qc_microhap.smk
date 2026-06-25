# qc_microhap.smk — microhap QC seam (STUB)
# --------------------------------------------------------------------------
# Contract (DESIGN §4):
#   When implemented, this file must emit
#     {outputs}/qc/snps.qc.vcf.gz   — a QC'd VCF with locus-level call-rate
#                                     QC (NOT bcftools -m2 -M2; microhap loci
#                                     are multiallelic by design)
#   plus a loci.tsv describing per-locus allele structure (locus_id, contig,
#   pos, n_alleles, allele_seqs, called_n, called_rate) so Stage 2's MOI seam
#   and Stage 3's structure-input prep can branch on real per-locus info
#   without re-parsing the VCF header for each rule.
#
# This file is intentionally a fast-fail stub: dry-running `input_type:
# microhap` shows the stub wired in; actually running it exits with a
# MicrohapNotImplemented message so no one accidentally treats the WGS
# downstream as valid for a microhap cohort.
# --------------------------------------------------------------------------


rule microhap_qc_stub:
    """
    STUB. Produces the same seam output as qc_wgs.smk's terminal rule,
    so Stage 2+ wiring is identical regardless of input_type. Running it
    fails fast with the contract message.

    WHAT: echo MicrohapNotImplemented >&2 ; exit 1
    WHY:  Forces explicit work to implement the microhap path; refuses to
          let a half-built path produce silently-wrong output.
    TUNABLES: (none until implemented)
    OUTPUT: {outputs}/qc/snps.qc.vcf.gz  (never actually written)
    TRY:    swap cohort.input_type to "wgs" to run the real WGS path; this
            rule disappears from the DAG.
    """
    input:
        vcf = QC_INPUT_VCF,
    output:
        vcf = f"{PATHS['outputs']}/qc/snps.qc.vcf.gz",
        idx = f"{PATHS['outputs']}/qc/snps.qc.vcf.gz.csi",
    log:
        f"{PATHS['logs']}/qc/microhap_qc_stub.log",
    message:
        "[qc:microhap] STUB — see qc_microhap.smk for contract"
    shell:
        r"""
        echo "MicrohapNotImplemented: this rule must emit a QC'd VCF (locus-level call-rate QC, multiallelic loci preserved — not -m2 -M2) and a loci.tsv describing per-locus allele structure. See agnostic/DESIGN.md §4. Microhap path not yet built." >&2
        exit 1
        """


# Downstream targets that the seam's WGS sibling provides
# (variant_count.txt + the four QC density figures) are not produced by the
# microhap stub yet. They'll be added when the real microhap QC is written:
# loci.tsv-derived equivalents of per-sample call-rate + per-locus allele
# count distributions belong here, NOT a port of the biallelic MAF script.
