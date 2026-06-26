# moi_microhap.smk — microhap MOI/Fws seam (STUB)
# --------------------------------------------------------------------------
# Contract (DESIGN §4, seam 2):
#   When implemented, this file must emit
#     {outputs}/moi/fws_MOI.tsv
#   with the seam's two columns (sample, Proportion). For microhap data
#   moimix's biallelic-BAF model does NOT apply; the path must call a
#   multiallelic-aware MOI estimator (e.g. THE REAL McCOIL, or a
#   heterozygosity-based score over locus allele counts). The cutoff
#   semantics also differ from the WGS Fws convention — a microhap-specific
#   `moi.fws_polyclonal_cutoff` may be required, but the SHAPE of the
#   downstream file (sample + score column) is fixed so the seam is
#   fork-free for Stage 3+.
#
# Stub: exit 1 with the contract message so an accidental microhap run
# never silently produces wrong Fws.
# --------------------------------------------------------------------------


rule microhap_moi_stub:
    """
    STUB. Emits the seam path so downstream rules wire up identically;
    actually running it fails fast with the contract message.

    WHAT: echo MicrohapNotImplemented >&2 ; exit 1
    WHY:  Refuses to let a half-built microhap path silently produce a
          biallelic-shaped Fws file from a non-biallelic dataset.
    TUNABLES: (none until implemented)
    OUTPUT: {outputs}/moi/fws_MOI.tsv  (never actually written)
    TRY:    switch cohort.input_type to "wgs" — this rule disappears, the
            real moimix chain wires in instead.
    """
    input:
        vcf = f"{PATHS['outputs']}/qc/snps.qc.vcf.gz",
    output:
        fws = f"{PATHS['outputs']}/moi/fws_MOI.tsv",
    log:
        f"{PATHS['logs']}/moi/microhap_moi_stub.log",
    message:
        "[moi:microhap] STUB — see moi_microhap.smk for contract"
    shell:
        r"""
        echo "MicrohapNotImplemented: must emit fws_MOI.tsv (sample, complexity score) from a multiallelic-aware MOI estimator — moimix's biallelic-BAF model does not apply. The Fws cutoff semantics differ for microhap; see DESIGN §4 seam 2." >&2
        exit 1
        """
