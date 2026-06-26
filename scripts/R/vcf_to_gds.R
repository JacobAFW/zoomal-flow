#!/usr/bin/env Rscript
# vcf_to_gds.R  (agnostic port — generic, identical logic to V1)
# --------------------------------------------------------------------------
# Convert a VCF to SeqArray GDS format. GDS is the on-disk representation
# moimix consumes — random-access by variant and sample, much faster than
# tabix VCF for the per-sample operations in moimix.
#
# Usage:
#   Rscript vcf_to_gds.R <in.vcf.gz> <out.gds>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(SeqArray)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript vcf_to_gds.R <in.vcf.gz> <out.gds>")
vcf_in  <- args[1]
gds_out <- args[2]

stopifnot(file.exists(vcf_in))
dir.create(dirname(gds_out), showWarnings = FALSE, recursive = TRUE)

message(sprintf("[vcf_to_gds] converting %s -> %s", vcf_in, gds_out))

SeqArray::seqVCF2GDS(
  vcf.fn         = vcf_in,
  out.fn         = gds_out,
  storage.option = "LZ4_RA",
  parallel       = TRUE
)

gds <- SeqArray::seqOpen(gds_out)
on.exit(SeqArray::seqClose(gds))
SeqArray::seqSummary(gds)

message(sprintf("[vcf_to_gds] done: %s", gds_out))
