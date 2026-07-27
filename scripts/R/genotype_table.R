#!/usr/bin/env Rscript
# genotype_table.R  (agnostic port — contig_map replaces the ordered_PKNH gsub)
# --------------------------------------------------------------------------
# Convert a VCF to hmmIBD's genotype-table format:
#
#   chrom (integer)  pos (integer)  sample1  sample2  ...
#
# where each sample column is 0/1/2/…/-1 for haploid calls (-1 = missing).
# For P. knowlesi and other Plasmodium, hmmIBD treats homozygous diploid
# genotypes as the haploid allele; heterozygous or missing map to -1.
#
# Contig → integer mapping is read from outputs/setup/contig_map.tsv
# (built by Stage 0's build_contig_map rule from the reference .fai) —
# replaces V1's hardcoded `sub("_v2","",sub("ordered_PKNH_","",chrom))`.
#
# Usage:
#   Rscript genotype_table.R <in.vcf.gz> <contig_map.tsv> <out.tsv>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript genotype_table.R <in.vcf.gz> <contig_map.tsv> <out.tsv>")
}
vcf_in    <- args[1]
cmap_in   <- args[2]
out_tsv   <- args[3]

sample_ids_raw <- system2("bcftools", c("query", "-l", vcf_in), stdout = TRUE)

# PLINK 2 --recode vcf writes sample IDs as FID_IID; when FID == IID (the
# --double-id default), the header carries doubled names like SAMPLE_SAMPLE.
# Undouble here so downstream (hmmIBD, clonal_clusters) sees clean IDs.
.undouble <- function(x) {
  vapply(x, function(s) {
    n <- nchar(s)
    if (n %% 2 == 1) {
      mid <- (n + 1) / 2
      if (substr(s, mid, mid) == "_" &&
          substr(s, 1, mid - 1) == substr(s, mid + 1, n)) {
        return(substr(s, 1, mid - 1))
      }
    }
    s
  }, character(1), USE.NAMES = FALSE)
}
sample_ids <- .undouble(sample_ids_raw)
message(sprintf("[genotype_table] samples: %d (%d undoubled)",
                length(sample_ids), sum(sample_ids != sample_ids_raw)))

cmap <- fread(cmap_in, header = FALSE, sep = "\t",
              col.names = c("chrom_str", "chrom"))

TAB <- "\t"
NL  <- "\n"
gt_format <- paste0("%CHROM", TAB, "%POS", "[", TAB, "%GT]", NL)
tmp_tsv <- tempfile(fileext = ".tsv")
status <- system2("bcftools",
                  c("query", "-f", shQuote(gt_format), vcf_in),
                  stdout = tmp_tsv)
if (status != 0) stop("bcftools query failed (exit ", status, ")")

gt <- fread(tmp_tsv, header = FALSE, sep = "\t", data.table = TRUE,
            colClasses = c("character", "integer",
                           rep("character", length(sample_ids))))
setnames(gt, c("chrom_str", "pos", sample_ids))
file.remove(tmp_tsv)

# Look up integer chrom. In the agnostic pipeline the VCF handed to this
# script has already gone through `plink2 --update-chr` (Stage 3), so its
# CHROM column is already the integer code. If we ever run this against
# a raw VCF the string names get mapped via contig_map.tsv.
vcf_chroms <- unique(gt$chrom_str)
suppressWarnings(as_int <- as.integer(vcf_chroms))
if (all(!is.na(as_int))) {
  # Already integer-coded.
  gt[, chrom := as.integer(chrom_str)]
} else {
  missing_before <- setdiff(vcf_chroms, cmap$chrom_str)
  if (length(missing_before) > 0) {
    stop("Contigs in VCF not present in contig_map.tsv: ",
         paste(missing_before, collapse = ", "))
  }
  gt[, chrom := cmap$chrom[match(chrom_str, cmap$chrom_str)]]
}
gt[, chrom_str := NULL]

# hmmIBD encoding: 0/0→0, 1/1→1, 2/2→2, …, else→-1.
encode_gt <- function(x) {
  out <- rep("-1", length(x))
  for (a in 0:9) {
    hom <- paste0(a, "/", a)
    out[x == hom] <- as.character(a)
  }
  out[x == "0|0"] <- "0"
  out[x == "1|1"] <- "1"
  out
}
gt[, (sample_ids) := lapply(.SD, encode_gt), .SDcols = sample_ids]

# Reorder to (chrom, pos, samples…) and sort.
setcolorder(gt, c("chrom", "pos", sample_ids))
setorder(gt, chrom, pos)
fwrite(gt, out_tsv, sep = "\t", quote = FALSE, na = "-1")

cat("\n=== Summary ===\n")
cat(sprintf("Variants: %d\n", nrow(gt)))
cat(sprintf("Samples:  %d\n", length(sample_ids)))
cat(sprintf("Chroms:   %s\n",
            paste(sort(unique(gt$chrom)), collapse = ",")))
