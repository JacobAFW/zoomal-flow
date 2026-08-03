#!/usr/bin/env Rscript
# run_rehh_ihs.R  (agnostic port — role-driven, no cohort hardcodes)
# --------------------------------------------------------------------------
# Genome-wide iHS scan via rehh. The case/control filters are R expressions
# evaluated against the merged canonical-role metadata + cluster assignment
# — so `case_filter` and `control_filter` reference the role columns
# (`group`, `geography`, `country`, `sample_id`), not the cohort-specific
# `Cluster`/`State` names V1 used.
#
# Example (Aceh vs other-Peninsular):
#   --case-filter    "group == 'Peninsular' & geography == 'Aceh'"
#   --control-filter "group == 'Peninsular' & geography != 'Aceh'"
#
# Sample-ID matching: PLINK 2's --recode vcf emits `<id>_<id>` doubled IDs.
# We build the keep list from clean IDs, and also try the doubled form so
# both work.
#
# Empty-result handling: a negative result is a finding; write header-only
# candidate_regions_iHS.tsv + a one-line ihs_summary.tsv rather than
# crashing.
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(rehh)
  library(data.table)
  library(R.utils)
})

parse_args <- function(av) {
  out <- list(); i <- 1
  while (i <= length(av)) {
    key <- av[i]
    if (!startsWith(key, "--")) stop("Bad arg: ", key)
    out[[sub("^--", "", key)]] <- av[i + 1]; i <- i + 2
  }
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("vcf","metadata","clusters","model","case-filter","control-filter",
              "threshold","p-threshold","window-size","overlap","min-extr-mrk","out-dir")
missing  <- setdiff(required, names(args))
if (length(missing) > 0) stop("Missing args: ", paste(missing, collapse = ", "))

model         <- args[["model"]]
case_filter   <- args[["case-filter"]]
ctrl_filter   <- args[["control-filter"]]
threshold     <- as.numeric(args[["threshold"]])
p_threshold   <- as.numeric(args[["p-threshold"]])
window_size   <- as.integer(args[["window-size"]])
overlap       <- as.integer(args[["overlap"]])
min_extr_mrk  <- as.integer(args[["min-extr-mrk"]])

out_dir <- file.path(args[["out-dir"]], model)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -- Merge canonical-role metadata with cluster assignment --------------
# admix_clusters.tsv's `Cluster` is authoritative (Stage 3 label). The
# metadata's `group` role (if present) is the user's raw column; drop it
# so the merged frame has one unambiguous `group` = Cluster.
meta <- read_tsv(args[["metadata"]], show_col_types = FALSE)
stopifnot("sample_id" %in% names(meta))
if ("group" %in% names(meta)) {
  meta <- meta %>% dplyr::select(-group)
}

clu <- read_tsv(args[["clusters"]], show_col_types = FALSE) %>%
  dplyr::select(sample_id = Sample, group = Cluster) %>% distinct()

m <- meta %>% left_join(clu, by = "sample_id")

case_ids  <- m %>% filter(!!rlang::parse_expr(case_filter))  %>% pull(sample_id)
ctrl_ids  <- m %>% filter(!!rlang::parse_expr(ctrl_filter))  %>% pull(sample_id)
if (length(case_ids)  < 5) stop("Case set < 5 samples — cannot run iHS.")
if (length(ctrl_ids)  < 5) stop("Control set < 5 samples — cannot run iHS.")
keep_ids <- unique(c(case_ids, ctrl_ids))
message(sprintf("[run_rehh_ihs] model=%s | case=%d | control=%d | total=%d",
                model, length(case_ids), length(ctrl_ids), length(keep_ids)))

# The VCF may carry PLINK's doubled sample IDs (`X_X`); include both forms
# in the keep set so --force-samples matches.
double_id <- function(x) paste0(x, "_", x)
keep_set <- unique(c(keep_ids, double_id(keep_ids)))

# -- Subset VCF ---------------------------------------------------------
keep_path <- file.path(out_dir, "rehh_keep.tsv")
writeLines(keep_set, keep_path)

bcftools <- Sys.which("bcftools")
if (!nzchar(bcftools)) stop("bcftools not on PATH — needed for VCF subsetting")

vcf_in  <- args[["vcf"]]
sub_vcf <- file.path(out_dir, "rehh_input.vcf.gz")
system2(bcftools, c("view","-S", keep_path, "--force-samples",
                    "-Oz","-o", sub_vcf, vcf_in), wait = TRUE)
system2(bcftools, c("index","-c", sub_vcf), wait = TRUE)

# Nuclear-contig detection is generic: read the ##contig lines from the
# header. The Stage-4 combined VCF has already been through Stage 3's
# plink2 --update-chr, so contig IDs are integers already.
hdr <- system2(bcftools, c("view","-h", sub_vcf), stdout = TRUE)
all_contigs <- hdr %>%
  grep("^##contig=<ID=", ., value = TRUE) %>%
  sub("^.*ID=", "", .) %>%
  sub(",.*$", "", .) %>%
  sub(">$", "", .)
# MIT/API-style haploid contigs were excluded upstream (reference.exclude_contigs
# in Stage 0), so all_contigs is already nuclear.
chroms <- all_contigs
message(sprintf("[run_rehh_ihs] %d contigs detected: %s",
                length(chroms), paste(chroms, collapse = ", ")))

# -- Per-chromosome scan_hh + accumulate wgscan -------------------------
wgscan <- NULL
for (chrom in chroms) {
  per_chrom_vcf <- file.path(out_dir, sprintf("rehh_%s.vcf.gz", chrom))
  system2(bcftools, c("view","-r", chrom, "-Oz","-o", per_chrom_vcf, sub_vcf),
          wait = TRUE)
  hh <- tryCatch({
    data2haplohh(hap_file = per_chrom_vcf,
                 chr.name = chrom,
                 polarize_vcf = FALSE,
                 min_perc_geno.mrk = 100,
                 min_maf = 0.05,
                 vcf_reader = "data.table")
  }, error = function(e) {
    message(sprintf("[run_rehh_ihs] %s: data2haplohh skipped (%s)",
                    chrom, conditionMessage(e)))
    NULL
  })
  if (is.null(hh)) next
  scan <- scan_hh(hh)
  wgscan <- if (is.null(wgscan)) scan else rbind(wgscan, scan)
}

empty_df <- function(cols) {
  as.data.frame(setNames(replicate(length(cols), character(0), simplify = FALSE), cols),
                stringsAsFactors = FALSE)
}

if (is.null(wgscan) || nrow(wgscan) == 0) {
  warning("[run_rehh_ihs] wgscan empty — writing empty outputs")
  cand_cols <- c("CHR","START","END","N_MRK","MEAN_MRK","MAX_MRK",
                 "N_EXTR_MRK","PERC_EXTR_MRK","Stat")
  ihs_cols  <- c("CHR","POSITION","IHS","LOGPVALUE")
  write_tsv(empty_df(cand_cols), file.path(out_dir, "candidate_regions_iHS.tsv"))
  write_tsv(empty_df(ihs_cols),  file.path(out_dir, "ihs_table.tsv"))
  writeLines(sprintf("no candidate regions at p < %g (wgscan empty)", p_threshold),
             file.path(out_dir, "ihs_summary.tsv"))
  saveRDS(NULL, file.path(out_dir, "wgs_ihs.rds"))
  quit(status = 0)
}

wgs_ihs <- ihh2ihs(na.omit(wgscan), min_maf = 0.05, freqbin = 0.025)
saveRDS(wgs_ihs, file.path(out_dir, "wgs_ihs.rds"))

# Chromosome ID normalisation: rehh returns whatever it was fed. If the
# CHR is not already numeric, best-effort coercion (strip common prefixes/
# suffixes then to numeric). Any that stay NA are left as character.
wgs_ihs2 <- wgs_ihs
wgs_ihs2$ihs <- wgs_ihs2$ihs %>%
  mutate(CHR_num = suppressWarnings(as.numeric(as.character(CHR)))) %>%
  mutate(CHR = ifelse(is.na(CHR_num), as.character(CHR), CHR_num)) %>%
  dplyr::select(-CHR_num) %>%
  as.data.frame()

cand <- tryCatch({
  calc_candidate_regions(wgs_ihs2,
                         threshold = threshold,
                         pval = TRUE,
                         window_size = window_size,
                         overlap = overlap,
                         min_n_extr_mrk = min_extr_mrk)
}, error = function(e) {
  message(sprintf("[run_rehh_ihs] calc_candidate_regions error: %s",
                  conditionMessage(e)))
  data.frame()
})

if (is.null(cand) || nrow(cand) == 0) {
  cand_cols <- c("CHR","START","END","N_MRK","MEAN_MRK","MAX_MRK",
                 "N_EXTR_MRK","PERC_EXTR_MRK","Stat")
  write_tsv(empty_df(cand_cols), file.path(out_dir, "candidate_regions_iHS.tsv"))
  writeLines(sprintf("no candidate regions at p < %g", p_threshold),
             file.path(out_dir, "ihs_summary.tsv"))
  cat(sprintf(
    "\n[run_rehh_ihs] %s: no candidate regions at p < %g — negative result recorded\n",
    model, p_threshold))
} else {
  cand <- as_tibble(cand) %>% add_column(Stat = "iHS")
  write_tsv(cand, file.path(out_dir, "candidate_regions_iHS.tsv"))
  cand %>% arrange(desc(N_EXTR_MRK)) %>% head(20) %>%
    write_tsv(file.path(out_dir, "ihs_summary.tsv"))
  cat(sprintf("\n[run_rehh_ihs] %s: %d candidate region(s) at threshold=%g, p<%g\n",
              model, nrow(cand), threshold, p_threshold))
}

wgs_ihs$ihs %>% as_tibble() %>%
  write_tsv(file.path(out_dir, "ihs_table.tsv"))
