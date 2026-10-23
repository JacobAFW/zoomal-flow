#!/usr/bin/env Rscript
# run_selection_scan.R  (agnostic — role-driven, no cohort hardcodes)
# --------------------------------------------------------------------------
# Stage 6. Cross-population selection scan (XP-EHH + Rsb) for a two-group
# contrast, plus a within-population iHS scan per group.
#
# This REPLACES scripts/R/run_rehh_ihs.R, which had four defects that together
# made its output meaningless. All four are fixed here; the detail and the
# measurements behind each are in CONTEXT-pipeline.md.
#
#   1. iHS cannot answer a two-group question. It is a WITHIN-population
#      statistic, and the old script pooled case and control into one keep-list
#      and scanned the pool, so the contrast never entered the statistic. The
#      cross-population statistics that exist for this are XP-EHH
#      (`ies2xpehh`) and Rsb (`ines2rsb`), and they take two SEPARATE scans.
#      Within-group iHS is still produced, as its own per-group output.
#
#   2. The input was the Stage-3 `cleaned` bfile — a 24,641-variant set
#      filtered cohort-wide for ADMIXTURE and PCA. At roughly one marker per
#      6.5 kb, EHH cannot decay before the integration runs out of chromosome,
#      so `scan_hh` returns NA almost everywhere. This reads the post-QC SNP
#      set instead, and assesses missingness WITHIN each group rather than
#      inheriting a cohort-wide figure.
#
#   3. `data2haplohh` reads a diploid-coded VCF and emits two haplotypes per
#      sample. For a haploid organism with monoclonal isolates that counts
#      every genome twice as a pair of identical haplotypes, which guarantees
#      haplotype homozygosity and biases every EHH integral upward. Genotypes
#      are haploidised here before rehh sees them, so nhap equals the genome
#      count.
#
#   4. `min_maf` was hardcoded at 0.05 and applied twice. Every filter is now
#      a parameter.
#
# AND the reporting contract that stops the stage being actively misleading:
# rehh standardizes these statistics to mean 0 / SD 1 across the markers in
# the scan and reads p-values off a normal. When two groups are similar that
# standardisation denominator collapses and trivial differences become extreme
# z-scores, so a weakly-differentiated contrast can return hundreds of
# "candidate regions" that are artefacts of a vanishing denominator rather
# than signal. So **no candidate count is emitted without its tail diagnostic**
# (kurtosis, max |z|, the expected maximum under normality, and the tail
# excess), and the summary carries a verdict string derived from them.
#
# Filters reference the canonical role columns (`group`, `geography`,
# `country`, `sample_id`) joined with the Stage-3 cluster assignment, exactly
# as the old script did — no cohort-specific names.
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(rehh)
  library(data.table)
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
required <- c("vcf", "metadata", "clusters", "model", "case-filter",
              "control-filter", "out-dir")
missing <- setdiff(required, names(args))
if (length(missing) > 0) stop("Missing args: ", paste(missing, collapse = ", "))

num <- function(key, default) {
  if (is.null(args[[key]])) default else as.numeric(args[[key]])
}
model        <- args[["model"]]
case_filter  <- args[["case-filter"]]
ctrl_filter  <- args[["control-filter"]]
min_callrate <- num("min-callrate", 0.9)
min_maf      <- num("min-maf", 0)
max_het_rate <- num("max-het-rate", 0.1)
ihs_min_maf  <- num("ihs-min-maf", 0.05)
threshold    <- num("threshold", 4)
window_size  <- num("window-size", 10000)
overlap      <- num("overlap", 1000)
min_extr_mrk <- num("min-extr-mrk", 3)
min_group_n  <- num("min-group-n", 5)

out_dir <- file.path(args[["out-dir"]], model)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bcftools <- Sys.which("bcftools")
if (!nzchar(bcftools)) stop("bcftools not on PATH — needed for VCF access")

# -- Merge canonical-role metadata with the cluster assignment -------------
meta <- read_tsv(args[["metadata"]], show_col_types = FALSE)
stopifnot("sample_id" %in% names(meta))
if ("group" %in% names(meta)) meta <- meta %>% dplyr::select(-group)

clu <- read_tsv(args[["clusters"]], show_col_types = FALSE) %>%
  dplyr::select(sample_id = Sample, group = Cluster) %>% distinct()
m <- meta %>% left_join(clu, by = "sample_id")

.drop_arg <- args[["exclude"]]
if (!is.null(.drop_arg) && nzchar(.drop_arg) &&
    !(.drop_arg %in% c("NULL", "None")) && file.exists(.drop_arg)) {
  d <- readLines(.drop_arg, warn = FALSE); d <- d[nzchar(d)]
  if (length(d) > 0) {
    n_before <- nrow(m)
    m <- m %>% filter(!(sample_id %in% d))
    message(sprintf("[selection] de-clonalized: %d of %d samples eligible",
                    nrow(m), n_before))
  }
}

case_ids <- m %>% filter(!!rlang::parse_expr(case_filter)) %>% pull(sample_id)
ctrl_ids <- m %>% filter(!!rlang::parse_expr(ctrl_filter)) %>% pull(sample_id)

.floor_msg <- function(side, n) sprintf(
  "%s set has %d samples (< %d) — cannot scan. If this is the de-clonalized arm, the full arm may well have cleared the minimum, which itself says the set was largely clonal replicates.",
  side, n, min_group_n)
if (length(case_ids) < min_group_n) stop(.floor_msg("Case", length(case_ids)))
if (length(ctrl_ids) < min_group_n) stop(.floor_msg("Control", length(ctrl_ids)))

groups <- list(case = case_ids, control = ctrl_ids)
message(sprintf("[selection] model=%s | case=%d | control=%d",
                model, length(case_ids), length(ctrl_ids)))

# PLINK 2's --recode vcf emits doubled sample IDs; accept either form.
double_id <- function(x) paste0(x, "_", x)
vcf_samples <- system2(bcftools, c("query", "-l", args[["vcf"]]), stdout = TRUE)
resolve <- function(ids) {
  hit <- ifelse(ids %in% vcf_samples, ids,
                ifelse(double_id(ids) %in% vcf_samples, double_id(ids), NA))
  if (anyNA(hit)) {
    stop(sprintf("%d sample(s) not present in the VCF, e.g. %s",
                 sum(is.na(hit)), paste(head(ids[is.na(hit)], 3), collapse = ", ")))
  }
  hit
}
groups <- lapply(groups, resolve)

keep_path <- file.path(out_dir, "scan_keep.tsv")
writeLines(unlist(groups, use.names = FALSE), keep_path)

# Contigs come from the VCF header — generic, no cohort assumption. MIT/API
# style contigs were excluded upstream by reference.exclude_contigs.
hdr <- system2(bcftools, c("view", "-h", args[["vcf"]]), stdout = TRUE)
chroms <- hdr %>% grep("^##contig=<ID=", ., value = TRUE) %>%
  sub("^.*ID=", "", .) %>% sub(",.*$", "", .) %>% sub(">$", "", .)

# -- Read one chromosome's genotypes for the two groups --------------------
# bcftools query streams POS + GT; data.table::fread parses it. Haploidising
# here is defect 3's fix: 0/0 -> 0, 1/1 -> 1, and a heterozygous call becomes
# missing, because a het in a haploid monoclonal isolate has no single
# haplotype and rounding it either way would invent one.
read_chrom <- function(chrom) {
  cmd <- sprintf("%s view -r %s -S %s --force-samples %s 2>/dev/null | %s query -f '%%POS[\\t%%GT]\\n'",
                 shQuote(bcftools), shQuote(chrom), shQuote(keep_path),
                 shQuote(args[["vcf"]]), shQuote(bcftools))
  # Contigs with no records (e.g. organellar contigs excluded upstream, whose
  # ##contig headers survive) make fread warn about a zero-size file. That is
  # expected, not a problem, so it is suppressed rather than printed 14 times.
  dt <- tryCatch(suppressWarnings(
                   fread(cmd = cmd, header = FALSE, sep = "\t",
                         colClasses = "character", showProgress = FALSE)),
                 error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  pos <- as.numeric(dt[[1]])
  gt <- as.matrix(dt[, -1, with = FALSE])          # markers x samples
  hap <- matrix(NA_integer_, nrow = nrow(gt), ncol = ncol(gt))
  hap[gt %in% c("0/0", "0|0", "0")] <- 0L
  hap[gt %in% c("1/1", "1|1", "1")] <- 1L
  het <- gt %in% c("0/1", "1/0", "0|1", "1|0")
  dim(het) <- dim(gt)
  list(pos = pos, hap = hap, het = het,
       samples = unlist(groups, use.names = FALSE))
}

# -- Scan ------------------------------------------------------------------
scans <- list(case = NULL, control = NULL)
manifest <- list()
idx <- list(case = seq_along(groups$case),
            control = length(groups$case) + seq_along(groups$control))

for (chrom in chroms) {
  blk <- read_chrom(chrom)
  if (is.null(blk)) next
  n_input <- length(blk$pos)

  keep <- rep(TRUE, n_input)
  for (g in names(idx)) {
    sub <- blk$hap[, idx[[g]], drop = FALSE]
    hets <- blk$het[, idx[[g]], drop = FALSE]
    called <- rowSums(!is.na(sub))
    het_rate <- ifelse(called > 0, rowSums(hets) / pmax(called, 1), 0)
    keep <- keep & (called >= min_callrate * length(idx[[g]])) &
      (het_rate <= max_het_rate)
  }
  # Polymorphism is assessed on the POOLED case+control sample, deliberately.
  # A site fixed in one group and variable in the other is the most
  # informative kind of site for a cross-population statistic, and a per-group
  # MAF floor would delete exactly those. Missingness stays per-group above.
  pooled <- blk$hap
  called_p <- rowSums(!is.na(pooled))
  af <- ifelse(called_p > 0, rowSums(pooled, na.rm = TRUE) / pmax(called_p, 1), NA)
  maf <- pmin(af, 1 - af)
  keep <- keep & !is.na(maf) & maf > 0 & (min_maf <= 0 | maf >= min_maf)

  manifest[[chrom]] <- data.frame(
    chrom = chrom, n_input = n_input, n_markers = sum(keep))
  if (sum(keep) < 2) next

  pos_k <- blk$pos[keep]
  for (g in names(idx)) {
    h <- t(blk$hap[keep, idx[[g]], drop = FALSE])   # haplotypes x markers
    storage.mode(h) <- "integer"
    colnames(h) <- sprintf("%s:%d", chrom, as.integer(pos_k))
    rownames(h) <- sprintf("%s_%d", g, seq_len(nrow(h)))
    hh <- new("haplohh", chr.name = as.character(match(chrom, chroms)),
              positions = as.numeric(pos_k), haplo = h)
    # polarized = FALSE: there is no outgroup, so ancestral state is unknown
    # and rehh works on major/minor. The old script did the same.
    s <- scan_hh(hh, polarized = FALSE, phased = TRUE)
    scans[[g]] <- if (is.null(scans[[g]])) s else rbind(scans[[g]], s)
  }
  message(sprintf("[selection]   %s: %d of %d markers scanned",
                  chrom, sum(keep), n_input))
}

manifest <- bind_rows(manifest)
write_tsv(manifest, file.path(out_dir, "marker_manifest.tsv"))

empty_out <- function(reason) {
  write_tsv(tibble(CONTIG = character(), CHR = character(), START = numeric(),
                   END = numeric(), N_MRK = integer()),
            file.path(out_dir, "candidate_regions_XPEHH.tsv"))
  write_tsv(tibble(statistic = character(), markers = integer(),
                   note = reason), file.path(out_dir, "tail_diagnostics.tsv"))
  writeLines(reason, file.path(out_dir, "scan_summary.tsv"))
  quit(status = 0)
}
if (is.null(scans$case) || nrow(scans$case) == 0) {
  empty_out(sprintf("no markers survived the filters (call rate >= %g, MAF >= %g) — nothing was tested",
                    min_callrate, min_maf))
}

code_to_name <- setNames(chroms, as.character(seq_along(chroms)))
add_contig <- function(df) {
  if (nrow(df) == 0) return(df)
  cbind(CONTIG = code_to_name[as.character(df$CHR)], df)
}

# -- The tail diagnostic ---------------------------------------------------
# rehh reads p-values off a standard normal after standardising. These four
# numbers say whether that assumption holds for THIS scan, and no candidate
# count is reported without them.
expected_max <- function(n) {
  if (n < 2) return(NA_real_)
  l <- log(n)
  sqrt(2 * l) - (log(l) + log(4 * pi)) / (2 * sqrt(2 * l))
}
tail_diag <- function(v, name) {
  v <- v[is.finite(v)]
  n <- length(v)
  if (n < 100) {
    return(tibble(statistic = name, markers = n, kurtosis = NA_real_,
                  max_abs_z = NA_real_, expected_max = NA_real_,
                  obs_gt5 = NA_integer_, exp_gt5 = NA_real_,
                  excess_5 = NA_real_,
                  verdict = "UNDERPOWERED — too few markers to assess"))
  }
  kurt <- mean((v - mean(v))^4) / stats::sd(v)^4
  mx <- max(abs(v)); emax <- expected_max(n)
  obs <- sum(abs(v) > 5); expd <- n * 2 * stats::pnorm(5, lower.tail = FALSE)
  exc <- if (expd > 0) obs / expd else NA_real_
  verdict <- if (kurt > 4) {
    sprintf("DEGENERATE — kurtosis %.1f; the standardisation denominator is too small, so neither the counts nor the p-values mean anything", kurt)
  } else if (mx <= emax + 0.75) {
    sprintf("FLAT — max |z| %.2f against an expected maximum of %.2f; indistinguishable from noise", mx, emax)
  } else if (!is.na(exc) && exc >= 10) {
    sprintf("EXCESS — well-behaved bulk (kurtosis %.2f) with %.0fx the normal expectation beyond |z|=5; NOT evidence of selection until compared with a matched-design null", kurt, exc)
  } else {
    sprintf("MARGINAL — %.0fx the normal expectation beyond |z|=5, kurtosis %.2f", exc, kurt)
  }
  tibble(statistic = name, markers = n, kurtosis = kurt, max_abs_z = mx,
         expected_max = emax, obs_gt5 = obs, exp_gt5 = expd, excess_5 = exc,
         verdict = verdict)
}

cand_for <- function(tbl, stat_col, name) {
  df <- as.data.frame(tbl)
  df <- df[is.finite(df[[stat_col]]), , drop = FALSE]
  if (nrow(df) == 0) return(list(n = 0L, top = NA_character_))
  cand <- tryCatch(
    calc_candidate_regions(df, threshold = threshold, pval = TRUE,
                           window_size = window_size, overlap = overlap,
                           min_n_extr_mrk = min_extr_mrk, ignore_sign = TRUE),
    error = function(e) data.frame())
  out <- if (is.null(cand) || nrow(cand) == 0) {
    tibble(CONTIG = character(), CHR = character(), START = numeric(),
           END = numeric(), N_MRK = integer())
  } else add_contig(cand)
  write_tsv(as_tibble(out),
            file.path(out_dir, sprintf("candidate_regions_%s.tsv", name)))
  best <- df[which.max(df$LOGPVALUE), , drop = FALSE]
  list(n = nrow(out),
       top = sprintf("%s:%d", code_to_name[as.character(best$CHR)],
                     as.integer(best$POSITION)))
}

diags <- list(); summary_rows <- list()

# -- Cross-population statistics: the actual two-group test ----------------
xp <- ies2xpehh(scans$case, scans$control, popname1 = "case",
                popname2 = "control", p.side = NA, verbose = FALSE)
rs <- ines2rsb(scans$case, scans$control, popname1 = "case",
               popname2 = "control", p.side = NA, verbose = FALSE)
write_tsv(as_tibble(add_contig(as.data.frame(xp))),
          file.path(out_dir, "xpehh_table.tsv"))
write_tsv(as_tibble(add_contig(as.data.frame(rs))),
          file.path(out_dir, "rsb_table.tsv"))

for (p in list(list(t = xp, pre = "XPEHH"), list(t = rs, pre = "RSB"))) {
  col <- grep(paste0("^", p$pre), colnames(p$t), value = TRUE)[1]
  d <- tail_diag(as.data.frame(p$t)[[col]], p$pre)
  cc <- cand_for(p$t, col, p$pre)
  diags[[p$pre]] <- d
  summary_rows[[p$pre]] <- tibble(
    statistic = p$pre, markers = d$markers, candidate_regions = cc$n,
    top_marker = cc$top, kurtosis = d$kurtosis, max_abs_z = d$max_abs_z,
    expected_max = d$expected_max, excess_beyond_z5 = d$excess_5,
    verdict = d$verdict)
}

# -- Within-population iHS, one per group (legitimate, but not a test) -----
for (g in names(scans)) {
  ihs <- ihh2ihs(scans[[g]], min_maf = ihs_min_maf, freqbin = 0.025,
                 p.side = NA, verbose = FALSE)
  write_tsv(as_tibble(add_contig(as.data.frame(ihs$ihs))),
            file.path(out_dir, sprintf("ihs_%s_table.tsv", g)))
  nm <- sprintf("iHS_%s", g)
  d <- tail_diag(as.data.frame(ihs$ihs)[["IHS"]], nm)
  cc <- cand_for(ihs$ihs, "IHS", nm)
  diags[[nm]] <- d
  summary_rows[[nm]] <- tibble(
    statistic = nm, markers = d$markers, candidate_regions = cc$n,
    top_marker = cc$top, kurtosis = d$kurtosis, max_abs_z = d$max_abs_z,
    expected_max = d$expected_max, excess_beyond_z5 = d$excess_5,
    verdict = d$verdict)
}

write_tsv(bind_rows(diags), file.path(out_dir, "tail_diagnostics.tsv"))

summary_tbl <- bind_rows(summary_rows) %>%
  mutate(model = model, n_case = length(groups$case),
         n_control = length(groups$control),
         markers_scanned = nrow(scans$case),
         source_variants = sum(manifest$n_input),
         min_callrate = min_callrate, min_maf = min_maf,
         ihs_min_maf = ihs_min_maf, threshold = threshold,
         window_size_bp = window_size) %>%
  dplyr::select(model, n_case, n_control, source_variants, markers_scanned,
                everything())
write_tsv(summary_tbl, file.path(out_dir, "scan_summary.tsv"))

cat(sprintf("\n[selection] %s | case %d / control %d | %s of %s source variants scanned\n",
            model, length(groups$case), length(groups$control),
            format(nrow(scans$case), big.mark = ","),
            format(sum(manifest$n_input), big.mark = ",")))
for (i in seq_len(nrow(summary_tbl))) {
  r <- summary_tbl[i, ]
  cat(sprintf("  %-12s markers %-9s regions %-5s | %s\n",
              r$statistic, format(r$markers, big.mark = ","),
              r$candidate_regions, r$verdict))
}
