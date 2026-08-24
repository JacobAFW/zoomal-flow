#!/usr/bin/env Rscript
# introgression_floor_derivation.R — derive the per-cluster support floor N
# by permutation, with an explicit false-window / FDR target (spec §9.5)
# --------------------------------------------------------------------------
# WHAT: builds a null distribution for per-window sample-support inside each
#       cluster, and reports the smallest support floor `N` at which chance
#       co-location essentially stops producing windows. That `N` becomes
#       `introgression.per_cluster_min_samples`.
# WHY:  Stage 5b established that the per-cluster filter — not the detection
#       rule — is what drives Stage 5's cluster-size behaviour, and that the
#       default should be an equal ABSOLUTE floor rather than a percentage of
#       cluster n. An absolute floor is only defensible if the number is
#       derived. Until now Stage 5 had no false-discovery control anywhere;
#       this is it.
#
# THE QUESTION. Filter 2 asks "is this window supported by enough samples of
# this cluster to be a shared event rather than a coincidence?". The failure
# mode is independent, sample-specific calls that happen to land on the same
# window. So: if every sample's calls were scattered at random over the
# windows where that sample could have been called, how much co-location would
# chance alone produce?
#
# THE NULL. Per cluster, and per sample within it:
#   - hold the sample's NUMBER of introgressed windows fixed (k_s), and
#   - hold the sample's ELIGIBLE window set fixed (E_s), then
#   - redraw which k_s of E_s it hits, uniformly without replacement.
# This destroys genomic co-location between samples while preserving per-sample
# call rate, per-sample testability, and cluster size — the three things that
# would otherwise confound the comparison between clusters.
#
# E_s IS COMPUTED, NOT ASSUMED. A sample can only be called in a window where
# the pipeline actually computed a distance for it, i.e. where it has more than
# `min_snps` non-missing calls. That set is per-sample (it follows marker
# density and per-sample missingness) and it is the denominator the whole
# argument rests on: assume "any 10 kb window in the genome" instead and the
# null is diluted, the tail shrinks, and N comes out too small. So this script
# recomputes the ELIGIBILITY table from the genotype table — window binning and
# a per-(sample, window) SNP count. That is much cheaper than detection: no
# consensus alleles, no distances, no density surface. Detection itself is NOT
# re-run; the observed calls are read from the cached per-pair files.
#
#   Caveat, stated because it moves N the anti-conservative way: the pipeline
#   also drops sample-windows whose two distances are exactly equal (V1's
#   `filter(DIST_X != DIST_Y)`), which needs the consensus. E_s here therefore
#   includes a small number of tied sample-windows the detector would have
#   dropped, making the universe marginally too large and the null marginally
#   too thin. The reported `eligible_tie_inflation` column quantifies it
#   against the observed call set.
#
# THE FILTER CHAIN IS APPLIED TO THE NULL TOO. Each permutation replicate goes
# through the dataset-level low-n filter (filter 1) exactly as the real data
# does before per-window support is measured. A null that skips the selection
# the data underwent would understate the tail.
#
# ONE GLOBAL N. The point of the floor is equal evidentiary support for every
# window, so N is a single number, not per-cluster. It is chosen as the
# smallest N meeting the target in EVERY cluster — which in practice means the
# largest cluster, where chance coincidence is easiest (more samples drawing
# from the same universe). The binding cluster is reported.
#
# CLI:
#   Rscript scripts/R/introgression_floor_derivation.R \
#     --genotype-table outputs/ibd/combined/hmmIBD_input.tsv \
#     --clusters       outputs/structure/admix_clusters.tsv \
#     --window-size    10000 \
#     --min-snps       5 \
#     --min-samples-per-window 2 \
#     --permutations   1000 \
#     --seed           20260824 \
#     --max-n          25 \
#     --target-false-windows 1 \
#     --fdr-target     0.05 \
#     --out            outputs/introgression/floor_derivation.tsv \
#     outputs/introgression/pairs/*.tsv
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

script_args <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- dirname(sub("^--file=", "",
                          script_args[grep("^--file=", script_args)][1]))
source(file.path(SCRIPT_DIR, "introgression_core.R"))

parse_args <- function(av) {
  flags <- list(); pos <- character(); i <- 1
  while (i <= length(av)) {
    if (startsWith(av[i], "--")) { flags[[sub("^--", "", av[i])]] <- av[i + 1]; i <- i + 2 }
    else { pos <- c(pos, av[i]); i <- i + 1 }
  }
  list(flags = flags, pos = pos)
}
parsed     <- parse_args(commandArgs(trailingOnly = TRUE))
args       <- parsed$flags
call_files <- parsed$pos

`%||%` <- function(a, b) if (is.null(a)) b else a
required <- c("genotype-table", "clusters", "out")
missing  <- setdiff(required, names(args))
if (length(missing) > 0) stop("Missing args: ", paste(missing, collapse = ", "))
if (length(call_files) == 0) stop("No per-pair call files given")

WINDOW    <- as.integer(args[["window-size"]] %||% 10000)
MIN_SNPS  <- as.integer(args[["min-snps"]] %||% 5)
MIN_SAMP  <- as.integer(args[["min-samples-per-window"]] %||% 2)
NPERM     <- as.integer(args[["permutations"]] %||% 1000)
SEED      <- as.integer(args[["seed"]] %||% 20260824)
MAX_N     <- as.integer(args[["max-n"]] %||% 25)
TARGET_FW <- as.numeric(args[["target-false-windows"]] %||% 1)
FDR_TGT   <- as.numeric(args[["fdr-target"]] %||% 0.05)
OUT       <- args[["out"]]
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

say <- function(...) message(sprintf("[floor_derivation] %s", sprintf(...)))
tick <- function() proc.time()[["elapsed"]]

# --------------------------------------------------------------------------
# Observed calls (cached — detection is NOT re-run)
# --------------------------------------------------------------------------
clusters_in <- read_tsv(args[["clusters"]], show_col_types = FALSE) %>%
  dplyr::select(SAMPLE = Sample, Cluster) %>% distinct()

calls <- call_files %>%
  map(~ read_tsv(.x, show_col_types = FALSE,
                 col_types = cols(WINDOW = col_character(), SAMPLE = col_character(),
                                  Cluster = col_character(), .default = col_guess()))) %>%
  bind_rows()
# One row per (sample, window): a sample appears in every pair its cluster is
# part of, and may hit the same window in more than one of them.
obs <- calls %>% distinct(Cluster, SAMPLE, WINDOW)
say("observed: %d raw calls -> %d distinct (sample, window) over %d windows, %d samples",
    nrow(calls), nrow(obs), n_distinct(obs$WINDOW), n_distinct(obs$SAMPLE))

# --------------------------------------------------------------------------
# Eligibility: which (sample, window) pairs the detector could have called.
# Window binning + a per-(sample, window) SNP count — the same `window_bin()`
# and the same strict `n > min_snps` cut the distance step uses.
# --------------------------------------------------------------------------
t0 <- tick()
gt <- fread(args[["genotype-table"]], header = TRUE, sep = "\t")
setnames(gt, 1:2, c("CHROM", "POS"))
gt[, CHROM := as.integer(CHROM)]
keep_samples <- intersect(clusters_in$SAMPLE, names(gt))
gt <- gt[, c("CHROM", "POS", keep_samples), with = FALSE]
gt[, BIN := window_bin(POS, WINDOW)]

elig <- melt(gt, id.vars = c("CHROM", "POS", "BIN"), variable.name = "SAMPLE",
             value.name = "SNP", variable.factor = FALSE)
elig <- elig[SNP >= 0]
elig <- elig[, .(n_snps = .N), by = .(SAMPLE, CHROM, BIN)][n_snps > MIN_SNPS]
elig[, WINDOW := window_id(CHROM, BIN)]
elig <- elig[, .(SAMPLE, WINDOW)]
rm(gt); invisible(gc(FALSE))
say("eligible (sample, window) pairs: %d over %d windows (%.0fs)",
    nrow(elig), uniqueN(elig$WINDOW), tick() - t0)

# Every observed call must be eligible — if not, the two steps disagree about
# window binning or the SNP cut and the null would be built on the wrong set.
elig_keys <- paste0(elig$SAMPLE, "\r", elig$WINDOW)
obs_keys  <- paste0(obs$SAMPLE,  "\r", obs$WINDOW)
orphans   <- setdiff(obs_keys, elig_keys)
if (length(orphans) > 0) {
  stop("[floor_derivation] ", length(orphans), " observed call(s) are not in the ",
       "eligible set — window binning or --min-snps disagrees with the cached calls")
}
say("all %d observed calls are inside the eligible set (consistency check passed)", nrow(obs))

# --------------------------------------------------------------------------
# Index everything to integers: permutation is index arithmetic from here.
# --------------------------------------------------------------------------
win_levels <- sort(unique(elig$WINDOW))
W <- length(win_levels)
elig[, w_idx := match(WINDOW, win_levels)]
setorder(elig, SAMPLE, w_idx)

samples <- clusters_in %>%
  filter(SAMPLE %in% unique(obs$SAMPLE)) %>%
  arrange(SAMPLE)
k_by_sample <- obs %>% count(SAMPLE, name = "k")
samples <- samples %>% left_join(k_by_sample, by = "SAMPLE") %>% filter(k > 0)

# Per-sample eligible window indices, in the same order as `samples`.
elig_split <- split(elig$w_idx, elig$SAMPLE)
E_list <- elig_split[samples$SAMPLE]
E_size <- vapply(E_list, length, integer(1))
if (any(E_size < samples$k)) {
  stop("[floor_derivation] a sample has fewer eligible windows than observed calls")
}

cluster_levels <- sort(unique(samples$Cluster))
c_idx <- match(samples$Cluster, cluster_levels)
n_by_cluster <- clusters_in %>% count(Cluster, name = "cluster_n")

# How much bigger is the eligible universe than it should be, because the
# consensus-tie drop is not reproduced here? Expressed as the median share of a
# sample's eligible windows — small = the caveat above is immaterial.
tie_inflation <- median(E_size)
say("universe: %d windows total; per-sample eligible median %d (min %d, max %d)",
    W, median(E_size), min(E_size), max(E_size))
say("per-sample observed calls k: median %d, max %d", median(samples$k), max(samples$k))

# --------------------------------------------------------------------------
# Observed support, through the same filter chain
# --------------------------------------------------------------------------
obs_w    <- match(obs$WINDOW, win_levels)
obs_c    <- match(obs$Cluster, cluster_levels)
pooled   <- tabulate(obs_w, nbins = W)          # distinct samples per window
keep_obs <- pooled > MIN_SAMP                   # filter 1, exactly as the pipeline
obs_support <- vapply(seq_along(cluster_levels), function(ci)
  tabulate(obs_w[obs_c == ci], nbins = W), integer(W))
obs_support[!keep_obs, ] <- 0L
say("observed windows surviving filter 1 (n > %d): %d of %d",
    MIN_SAMP, sum(keep_obs), W)

# --------------------------------------------------------------------------
# The permutation
# --------------------------------------------------------------------------
CAND_N <- seq_len(MAX_N)
null_hits <- matrix(0, nrow = MAX_N, ncol = length(cluster_levels),
                    dimnames = list(NULL, cluster_levels))
# Support histogram per cluster, over windows that survive filter 1 in that
# replicate. Bounded by cluster size, so a fixed-width accumulator is exact.
max_sup <- max(as.integer(n_by_cluster$cluster_n)) + 1L
sup_hist <- matrix(0, nrow = max_sup + 1L, ncol = length(cluster_levels),
                   dimnames = list(0:max_sup, cluster_levels))

set.seed(SEED)
t1 <- tick()
n_s <- nrow(samples)
k_v <- samples$k
for (p in seq_len(NPERM)) {
  draws <- vector("list", n_s)
  for (s in seq_len(n_s)) {
    # Uniform without replacement from THIS sample's eligible windows.
    draws[[s]] <- E_list[[s]][sample.int(E_size[s], k_v[s])]
  }
  w_all <- unlist(draws, use.names = FALSE)
  c_all <- rep.int(c_idx, k_v)

  pooled_p <- tabulate(w_all, nbins = W)
  keep_p   <- pooled_p > MIN_SAMP               # filter 1 on the null replicate

  for (ci in seq_along(cluster_levels)) {
    sup <- tabulate(w_all[c_all == ci], nbins = W)
    sup[!keep_p] <- 0L
    sup_hist[, ci] <- sup_hist[, ci] + tabulate(sup + 1L, nbins = max_sup + 1L)
    # Windows reaching each candidate N = reverse-cumulative of the support
    # histogram. nbins must cover the realised maximum: tabulate() silently
    # DISCARDS values above nbins, which would drop the very tail the whole
    # derivation is about whenever MAX_N sits below it.
    tb <- tabulate(sup, nbins = max(MAX_N, max(sup)))
    null_hits[, ci] <- null_hits[, ci] + rev(cumsum(rev(tb)))[seq_len(MAX_N)]
  }
  if (p %% 100 == 0) say("  permutation %d/%d (%.0fs)", p, NPERM, tick() - t1)
}
say("permutation done: %d replicates in %.1f min", NPERM, (tick() - t1) / 60)

# --------------------------------------------------------------------------
# Summarise
# --------------------------------------------------------------------------
# Null support distribution per cluster, over the filter-1-surviving windows
# (support 0 means "did not survive filter 1 in that replicate" and is excluded
# so the summary describes the population filter 2 actually sees).
null_summary <- map_dfr(seq_along(cluster_levels), function(ci) {
  h <- sup_hist[, ci]
  vals <- as.integer(names(h))
  h[1] <- 0                                     # drop support 0
  tot <- sum(h)
  if (tot == 0) return(tibble(Cluster = cluster_levels[ci], null_mean = NA_real_,
                              null_p95 = NA_real_, null_p99 = NA_real_, null_max = NA_real_))
  cw <- cumsum(h) / tot
  q  <- function(pr) vals[which(cw >= pr)[1]]
  tibble(Cluster   = cluster_levels[ci],
         null_mean = sum(vals * h) / tot,
         null_p95  = q(0.95), null_p99 = q(0.99),
         null_max  = max(vals[h > 0]))
})

tbl <- expand_grid(Cluster = cluster_levels, N = CAND_N) %>%
  mutate(ci = match(Cluster, cluster_levels)) %>%
  mutate(
    expected_null_windows = null_hits[cbind(N, ci)] / NPERM,
    observed_windows      = map2_int(ci, N, ~ sum(obs_support[, .x] >= .y)),
    fdr                   = if_else(observed_windows > 0,
                                    expected_null_windows / observed_windows, NA_real_)
  ) %>%
  dplyr::select(-ci) %>%
  left_join(n_by_cluster, by = "Cluster") %>%
  left_join(null_summary, by = "Cluster") %>%
  relocate(cluster_n, .after = Cluster)

# The chosen N: smallest N meeting BOTH targets in EVERY cluster.
meets <- tbl %>%
  group_by(N) %>%
  summarise(worst_expected = max(expected_null_windows),
            worst_fdr      = max(fdr, na.rm = TRUE),
            ok = all(expected_null_windows < TARGET_FW) &&
                 all(fdr < FDR_TGT | is.na(fdr)),
            .groups = "drop")
chosen <- meets %>% filter(ok) %>% slice_min(N, n = 1, with_ties = FALSE) %>% pull(N)
if (length(chosen) == 0) {
  say("WARNING: no N <= %d meets the targets; raise --max-n", MAX_N)
  chosen <- NA_integer_
}
binding <- if (is.na(chosen)) NA_character_ else {
  tbl %>% filter(N == chosen) %>% slice_max(expected_null_windows, n = 1,
                                            with_ties = FALSE) %>% pull(Cluster)
}

tbl <- tbl %>%
  mutate(chosen_N = chosen,
         is_chosen = !is.na(chosen) & N == chosen,
         binding_cluster = binding,
         universe_windows = W,
         eligible_windows_per_sample_median = tie_inflation,
         permutations = NPERM, seed = SEED,
         target_false_windows = TARGET_FW, fdr_target = FDR_TGT,
         min_samples_per_window = MIN_SAMP) %>%
  arrange(Cluster, N)
write_tsv(tbl, OUT)

cat("\n=============== per-cluster floor derivation ===============\n")
cat(sprintf("null: %d permutations, seed %d; each sample redraws its %s calls\n",
            NPERM, SEED, "observed number of"))
cat(sprintf("      uniformly from its own eligible window set (universe %d windows)\n", W))
cat(sprintf("target: expected null windows < %g per cluster AND FDR < %g\n\n",
            TARGET_FW, FDR_TGT))
cat("--- null per-window support, per cluster (post filter 1) ---\n")
print(as.data.frame(null_summary %>% left_join(n_by_cluster, by = "Cluster") %>%
                      relocate(cluster_n, .after = Cluster) %>%
                      mutate(across(where(is.numeric), ~ round(.x, 2)))))
cat("\n--- expected null windows / observed windows / FDR, by N ---\n")
print(as.data.frame(
  tbl %>% filter(N <= min(MAX_N, ifelse(is.na(chosen), MAX_N, chosen + 4))) %>%
    dplyr::select(Cluster, N, expected_null_windows, observed_windows, fdr) %>%
    mutate(expected_null_windows = round(expected_null_windows, 3),
           fdr = signif(fdr, 3))))
cat("\n")
if (is.na(chosen)) {
  cat("NO N MEETS THE TARGET within --max-n; nothing to recommend.\n")
} else {
  cat(sprintf("CHOSEN N = %d   (binding cluster: %s, n = %d)\n", chosen, binding,
              n_by_cluster$cluster_n[n_by_cluster$Cluster == binding]))
  cat(sprintf("  worst-case expected null windows at N: %.3f  (target < %g)\n",
              meets$worst_expected[meets$N == chosen], TARGET_FW))
  cat(sprintf("  worst-case FDR at N: %.4f  (target < %g)\n",
              meets$worst_fdr[meets$N == chosen], FDR_TGT))
}
cat(sprintf("\nwrote %s\n", OUT))
cat("============================================================\n")
