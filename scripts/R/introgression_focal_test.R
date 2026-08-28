#!/usr/bin/env Rscript
# introgression_focal_test.R — is a window's introgression ENRICHED in a focal
# subgroup, relative to the rest of its own cluster? (spec §9.6)
# --------------------------------------------------------------------------
# WHAT: for one focal subgroup (a value of a metadata role, e.g. geography ==
#       "Aceh"), test every window that its parent cluster carries: are more
#       of the subgroup's members called introgressed there than a random
#       subgroup of the same size would be? Per-window p, BH-adjusted across
#       the cluster's windows, called at `focal_fdr`.
#
# WHY:  Stage 5c derived a per-cluster support floor from a cluster-wide FDR
#       null, and showed it is structurally incompatible with a focal-group
#       claim: a 10-sample subgroup inside a 35-sample cluster can never reach
#       a floor calibrated on a 410-sample cluster (spec §9.5). That is not a
#       threshold to tune — it is the wrong test. The cluster-level question
#       ("which windows show gene flow between clusters") and the focal-level
#       question ("is this window's introgression concentrated in this
#       subgroup") are two different tests at two different scales, and each
#       needs its own null.
#
#       This is the focal-level test, and it is scaled to the subgroup. It is
#       a general design, not a rescue for one cohort: any cohort with a focal
#       subgroup — a district, a host species, a time window — needs its focal
#       contrast scaled to that subgroup rather than filtered through its
#       parent cluster's floor. The two tests are independent; nothing here
#       changes, weakens or reinterprets `per_cluster_min_samples`.
#
# THE STATISTIC. Per window w, `focal_support` = the number of FOCAL samples
# called introgressed at w. `background_support` = the number of the cluster's
# other members called at w.
#
# THE NULL (size-preserving label permutation). Hold each window's set of
# introgressed samples inside the cluster FIXED, and randomly reassign which
# `n_focal` of the cluster's members carry the focal label. Recompute
# per-window focal support. The question this poses is: *given the cluster's
# own introgression pattern, would a random subgroup of this size show this
# much support at this window by chance?*
#
#   - a window introgressed cluster-wide shows NO focal enrichment (a random
#     subgroup gets the same support);
#   - a window carried mostly by the focal group does.
#
# This SUBSUMES V1's "windows unique to the focal group" as the degenerate case
# where background support is exactly 0, but unlike "unique" it is not
# destroyed by a single background sample carrying the window.
#
# TWO P-VALUES, ONE NULL — read this before citing a number. Permuting the
# focal/rest label with the window's called set held fixed makes each window's
# marginal null exactly hypergeometric: focal_support ~ Hyper(cluster_n,
# n_called, n_focal). So the permutation and `phyper()` describe the SAME null,
# and `phyper()` evaluates it exactly instead of sampling it. Both are
# reported:
#   `p_perm` — the Monte-Carlo estimate, (b + 1) / (B + 1), `--permutations`
#              replicates at `--seed`. Its resolution floor is 1/(B+1), which
#              at B = 1000 is ~1e-3 — after BH across several hundred windows
#              that floor alone can make significance unreachable.
#   `p_raw`  — the exact tail of the same null. This is what BH adjusts, and
#              `p_adj` is what the call is made on.
# `p_perm` is kept because it is the specified procedure and because it is a
# live check on the analytic form: the run prints the largest discrepancy
# between the two over windows where the Monte-Carlo has resolution. A large
# discrepancy means the two disagree and the result should not be trusted.
#
# WHICH CALLS THE TEST SEES. The cached per-pair calls (detection is NOT
# re-run), restricted to the parent cluster's samples, with the two ARTIFACT
# masks applied and NO support floor of any kind:
#   - gene-family mask (filter 3) and hypervariable/multi-cluster mask
#     (filter 4) are read from the files `introgression_aggregate.R` writes.
#     Both excise things that are not introgression at all, so removing them
#     is not a support judgement and the focal test inherits them as-is. Note
#     that "hypervariable" is a cluster-scale designation, computed by the
#     aggregate step downstream of the cluster floor — it is taken from the
#     cluster-scale analysis on purpose, because that is the scale at which
#     "this window is variable in more than one cluster" is a meaningful claim.
#   - `min_samples_per_window` (filter 1) and `per_cluster_min_samples`
#     (filter 2) are NOT applied. Both are support floors, and this test has
#     its own multiple-testing control. Applying the cluster floor here is the
#     mis-posed test §9.5 diagnosed.
#
# THE TEST UNIVERSE (the BH denominator) is every window with at least one
# call among the cluster's samples after the masks. Windows the focal group
# does not touch are still tested (they score p = 1) — that is conservative
# and it is the honest denominator, since the test asks "which of this
# cluster's windows are focal-enriched".
#
# STATED LIMITATION — call-rate imbalance, and which way it bites. The label
# permutation treats every cluster member as equally likely to be called
# anywhere, so it does NOT preserve per-sample call rate: the null assumes the
# focal group's call rate is the cluster average. When it is not, the direction
# decides whether the result is threatened or reinforced:
#   focal called MORE often than average — the null under-states the focal
#     support a random subgroup would show, so the test OVER-calls. This is the
#     dangerous direction: a focal group with deeper coverage looks enriched
#     everywhere. Treat calls as provisional; a rate-stratified null is the fix.
#   focal called LESS often than average — the null over-states it, so the test
#     UNDER-calls. Surviving windows cleared a higher bar than their p-value
#     implies, and an empty result may be a power limit rather than an absence.
# The run prints the focal vs background per-sample call-rate summary and names
# the direction. Read it before reading the result.
#
# CLI (flags first, then one or more per-pair call TSVs as bare arguments):
#   Rscript scripts/R/introgression_focal_test.R \
#     --clusters     outputs/structure/admix_clusters.tsv \
#     --metadata     outputs/metadata/samples.tsv \
#     --contig-map   outputs/setup/contig_map.tsv \
#     --focal-group  Aceh \
#     --focal-role   geography \
#     --gene-family-mask outputs/introgression/gene_family_masked_windows.tsv \
#     --hypervariable-mask outputs/introgression/hypervariable_masked_windows.tsv \
#     --window-size  10000 \
#     --permutations 1000 \
#     --seed         20260828 \
#     --fdr          0.05 \
#     --out-enriched outputs/introgression/focal_Aceh_enriched_windows.tsv \
#     --out-tests    outputs/introgression/focal_Aceh_window_tests.tsv \
#     outputs/introgression/pairs/*.tsv
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

script_args <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- dirname(sub("^--file=", "",
                          script_args[grep("^--file=", script_args)][1]))
source(file.path(SCRIPT_DIR, "introgression_focal_core.R"))

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
required <- c("clusters", "metadata", "contig-map", "focal-group",
              "out-enriched", "out-tests")
missing  <- setdiff(required, names(args))
if (length(missing) > 0) stop("Missing args: ", paste(missing, collapse = ", "))

FOCAL     <- args[["focal-group"]]
ROLE      <- args[["focal-role"]] %||% "geography"
WINDOW_SZ <- as.integer(args[["window-size"]] %||% 10000)
NPERM     <- as.integer(args[["permutations"]] %||% 1000)
SEED      <- as.integer(args[["seed"]] %||% 20260828)
FDR       <- as.numeric(args[["fdr"]] %||% 0.05)
OUT_ENR   <- args[["out-enriched"]]
OUT_TESTS <- args[["out-tests"]]
for (p in c(OUT_ENR, OUT_TESTS)) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)

say <- function(...) message(sprintf("[introgression_focal_test] %s", sprintf(...)))

is_null_arg <- function(x) is.null(x) || is.na(x) || !nzchar(x) ||
                           x %in% c("NULL", "None", "null")

# Column contracts, so an empty result is a well-formed table rather than a
# missing file. A cohort with no focal enrichment is a result, not a crash.
TEST_COLS <- c("WINDOW", "CHROM", "CONTIG", "BIN", "START", "END",
               "n_called_cluster", "focal_support", "background_support",
               "focal_frac", "background_frac", "p_perm", "p_raw", "p_adj",
               "enriched", "Cluster", "focal_group", "focal_role", "n_focal",
               "n_background", "n_windows_tested", "permutations", "seed",
               "fdr_target")

write_empty_and_quit <- function(reason) {
  say("%s — writing header-only outputs", reason)
  empty <- as_tibble(setNames(
    replicate(length(TEST_COLS), character(0), simplify = FALSE), TEST_COLS))
  write_tsv(empty, OUT_ENR)
  write_tsv(empty, OUT_TESTS)
  quit(status = 0)
}

if (length(call_files) == 0) write_empty_and_quit("no per-pair call files given")

# --------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------
clusters_in <- read_tsv(args[["clusters"]], show_col_types = FALSE) %>%
  dplyr::select(SAMPLE = Sample, Cluster) %>%
  distinct()

meta <- read_tsv(args[["metadata"]], show_col_types = FALSE)
if (!("sample_id" %in% names(meta))) stop("metadata has no sample_id column")
if (!(ROLE %in% names(meta))) {
  write_empty_and_quit(sprintf("role '%s' absent from the metadata — cannot resolve focal_group", ROLE))
}

cmap <- read_tsv(args[["contig-map"]], col_names = c("CONTIG", "CHROM"),
                 show_col_types = FALSE)

membership <- meta %>%
  dplyr::select(SAMPLE = sample_id, role_value = all_of(ROLE)) %>%
  inner_join(clusters_in, by = "SAMPLE") %>%
  mutate(side = if_else(!is.na(role_value) & role_value == FOCAL, "focal", "background"))

focal_members <- membership %>% filter(side == "focal")
if (nrow(focal_members) == 0) {
  write_empty_and_quit(sprintf("no clustered samples with %s == '%s'", ROLE, FOCAL))
}

# -- scope: the focal group must sit inside ONE cluster --------------------
# A subgroup that straddles clusters has no single "rest of its cluster" to be
# contrasted against, so the test is undefined. The headline step resolves the
# same ambiguity by majority vote; here it is an error, because a permutation
# null built on the wrong comparison set is a wrong p-value rather than a
# rough one.
focal_by_cluster <- focal_members %>% count(Cluster, name = "n_focal") %>%
  arrange(desc(n_focal))
if (nrow(focal_by_cluster) > 1) {
  say("'%s' spans %d clusters: %s", FOCAL, nrow(focal_by_cluster),
      paste(sprintf("%s=%d", focal_by_cluster$Cluster, focal_by_cluster$n_focal),
            collapse = ", "))
  stop("[introgression_focal_test] focal group '", FOCAL, "' straddles ",
       nrow(focal_by_cluster), " clusters — the focal-vs-rest-of-cluster ",
       "contrast is undefined. Restrict the group, or run the test per cluster.")
}
focal_cluster <- focal_by_cluster$Cluster[1]

cohort <- membership %>% filter(Cluster == focal_cluster)
cluster_samples <- sort(cohort$SAMPLE)
N_CLUSTER  <- length(cluster_samples)
focal_set  <- cohort %>% filter(side == "focal") %>% pull(SAMPLE)
N_FOCAL    <- length(focal_set)
N_BACKGRND <- N_CLUSTER - N_FOCAL
say("cluster '%s': %d members = %d focal ('%s' by role '%s') + %d background",
    focal_cluster, N_CLUSTER, N_FOCAL, FOCAL, ROLE, N_BACKGRND)
if (N_BACKGRND == 0) {
  write_empty_and_quit(sprintf(
    "the focal group IS the whole cluster '%s' — there is nothing to contrast against",
    focal_cluster))
}

# --------------------------------------------------------------------------
# Calls: cached per-pair files, cluster-scoped, artifact masks applied,
# NO support floor. See the header for why filters 1 and 2 are excluded.
# --------------------------------------------------------------------------
calls <- call_files %>%
  map(~ read_tsv(.x, show_col_types = FALSE,
                 col_types = cols(WINDOW = col_character(), SAMPLE = col_character(),
                                  Cluster = col_character(), CHROM = col_integer(),
                                  BIN = col_double(), .default = col_guess()))) %>%
  bind_rows()
if (nrow(calls) == 0) write_empty_and_quit("no calls in the cached per-pair files")

read_mask <- function(path, label) {
  if (is_null_arg(path) || !file.exists(path)) {
    say("%s mask: not supplied — skipped", label); return(character(0))
  }
  w <- read_tsv(path, show_col_types = FALSE)$WINDOW
  w <- as.character(w[!is.na(w)])
  say("%s mask: %d window(s)", label, length(w))
  w
}
masked <- unique(c(read_mask(args[["gene-family-mask"]],   "gene-family"),
                   read_mask(args[["hypervariable-mask"]], "hypervariable")))

cluster_calls <- calls %>%
  filter(SAMPLE %in% cluster_samples, !(WINDOW %in% masked)) %>%
  distinct(WINDOW, CHROM, BIN, SAMPLE)
if (nrow(cluster_calls) == 0) {
  write_empty_and_quit(sprintf("no unmasked calls among cluster '%s' samples", focal_cluster))
}

# -- call-rate diagnostic (see the header's stated limitation) -------------
per_sample <- cluster_calls %>% count(SAMPLE, name = "n_windows") %>%
  right_join(tibble(SAMPLE = cluster_samples), by = "SAMPLE") %>%
  mutate(n_windows = replace_na(n_windows, 0L),
         side = if_else(SAMPLE %in% focal_set, "focal", "background"))
rate <- per_sample %>% group_by(side) %>%
  summarise(n_samples = dplyr::n(), mean = mean(n_windows), median = median(n_windows),
            min = min(n_windows), max = max(n_windows), .groups = "drop")
rate_ratio <- {
  m <- setNames(rate$mean, rate$side)
  if (all(c("focal", "background") %in% names(m)) && m[["background"]] > 0)
    m[["focal"]] / m[["background"]] else NA_real_
}

# --------------------------------------------------------------------------
# Incidence matrix: windows x cluster members
# --------------------------------------------------------------------------
win_levels <- sort(unique(cluster_calls$WINDOW))
W <- length(win_levels)
M <- matrix(0L, nrow = W, ncol = N_CLUSTER,
            dimnames = list(win_levels, cluster_samples))
M[cbind(match(cluster_calls$WINDOW, win_levels),
        match(cluster_calls$SAMPLE, cluster_samples))] <- 1L

focal_cols <- match(focal_set, cluster_samples)
say("test universe: %d window(s) with >= 1 call in cluster '%s' (%d call(s) total)",
    W, focal_cluster, nrow(cluster_calls))

# --------------------------------------------------------------------------
# The null — both forms, in scripts/R/introgression_focal_core.R
# --------------------------------------------------------------------------
stats <- focal_enrichment(M, focal_cols, nperm = NPERM, seed = SEED)
p_raw <- stats$p_exact
p_adj <- p.adjust(p_raw, method = "BH")

# Agreement check between the sampled and the exact form of the null, over the
# windows where the Monte Carlo has resolution to speak (p above ~5 counts).
resolvable <- stats$p_perm > 5 / (NPERM + 1)
max_disc <- if (any(resolvable))
  max(abs(stats$p_perm[resolvable] - p_raw[resolvable])) else NA_real_

# --------------------------------------------------------------------------
# Assemble + write
# --------------------------------------------------------------------------
coords <- cluster_calls %>% distinct(WINDOW, CHROM, BIN) %>%
  mutate(START = as.integer(BIN - WINDOW_SZ %/% 2),
         END   = as.integer(BIN - WINDOW_SZ %/% 2 + WINDOW_SZ - 1)) %>%
  left_join(cmap, by = "CHROM")

tests <- stats %>%
  dplyr::select(WINDOW, n_called_cluster, focal_support, background_support, p_perm) %>%
  mutate(p_raw = p_raw, p_adj = p_adj) %>%
  mutate(focal_frac      = focal_support / N_FOCAL,
         background_frac = background_support / N_BACKGRND,
         enriched        = p_adj < FDR) %>%
  left_join(coords, by = "WINDOW") %>%
  mutate(Cluster = focal_cluster, focal_group = FOCAL, focal_role = ROLE,
         n_focal = N_FOCAL, n_background = N_BACKGRND, n_windows_tested = W,
         permutations = NPERM, seed = SEED, fdr_target = FDR) %>%
  dplyr::select(all_of(TEST_COLS)) %>%
  arrange(p_adj, desc(focal_support), CHROM, BIN)

write_tsv(tests, OUT_TESTS)
enriched <- tests %>% filter(enriched)
write_tsv(enriched, OUT_ENR)

# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------
cat("\n============ focal-group enrichment test ============\n")
cat(sprintf("focal group : %s == '%s'  (n = %d)\n", ROLE, FOCAL, N_FOCAL))
cat(sprintf("contrast    : the rest of cluster '%s'  (n = %d)\n", focal_cluster, N_BACKGRND))
cat(sprintf("null        : size-preserving focal/background label permutation,\n"))
cat(sprintf("              %d replicates, seed %d; exact tail via phyper()\n", NPERM, SEED))
cat(sprintf("universe    : %d window(s), BH-adjusted, called at FDR < %g\n", W, FDR))
cat(sprintf("filters     : gene-family + hypervariable masks applied; NO support\n"))
cat(sprintf("              floor (the cluster floor is a separate, cluster-scaled test)\n\n"))

cat("--- per-sample call rate, focal vs background (read this first) ---\n")
print(as.data.frame(rate %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))
if (!is.na(rate_ratio)) {
  cat(sprintf("focal / background mean call rate = %.2f\n", rate_ratio))
  # Direction matters, and it decides whether the imbalance is a threat to the
  # result or a reason to trust it more. The null draws a random subgroup from
  # the whole cluster, so it assumes the focal group's call rate is the cluster
  # average.
  #   ratio > 1 (focal called MORE often): the null under-states expected focal
  #     support, the test over-calls — ANTI-CONSERVATIVE, treat calls as
  #     provisional.
  #   ratio < 1 (focal called LESS often): the null over-states it, the test
  #     under-calls — CONSERVATIVE, so surviving calls are if anything harder
  #     to reach than the p-value says, and an empty result may be a power
  #     problem rather than an absence.
  if (rate_ratio > 1.25) {
    cat("  ⚠ MATERIAL IMBALANCE, ANTI-CONSERVATIVE DIRECTION. The label permutation\n")
    cat("    does not preserve per-sample call rate, and the focal group is called\n")
    cat("    more often than the cluster average — so it is over-represented at\n")
    cat("    EVERY window and the test will read that as enrichment. Treat the calls\n")
    cat("    below as provisional; a rate-stratified null is the fix. See the stated\n")
    cat("    limitation in this script's header.\n")
  } else if (rate_ratio < 0.8) {
    cat("  MATERIAL IMBALANCE, CONSERVATIVE DIRECTION. The label permutation does not\n")
    cat("    preserve per-sample call rate, and the focal group is called LESS often\n")
    cat("    than the cluster average — so the null over-states the focal support a\n")
    cat("    random subgroup would show, and the test UNDER-calls. Any window that\n")
    cat("    survives cleared a bar higher than its p-value implies; an empty result\n")
    cat("    may be a power limit rather than an absence. See the header.\n")
  }
}

cat(sprintf("\nnull form check: max |p_perm - p_raw| over resolvable windows = %s\n",
            if (is.na(max_disc)) "n/a (no window above the MC resolution floor)"
            else sprintf("%.4f", max_disc)))
if (!is.na(max_disc) && max_disc > 0.05) {
  cat("  ⚠ the sampled and exact forms of the null disagree — investigate before citing.\n")
}

cat(sprintf("\nRESULT: %d of %d window(s) focal-enriched in '%s' at FDR < %g\n",
            nrow(enriched), W, FOCAL, FDR))
if (nrow(enriched) > 0) {
  enriched %>% head(15) %>%
    mutate(line = sprintf("  %-14s %s:%d-%d  focal %d/%d  background %d/%d  p_raw=%.3g  p_adj=%.3g",
                          WINDOW, CONTIG, START, END, focal_support, N_FOCAL,
                          background_support, N_BACKGRND, p_raw, p_adj)) %>%
    pull(line) %>% cat(sep = "\n")
  cat("\n")
} else {
  cat("  An empty result is a result: no window in this cluster carries more focal\n")
  cat("  support than a random subgroup of the same size would.\n")
  best <- tests %>% slice_min(p_adj, n = 3, with_ties = FALSE)
  cat("  Closest windows (not significant):\n")
  best %>%
    mutate(line = sprintf("  %-14s %s:%d-%d  focal %d/%d  background %d/%d  p_raw=%.3g  p_adj=%.3g",
                          WINDOW, CONTIG, START, END, focal_support, N_FOCAL,
                          background_support, N_BACKGRND, p_raw, p_adj)) %>%
    pull(line) %>% cat(sep = "\n")
  cat("\n")
}
cat(sprintf("\nwrote %s\n     %s\n", OUT_ENR, OUT_TESTS))
cat("=====================================================\n")
