#!/usr/bin/env Rscript
# introgression_headline.R — DESCRIPTIVE set difference: windows called on
# only one side of a focal/rest split (spec §6, §9.6)
# --------------------------------------------------------------------------
# ⚠ SUPERSEDED AS A RESULT — read before citing anything from this script.
#   This is NOT a statistical test. "Unique to the focal group" is a raw set
#   difference with no null, no p-value and no multiple-testing control, and it
#   is brittle in the worst way: one background sample carrying the window
#   flips a window out of the set entirely. It was V1's headline; it is now a
#   descriptive companion table only.
#
#   The result to cite is scripts/R/introgression_focal_test.R —
#   focal_<group>_enriched_windows.tsv — which asks the same question properly:
#   is a window's introgression ENRICHED in the focal subgroup relative to the
#   rest of its cluster, under a size-preserving label-permutation null scaled
#   to the subgroup, BH-adjusted at `focal_fdr`. Enrichment subsumes "unique"
#   as the degenerate case where background support is 0, and it survives a
#   single background carrier. See spec §9.6 for why the two were separated.
#
# WHAT: split one cluster's samples into "the focal group" and "the rest of
#       that cluster", then find the windows called in exactly one of the two
#       sides.
# WHY:  kept because the set difference is still worth eyeballing beside the
#       test, and because it is what V1 reported, so the two are comparable.
#       V1 hardcoded `State %in% c("Aceh", "Peninsular")` and compared two
#       literal geography values (legacy introgression_multi_cluster.R:583-632).
#       Here `focal_group` is a config value matched against the role named by
#       `--focal-role` (default `geography`), and the comparison set is
#       derived: whichever cluster the focal samples belong to, minus the focal
#       samples themselves. Any cohort, any group, no name literals.
#
# Ambiguity handled explicitly: if the focal group's samples span more than one
# cluster, the cluster holding the most of them is used and the split is
# logged. That is a judgement call the user should see, so it goes to the log
# and to focal_group_membership.tsv rather than being silently resolved.
#
# CLI:
#   Rscript introgression_headline.R \
#     --calls       outputs/introgression/introgressed_windows_filtered.tsv \
#     --clusters    outputs/structure/admix_clusters.tsv \
#     --metadata    outputs/metadata/samples.tsv \
#     --focal-group Aceh \
#     --focal-role  geography \
#     --out-dir     outputs/introgression \
#     --out-unique  outputs/introgression/unique_windows_in_Aceh_with_freq_and_coords.tsv \
#     --out-per-chrom outputs/introgression/unique_windows_per_chrom_Aceh_vs_rest.tsv
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
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
required <- c("calls", "clusters", "metadata", "focal-group", "out-dir",
              "out-unique", "out-per-chrom")
missing  <- setdiff(required, names(args))
if (length(missing) > 0) stop("Missing args: ", paste(missing, collapse = ", "))

focal    <- args[["focal-group"]]
role     <- if (is.null(args[["focal-role"]])) "geography" else args[["focal-role"]]
out_dir  <- args[["out-dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

UNIQUE_COLS   <- c("WINDOW", "CHROM", "CONTIG", "BIN", "START", "END", "n_samples")
PER_CHROM_COLS<- c("CHROM", "CONTIG", "side", "n_windows")

write_empty_and_quit <- function(reason) {
  message("[introgression_headline] ", reason, " — writing header-only outputs")
  empty <- function(cols) as_tibble(setNames(
    replicate(length(cols), character(0), simplify = FALSE), cols))
  write_tsv(empty(UNIQUE_COLS),    args[["out-unique"]])
  write_tsv(empty(PER_CHROM_COLS), args[["out-per-chrom"]])
  quit(status = 0)
}

calls <- read_tsv(args[["calls"]], show_col_types = FALSE,
                  col_types = cols(WINDOW = col_character(),
                                   SAMPLE = col_character(),
                                   Cluster = col_character(),
                                   .default = col_guess()))
if (nrow(calls) == 0) write_empty_and_quit("no filtered calls")

meta <- read_tsv(args[["metadata"]], show_col_types = FALSE)
if (!(role %in% names(meta))) {
  write_empty_and_quit(sprintf("role '%s' absent — cannot resolve focal_group", role))
}
clusters_in <- read_tsv(args[["clusters"]], show_col_types = FALSE) %>%
  dplyr::select(SAMPLE = Sample, Cluster) %>%
  distinct()

membership <- meta %>%
  dplyr::select(SAMPLE = sample_id, role_value = all_of(role)) %>%
  inner_join(clusters_in, by = "SAMPLE") %>%
  mutate(side = ifelse(!is.na(role_value) & role_value == focal, "focal", "rest"))

focal_samples <- membership %>% filter(side == "focal")
if (nrow(focal_samples) == 0) {
  write_empty_and_quit(sprintf("no samples with %s == '%s'", role, focal))
}

# Which cluster do the focal samples belong to?
focal_cluster_counts <- focal_samples %>% count(Cluster, name = "n_focal") %>%
  arrange(desc(n_focal))
focal_cluster <- focal_cluster_counts$Cluster[1]
if (nrow(focal_cluster_counts) > 1) {
  message(sprintf(paste0("[introgression_headline] '%s' samples span %d clusters (%s); ",
                         "using '%s' (n=%d) as the comparison cluster"),
                  focal, nrow(focal_cluster_counts),
                  paste(sprintf("%s:%d", focal_cluster_counts$Cluster,
                                focal_cluster_counts$n_focal), collapse = ", "),
                  focal_cluster, focal_cluster_counts$n_focal[1]))
}

cohort <- membership %>% filter(Cluster == focal_cluster)
n_focal <- sum(cohort$side == "focal")
n_rest  <- sum(cohort$side == "rest")
message(sprintf("[introgression_headline] cluster '%s': %d focal ('%s') vs %d rest",
                focal_cluster, n_focal, focal, n_rest))
write_tsv(cohort, file.path(out_dir, "focal_group_membership.tsv"))
if (n_rest == 0) {
  message("[introgression_headline] WARNING: the focal group IS the whole cluster — ",
          "every window is trivially 'unique'; interpret with care")
}

# -- windows unique to one side (V1 steps 1-5, generalised) ----------------
side_calls <- calls %>%
  inner_join(cohort %>% dplyr::select(SAMPLE, side), by = "SAMPLE")
if (nrow(side_calls) == 0) {
  write_empty_and_quit(sprintf("no calls among cluster '%s' samples", focal_cluster))
}

window_by_side <- side_calls %>% distinct(WINDOW, side)
unique_windows <- window_by_side %>%
  count(WINDOW, name = "n_sides") %>%
  filter(n_sides == 1) %>%
  left_join(window_by_side, by = "WINDOW")

coords <- side_calls %>%
  distinct(WINDOW, CHROM, CONTIG, BIN, START, END)

# Per-chromosome counts, both sides (V1's unique_windows_per_chr_*.tsv).
unique_windows %>%
  left_join(coords, by = "WINDOW") %>%
  count(CHROM, CONTIG, side, name = "n_windows") %>%
  arrange(side, CHROM) %>%
  write_tsv(args[["out-per-chrom"]])

# The headline table: focal-unique windows with sample frequency + coordinates.
focal_unique <- unique_windows %>% filter(side == "focal") %>% pull(WINDOW)

headline <- side_calls %>%
  filter(side == "focal", WINDOW %in% focal_unique) %>%
  distinct(SAMPLE, WINDOW) %>%
  count(WINDOW, name = "n_samples") %>%
  left_join(coords, by = "WINDOW") %>%
  dplyr::select(all_of(UNIQUE_COLS)) %>%
  arrange(desc(n_samples), CHROM, BIN)

write_tsv(headline, args[["out-unique"]])

cat(sprintf(paste0("\n[introgression_headline] DESCRIPTIVE (not a test): %d window(s) ",
                   "called only in '%s', not in the rest of cluster '%s' ",
                   "(%d focal vs %d rest samples)\n"),
            nrow(headline), focal, focal_cluster, n_focal, n_rest))
cat(paste0("  This is a set difference, not a result. No null, no p-value, no ",
           "multiple-testing\n  control, and one background carrier flips a window ",
           "out of the set. Cite the\n  focal ENRICHMENT test instead: ",
           "focal_", focal, "_enriched_windows.tsv (spec \u00a79.6).\n"))
if (nrow(headline) > 0) {
  headline %>% head(10) %>%
    mutate(line = sprintf("  %-12s %s:%d-%d  n=%d", WINDOW, CONTIG, START, END, n_samples)) %>%
    pull(line) %>% cat(sep = "\n")
  cat("\n")
}
