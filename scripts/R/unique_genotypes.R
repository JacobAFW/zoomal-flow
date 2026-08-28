#!/usr/bin/env Rscript
# unique_genotypes.R — collapse clonal groups to one representative each
# (docs/clonality.md)
# --------------------------------------------------------------------------
# WHAT: from Stage 4's `clonal_clusters.tsv`, write the de-clonalized sample
#       keep-list: one representative per `clonal_group`, plus every sample
#       that belongs to no group. Also writes the full-set keep-list, so the
#       `full` and `unique` sample sets are the same KIND of object and the
#       downstream rules are genuinely one code path with a swapped file.
#       (The full-set list is normally written by the `all_genotypes` rule
#       instead, straight from the .fam: it MUST NOT depend on Stage 4, or the
#       DAG would cycle — `full` structure feeds Stage 4, which feeds the
#       unique keep-list. `--out-full` here is for standalone use.)
#
# WHY:  clonal samples are near-identical genotypes. Counting each as
#       independent evidence is pseudo-replication, and it biases every
#       frequency- or count-based analysis — allele frequencies, ADMIXTURE
#       ancestry proportions, per-group polyclonality rates, introgression
#       window support, iHS case/control frequencies. A clonal block of 3 in a
#       35-sample cluster is one genotype speaking three times. Plasmodium
#       cohorts are frequently clonal, so a pop-gen pipeline that ignores this
#       is scientifically incomplete.
#
# THE REPRESENTATIVE RULE. Deterministic, and stated so a collaborator can
# reproduce the choice by hand:
#   1. lowest per-sample missingness (`F_MISS` from PLINK's .smiss), then
#   2. alphabetical sample id, as a tie-break.
# This mirrors the Stage-3 duplicate pick (`find_duplicates.R`), which keeps
# the lowest-missingness member of a replicate pair — same question, same
# answer, so the pipeline is not making two different "which copy is best"
# judgements. Missingness is optional: without a .smiss the rule degrades to
# alphabetical only, which is still deterministic, and it is logged loudly.
#
# STRADDLING CLONAL GROUPS. A clonal group whose members sit in more than one
# ADMIXTURE cluster is a real possibility (clonality is derived per cluster
# here, so it should not arise on this pipeline's own output — but the script
# does not assume its input came from this pipeline). The representative is
# taken from the MAJORITY cluster; a tie between clusters is broken by the
# alphabetically-first cluster name and logged. The count of straddling groups
# is reported so it is never silent.
#
# WHAT THIS IS NOT. It is not a re-derivation of clonality, and it is not
# iterated: the keep-list is computed ONCE from the full-set clonal groups and
# used as-is. Re-deriving clonality on the de-clonalized set would be circular
# and is the field norm to avoid — clonal groups are robust to it. The
# full-set clusters stay the operational backbone; see docs/clonality.md.
#
# CLI:
#   Rscript scripts/R/unique_genotypes.R \
#     --clonal-clusters outputs/ibd/clonal_clusters.tsv \
#     --smiss           outputs/structure/Pk.smiss      # optional \
#     --all-samples     outputs/structure/Pk.fam        # defines the FULL set \
#     --out-unique      outputs/clonality/unique_genotypes.txt \
#     --out-exclude     outputs/clonality/exclude_unique.txt \
#     --out-full        outputs/clonality/all_genotypes.txt   # optional \
#     --out-audit       outputs/clonality/declonalization_audit.tsv
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

parse_args <- function(av) {
  out <- list(); i <- 1
  while (i <= length(av)) {
    if (!startsWith(av[i], "--")) stop("Bad arg: ", av[i])
    out[[sub("^--", "", av[i])]] <- av[i + 1]; i <- i + 2
  }
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
`%||%` <- function(a, b) if (is.null(a)) b else a
required <- c("clonal-clusters", "all-samples", "out-unique", "out-exclude", "out-audit")
missing  <- setdiff(required, names(args))
if (length(missing) > 0) stop("Missing args: ", paste(missing, collapse = ", "))

is_null_arg <- function(x) is.null(x) || is.na(x) || !nzchar(x) ||
                           x %in% c("NULL", "None", "null")

for (p in c(args[["out-unique"]], args[["out-full"]], args[["out-exclude"]],
            args[["out-audit"]])) {
  if (!is.null(p)) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
}
say <- function(...) message(sprintf("[unique_genotypes] %s", sprintf(...)))

# --------------------------------------------------------------------------
# The FULL sample set — the universe both keep-lists are drawn from.
# Read from the .fam so it is exactly the set the frequency stages can see,
# not the metadata table (which may list samples that never made it through QC).
# --------------------------------------------------------------------------
fam <- read.table(args[["all-samples"]], header = FALSE, stringsAsFactors = FALSE)
all_samples <- sort(unique(as.character(fam$V1)))
say("full sample set: %d samples (from %s)", length(all_samples), args[["all-samples"]])

write_keep <- function(samples, path) {
  writeLines(sort(unique(samples)), path)
}
# The FULL keep-list is deliberately derivable WITHOUT Stage 4 — see the
# `all_genotypes` rule. Written here only if a path is given, so this script
# stays usable standalone.
if (!is.null(args[["out-full"]])) write_keep(all_samples, args[["out-full"]])

# --------------------------------------------------------------------------
# Clonal groups
# --------------------------------------------------------------------------
clonal <- read_tsv(args[["clonal-clusters"]], show_col_types = FALSE,
                   col_types = cols(.default = col_character()))
stopifnot(all(c("sample", "cluster", "clonal_group") %in% names(clonal)))

clonal <- clonal %>%
  mutate(clonal_group = if_else(clonal_group %in% c("NA", ""), NA_character_,
                                clonal_group)) %>%
  # Only samples that survive into the analysis set can be kept or dropped.
  filter(sample %in% all_samples)

grouped <- clonal %>% filter(!is.na(clonal_group))
n_groups <- n_distinct(grouped$clonal_group)
say("clonal groups: %d, covering %d samples", n_groups, nrow(grouped))

if (n_groups == 0) {
  say("no clonal groups — the unique set equals the full set")
  write_keep(all_samples, args[["out-unique"]])
  writeLines(character(0), args[["out-exclude"]])
  write_tsv(tibble(clonal_group = character(), cluster = character(),
                   n_members = integer(), representative = character(),
                   dropped = character(), straddles_clusters = logical(),
                   representative_f_miss = numeric()),
            args[["out-audit"]])
  cat(sprintf("\n[unique_genotypes] no clonal signal: unique set = full set (%d samples)\n",
              length(all_samples)))
  quit(status = 0)
}

# --------------------------------------------------------------------------
# Missingness, for the representative rule
# --------------------------------------------------------------------------
smiss_path <- args[["smiss"]]
if (is_null_arg(smiss_path) || !file.exists(smiss_path)) {
  say("WARNING: no .smiss supplied — the representative rule degrades to")
  say("         ALPHABETICAL ONLY. Deterministic, but not missingness-aware.")
  miss <- tibble(sample = all_samples, F_MISS = NA_real_)
} else {
  miss <- read_table(smiss_path, show_col_types = FALSE) %>%
    rename(any_of(c(FID = "#FID"))) %>%
    dplyr::select(sample = IID, F_MISS) %>%
    mutate(F_MISS = as.numeric(F_MISS))
  say("missingness read for %d samples (median F_MISS %.4f)",
      nrow(miss), median(miss$F_MISS, na.rm = TRUE))
}

# --------------------------------------------------------------------------
# Pick one representative per group
# --------------------------------------------------------------------------
# Straddling: restrict the candidates to the group's MAJORITY cluster before
# the missingness pick. Ties between clusters go to the alphabetically-first
# cluster name, so the result never depends on row order.
majority_cluster <- grouped %>%
  count(clonal_group, cluster, name = "n_in_cluster") %>%
  arrange(clonal_group, desc(n_in_cluster), cluster) %>%
  group_by(clonal_group) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(clonal_group, majority_cluster = cluster)

straddling <- grouped %>%
  distinct(clonal_group, cluster) %>%
  count(clonal_group, name = "n_clusters") %>%
  filter(n_clusters > 1)
if (nrow(straddling) > 0) {
  say("WARNING: %d clonal group(s) straddle more than one cluster — the",
      nrow(straddling))
  say("         representative is taken from the majority cluster: %s",
      paste(straddling$clonal_group, collapse = ", "))
} else {
  say("no straddling clonal groups (every group sits in one cluster)")
}

reps <- grouped %>%
  left_join(majority_cluster, by = "clonal_group") %>%
  left_join(miss, by = "sample") %>%
  filter(cluster == majority_cluster) %>%
  # NA missingness sorts last, so a sample with a known F_MISS always wins
  # over one without; alphabetical breaks the remaining ties.
  arrange(clonal_group, F_MISS, sample) %>%
  group_by(clonal_group) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(clonal_group, representative = sample,
                representative_cluster = cluster, representative_f_miss = F_MISS)

dropped <- grouped %>%
  anti_join(reps, by = c("sample" = "representative")) %>%
  pull(sample)

unique_samples <- setdiff(all_samples, dropped)
write_keep(unique_samples, args[["out-unique"]])
# The EXCLUSION list is what the frequency stages actually consume — see
# sampleset_exclude() in the Snakefile for why a drop-list and not a keep-list.
write_keep(dropped, args[["out-exclude"]])

# --------------------------------------------------------------------------
# Audit: one row per clonal group, plus the per-cluster effective n
# --------------------------------------------------------------------------
audit <- grouped %>%
  group_by(clonal_group) %>%
  summarise(cluster    = paste(sort(unique(cluster)), collapse = ";"),
            n_members  = dplyr::n(),
            members    = paste(sort(sample), collapse = ";"),
            .groups    = "drop") %>%
  left_join(reps, by = "clonal_group") %>%
  left_join(straddling, by = "clonal_group") %>%
  mutate(straddles_clusters = !is.na(n_clusters),
         dropped = map2_chr(members, representative,
                            ~ paste(setdiff(strsplit(.x, ";")[[1]], .y), collapse = ";")),
         n_dropped = n_members - 1L) %>%
  dplyr::select(clonal_group, cluster, n_members, n_dropped, representative,
                representative_cluster, representative_f_miss, dropped,
                straddles_clusters) %>%
  arrange(desc(n_members), clonal_group)
write_tsv(audit, args[["out-audit"]])

# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------
per_cluster <- clonal %>%
  mutate(kept = sample %in% unique_samples) %>%
  group_by(cluster) %>%
  summarise(n_full = dplyr::n(), n_unique = sum(kept), .groups = "drop") %>%
  mutate(n_dropped = n_full - n_unique) %>%
  arrange(desc(n_full))

cat("\n================ de-clonalization ================\n")
cat(sprintf("representative rule: lowest F_MISS, ties broken alphabetically%s\n",
            if (all(is.na(miss$F_MISS))) "  [NO .smiss — ALPHABETICAL ONLY]" else ""))
cat(sprintf("clonal groups collapsed: %d  (%d samples -> %d representatives)\n",
            n_groups, nrow(grouped), n_groups))
cat(sprintf("straddling groups: %d\n", nrow(straddling)))
cat(sprintf("\nfull set:   %d samples\nunique set: %d samples  (%d dropped)\n",
            length(all_samples), length(unique_samples), length(dropped)))
cat("\n--- effective n per cluster ---\n")
print(as.data.frame(per_cluster))
cat("\n--- clonal groups (largest first) ---\n")
print(as.data.frame(audit %>% dplyr::select(clonal_group, cluster, n_members,
                                            representative, n_dropped) %>%
                      head(20)))
cat(sprintf("\nwrote %s\n      %s\n      %s\n", args[["out-unique"]],
            args[["out-exclude"]], args[["out-audit"]]))
cat("==================================================\n")
