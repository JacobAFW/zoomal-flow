#!/usr/bin/env Rscript
# clonal_clusters.R  (agnostic — connected-components, no hardcoded tribble)
# --------------------------------------------------------------------------
# From every per-cluster hmmIBD hmm_fract.txt, identify clonal sample pairs
# (fract_sites_IBD >= threshold) and derive clonal groups as the connected
# components of the resulting IBD graph. Joins per-sample rows to canonical
# `geography` and `date` roles where present; falls back silently when a
# role is absent.
#
# Two outputs:
#   <out_clusters_tsv>   per-sample table (sample, cluster, clonal_group,
#                        + geography/date if roles present)
#   <focal_pairs_tsv>    per-pair table for a single focal cluster (V1's
#                        Pen_clones_IBD.txt shape). Only written if a real
#                        path is given; pass "NULL" to skip.
#
# Sample-ID cleanup: PLINK 2 --double-id / --recode vcf preserves sample IDs
# verbatim in modern versions. The V1 doubled-ID undo (`.undouble`) and
# lane-suffix strip were needed for legacy inputs; this port keeps a safe
# undouble as a no-op for undoubled IDs but drops the lane-suffix strip
# (input samples are already clean — see data/INPUT_PREP.md for the Indo
# cohort's pre-cleaning).
#
# Usage:
#   Rscript clonal_clusters.R <threshold> <metadata.tsv> \
#       <out_clusters_tsv> <focal_pairs_tsv|NULL> <focal_cluster|NULL> \
#       <fract1:cluster1> <fract2:cluster2> ...
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop(paste(
    "Usage: Rscript clonal_clusters.R <threshold> <metadata.tsv>",
    "<out_clusters_tsv> <focal_pairs_tsv|NULL> <focal_cluster|NULL>",
    "<fract1:cluster1> [<fract2:cluster2> ...]"
  ))
}
thresh         <- as.numeric(args[1])
meta_in        <- args[2]
out_clusters   <- args[3]
focal_pairs_out<- args[4]
focal_cluster  <- args[5]
fract_specs    <- args[-(1:5)]

meta <- read_tsv(meta_in, show_col_types = FALSE)
stopifnot("sample_id" %in% names(meta))
meta_keep <- meta %>%
  select(any_of(c("sample_id", "country", "geography", "date")))

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

# Parse each `path:cluster` spec, load its fract, tag rows with the cluster.
all_ibd <- map_dfr(fract_specs, function(spec) {
  bits <- strsplit(spec, ":", fixed = TRUE)[[1]]
  if (length(bits) != 2) stop("Bad spec (want path:cluster): ", spec)
  path <- bits[1]; cluster <- bits[2]
  read.delim(path, stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    select(sample1, sample2, fract_sites_IBD) %>%
    mutate(
      cluster = cluster,
      sample1 = .undouble(sample1),
      sample2 = .undouble(sample2)
    )
})

clonal_pairs <- all_ibd %>%
  filter(fract_sites_IBD >= thresh)
message(sprintf("[clonal_clusters] %d clonal pairs across %d cluster(s)",
                nrow(clonal_pairs), n_distinct(all_ibd$cluster)))

# Per-sample table: connected components of the clonal-pair graph, one
# component id per (cluster, component). Samples not in any clonal pair
# still get a row but with clonal_group = NA.
sample_cluster <- bind_rows(
  all_ibd %>% select(sample = sample1, cluster),
  all_ibd %>% select(sample = sample2, cluster)
) %>% distinct()

component_map <- clonal_pairs %>%
  group_by(cluster) %>%
  group_map(function(pairs_df, keys) {
    if (nrow(pairs_df) == 0) return(NULL)
    g <- graph_from_data_frame(pairs_df %>% select(sample1, sample2),
                               directed = FALSE)
    memb <- components(g)$membership
    tibble(
      sample       = names(memb),
      cluster      = keys$cluster,
      clonal_group = paste0(keys$cluster, "_", unname(memb))
    )
  }) %>% bind_rows()

# Zero clonal pairs across every cluster → empty tibble. Give it the
# expected columns so the left_join below doesn't try to look up
# non-existent join keys.
if (nrow(component_map) == 0) {
  component_map <- tibble(
    sample       = character(0),
    cluster      = character(0),
    clonal_group = character(0),
  )
}

per_sample <- sample_cluster %>%
  left_join(component_map, by = c("sample", "cluster")) %>%
  left_join(meta_keep, by = c("sample" = "sample_id"))

write_tsv(per_sample, out_clusters)
message(sprintf("[clonal_clusters] wrote %s (%d rows)",
                out_clusters, nrow(per_sample)))

# Focal-cluster pair file (V1 shape) if requested. Undoubled IDs already
# applied above; the join to metadata happens per side.
if (focal_pairs_out != "NULL" && focal_cluster != "NULL") {
  focal <- clonal_pairs %>%
    filter(cluster == focal_cluster) %>%
    select(sample1, sample2, fract_sites_IBD)

  # V1 emitted (Country/State/District/EnrolDate) per side. Agnostic emits
  # whatever roles are present, prefixed 1/2.
  side_join <- function(df, side) {
    m <- meta_keep %>%
      rename_with(~ paste0(.x, side),
                  .cols = setdiff(names(meta_keep), "sample_id"))
    df %>%
      left_join(m, by = setNames("sample_id", paste0("sample", side)))
  }
  focal_full <- focal %>% side_join("1") %>% side_join("2")
  write_tsv(focal_full, focal_pairs_out)
  message(sprintf("[clonal_clusters] focal=%s: wrote %s (%d pairs)",
                  focal_cluster, focal_pairs_out, nrow(focal_full)))
}
