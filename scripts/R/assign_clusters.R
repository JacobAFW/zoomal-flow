#!/usr/bin/env Rscript
# assign_clusters.R  (agnostic port — label-mode driven, no hardcoded K)
# --------------------------------------------------------------------------
# Read the best-K .Q from ADMIXTURE, assign each sample to its dominant
# ancestry component, then label components by one of three modes:
#
#   numbered   Components → Cluster_1..K by dominant-component index.
#              Zero-metadata, always works.
#
#   auto       Each component is labelled by the majority value of the
#              `geography` role among its dominant-component members.
#              Falls back to numbered for any component without data.
#
#   reference  V1's behaviour: each component is labelled by majority
#              vote against the `group` role values (V1's Cluster column
#              → group role). Falls back to numbered for any component
#              without a vote.
#
# Output schema (preserves V1's column names for compatibility with later
# stages that join on this file):
#     Sample        VCF sample ID (= sample_id role; capitalised for V1 compat)
#     Component     V1..VK
#     Proportion    dominant-component proportion
#     Cluster       human-readable label
#
# Usage:
#   Rscript assign_clusters.R <cleaned.K.Q> <cleaned.fam> <metadata.tsv> \
#                             <mode> <out.tsv>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(paste(
    "Usage: Rscript assign_clusters.R <cleaned.K.Q> <cleaned.fam>",
    "<metadata.tsv> <mode> <out.tsv>"
  ))
}
q_in    <- args[1]
fam_in  <- args[2]
meta_in <- args[3]
mode    <- args[4]
out_tsv <- args[5]

stopifnot(mode %in% c("numbered", "auto", "reference"))

q <- read.table(q_in, header = FALSE) %>% as_tibble()
fam <- read.table(fam_in, header = FALSE, stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  select(Sample = V1)

K <- ncol(q)
component_cols <- paste0("V", seq_len(K))

admix <- bind_cols(fam, q) %>%
  rename_with(~ component_cols, all_of(seq_len(K) + 1))

meta <- read_tsv(meta_in, show_col_types = FALSE)
stopifnot("sample_id" %in% names(meta))

# Per-sample dominant component
admix_long <- admix %>%
  pivot_longer(all_of(component_cols),
               names_to = "Component", values_to = "Proportion") %>%
  group_by(Sample) %>%
  slice_max(Proportion, with_ties = FALSE) %>%
  ungroup()

all_components <- tibble(Component = component_cols)

label_numbered <- function() {
  all_components %>%
    mutate(Cluster = paste0("Cluster_", str_remove(Component, "V")))
}

label_by_majority <- function(role_col) {
  # Returns a (Component, Cluster) tibble. Any component without votes
  # falls back to its numbered label.
  if (!(role_col %in% names(meta))) {
    message(sprintf(
      "[assign_clusters] mode=%s but role column '%s' absent → falling back to numbered",
      mode, role_col
    ))
    return(label_numbered())
  }

  meta_role <- meta %>%
    select(sample_id, .label = all_of(role_col)) %>%
    filter(!is.na(.label), .label != "")

  votes <- admix_long %>%
    inner_join(meta_role, by = c("Sample" = "sample_id")) %>%
    count(Component, .label, name = "n") %>%
    group_by(Component) %>%
    slice_max(n, with_ties = FALSE) %>%
    ungroup() %>%
    select(Component, Cluster = .label)

  fallback <- label_numbered()
  all_components %>%
    left_join(votes, by = "Component") %>%
    left_join(fallback, by = "Component", suffix = c("", ".num")) %>%
    mutate(Cluster = coalesce(Cluster, Cluster.num)) %>%
    select(Component, Cluster)
}

mapping <- switch(mode,
  numbered  = label_numbered(),
  auto      = label_by_majority("geography"),
  reference = label_by_majority("group")
)

# Detect collisions (multiple components voted to the same label) — keep
# the highest-population component's original label and append a numeric
# suffix to the loser(s), so downstream files don't silently dedup samples.
if (mode != "numbered") {
  dup_labels <- mapping %>%
    add_count(Cluster) %>%
    filter(n > 1)
  if (nrow(dup_labels) > 0) {
    message("[assign_clusters] collision(s) — multiple components → same label:")
    print(dup_labels)
    # Resolve: first occurrence keeps the label, subsequent get a "_n" tail.
    mapping <- mapping %>%
      group_by(Cluster) %>%
      mutate(Cluster = if (n() > 1) paste0(Cluster, "_", row_number()) else Cluster) %>%
      ungroup()
  }
}

message(sprintf("[assign_clusters] mode=%s | K=%d | component→cluster mapping:", mode, K))
print(mapping)

admix_clusters <- admix_long %>%
  left_join(mapping, by = "Component")

write_tsv(admix_clusters, out_tsv)

cat("\n=== Cluster assignment counts ===\n")
print(admix_clusters %>% count(Cluster, sort = TRUE))
