#!/usr/bin/env Rscript
# declonalization_comparison.R — full vs unique, one stage at a time
# (docs/clonality.md)
# --------------------------------------------------------------------------
# WHAT: for one stage, read the SAME key numbers out of the `full` and the
#       `unique` output namespaces and put them side by side, with the delta.
#       One row per metric.
#
# WHY:  this is the artifact that makes "both outputs" worth producing. A
#       reader who does not want to blind-trust the pipeline can see, per
#       claim, how much of it survives collapsing clonal replicates. A result
#       that barely moves is a result clonality did not drive; one that halves
#       was substantially pseudo-replication. Neither is automatically wrong —
#       but the difference should never be invisible.
#
# READ IT LIKE THIS. `n_samples_removed` at the top of every table is the size
# of the de-clonalization. Then, per metric: `full`, `unique`, `delta`
# (unique − full) and `pct_change`. Counts that fall are the usual direction —
# fewer independent genotypes support fewer windows, fewer candidate regions.
# A metric that RISES is worth a look: it usually means a threshold that
# scales with cluster size (introgression's per-cluster support floor is a
# fraction/floor over cluster n) got easier to clear once the clonal padding
# was removed.
#
# DELIBERATELY DUMB. This script computes nothing scientific — it counts rows
# and reads numbers that the stages already wrote. Every statistic keeps its
# own null, its own FDR and its own script; putting the comparison in one
# place means the comparison logic cannot drift between stages.
#
# CLI:
#   Rscript scripts/R/declonalization_comparison.R \
#     --stage      structure|introgression|moi|selection \
#     --full-dir   outputs/structure/full \
#     --unique-dir outputs/structure/unique \
#     --audit      outputs/clonality/declonalization_audit.tsv \
#     --full-keep  outputs/clonality/all_genotypes.txt \
#     --unique-keep outputs/clonality/unique_genotypes.txt \
#     --out        outputs/structure/declonalization_comparison.tsv
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
required <- c("stage", "full-dir", "unique-dir", "out")
missing  <- setdiff(required, names(args))
if (length(missing) > 0) stop("Missing args: ", paste(missing, collapse = ", "))

STAGE  <- args[["stage"]]
FULL   <- args[["full-dir"]]
UNIQ   <- args[["unique-dir"]]
OUT    <- args[["out"]]
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
say <- function(...) message(sprintf("[declonalization_comparison] %s", sprintf(...)))

# Reading helpers that tolerate an absent file: a stage arm may legitimately
# not exist (a role not configured, a model that produced nothing). A missing
# input becomes NA, never a crash — the comparison is a report, not a gate.
read_or_null <- function(path) {
  if (!file.exists(path)) { say("absent: %s", path); return(NULL) }
  read_tsv(path, show_col_types = FALSE, progress = FALSE)
}
n_rows <- function(path) { d <- read_or_null(path); if (is.null(d)) NA_integer_ else nrow(d) }
n_distinct_col <- function(path, col) {
  d <- read_or_null(path)
  if (is.null(d) || !(col %in% names(d))) return(NA_integer_)
  dplyr::n_distinct(d[[col]])
}
count_lines <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NA_integer_)
  length(readLines(path, warn = FALSE))
}

metrics <- tibble(metric = character(), scope = character(),
                  full = numeric(), unique = numeric())
add <- function(metric, full, uniq, scope = "cohort") {
  metrics <<- bind_rows(metrics, tibble(metric = metric, scope = scope,
                                        full = as.numeric(full),
                                        unique = as.numeric(uniq)))
}

# --------------------------------------------------------------------------
# The size of the de-clonalization — the same first rows in every stage table
# --------------------------------------------------------------------------
n_full   <- count_lines(args[["full-keep"]])
n_unique <- count_lines(args[["unique-keep"]])
add("n_samples_analysed", n_full, n_unique)
if (!is.na(n_full) && !is.na(n_unique)) {
  add("n_samples_removed", 0, n_full - n_unique)
}
audit <- if (!is.null(args[["audit"]])) read_or_null(args[["audit"]]) else NULL
if (!is.null(audit)) {
  add("n_clonal_groups_collapsed", 0, nrow(audit))
}

# --------------------------------------------------------------------------
# Per-stage metrics
# --------------------------------------------------------------------------
if (STAGE == "structure") {
  # Best K, cluster count, and the per-cluster sizes. K changing is the
  # loudest possible signal that clonal padding was holding up a component.
  bk <- function(dir) {
    p <- file.path(dir, "best_k.txt")
    if (!file.exists(p)) return(NA_real_)
    as.numeric(trimws(readLines(p, warn = FALSE)[1]))
  }
  add("best_k", bk(FULL), bk(UNIQ))

  clu <- function(dir) read_or_null(file.path(dir, "admix_clusters.tsv"))
  cf <- clu(FULL); cu <- clu(UNIQ)
  add("n_clusters",
      if (is.null(cf)) NA else dplyr::n_distinct(cf$Cluster),
      if (is.null(cu)) NA else dplyr::n_distinct(cu$Cluster))
  add("n_samples_clustered",
      if (is.null(cf)) NA else dplyr::n_distinct(cf$Sample),
      if (is.null(cu)) NA else dplyr::n_distinct(cu$Sample))

  # Per-cluster size, joined on the cluster LABEL. Labels are derived from the
  # same role map in both arms, so they are comparable; a label present in one
  # arm and not the other shows up as NA and is itself the finding.
  size <- function(d) if (is.null(d)) tibble(Cluster = character(), n = integer())
                      else d %>% distinct(Sample, Cluster) %>% count(Cluster)
  sizes <- full_join(size(cf) %>% rename(full = n),
                     size(cu) %>% rename(unique = n), by = "Cluster")
  for (i in seq_len(nrow(sizes))) {
    add(sprintf("cluster_n[%s]", sizes$Cluster[i]), sizes$full[i], sizes$unique[i],
        scope = "cluster")
  }

  # PC1/PC2 variance — how much of the structure signal the clones were
  # carrying. A large shift means the leading axes were partly clonal.
  pcv <- function(dir) {
    d <- read_or_null(file.path(dir, "pca_variance.tsv"))
    if (is.null(d)) return(c(NA_real_, NA_real_))
    c(d$variance_percent[1], d$variance_percent[2])
  }
  vf <- pcv(FULL); vu <- pcv(UNIQ)
  add("pc1_variance_percent", vf[1], vu[1])
  add("pc2_variance_percent", vf[2], vu[2])

} else if (STAGE == "introgression") {
  add("n_introgressed_windows",
      n_distinct_col(file.path(FULL, "introgressed_windows_filtered.tsv"), "WINDOW"),
      n_distinct_col(file.path(UNIQ, "introgressed_windows_filtered.tsv"), "WINDOW"))
  add("n_window_calls",
      n_rows(file.path(FULL, "introgressed_windows_filtered.tsv")),
      n_rows(file.path(UNIQ, "introgressed_windows_filtered.tsv")))

  # Per-cluster window counts, from the stage's own summary table.
  wbc <- function(dir) {
    d <- read_or_null(file.path(dir, "windows_by_cluster.tsv"))
    if (is.null(d)) return(tibble(Cluster = character(), n = integer()))
    d %>% distinct(Cluster, WINDOW) %>% count(Cluster)
  }
  w <- full_join(wbc(FULL) %>% rename(full = n),
                 wbc(UNIQ) %>% rename(unique = n), by = "Cluster")
  for (i in seq_len(nrow(w))) {
    add(sprintf("windows[%s]", w$Cluster[i]), w$full[i], w$unique[i], scope = "cluster")
  }

  # The focal enrichment result, when a focal group is configured. Files are
  # named after the group, so glob rather than hardcode.
  enr <- function(dir) {
    f <- list.files(dir, pattern = "^focal_.*_enriched_windows\\.tsv$", full.names = TRUE)
    if (length(f) == 0) return(NA_integer_)
    sum(vapply(f, n_rows, integer(1)), na.rm = TRUE)
  }
  uw <- function(dir) {
    f <- list.files(dir, pattern = "^unique_windows_in_.*\\.tsv$", full.names = TRUE)
    if (length(f) == 0) return(NA_integer_)
    sum(vapply(f, n_rows, integer(1)), na.rm = TRUE)
  }
  add("n_focal_enriched_windows", enr(FULL), enr(UNIQ))
  add("n_focal_unique_windows_descriptive", uw(FULL), uw(UNIQ))

} else if (STAGE == "moi") {
  # Polyclonality rate per group, for whichever role tables exist. This is the
  # metric clonality bites hardest: a clonal block from one province counts
  # that province's polyclonality once per replicate.
  for (role in c("country", "geography")) {
    f <- read_or_null(file.path(FULL, sprintf("fws_by_%s.tsv", role)))
    u <- read_or_null(file.path(UNIQ, sprintf("fws_by_%s.tsv", role)))
    if (is.null(f) && is.null(u)) next
    key <- role
    pct_col <- intersect(c("pct_polyclonal", "percent_polyclonal", "prop_polyclonal"),
                         union(names(f), names(u)))
    n_col   <- intersect(c("n", "n_samples", "total"), union(names(f), names(u)))
    if (length(pct_col) == 0) {
      say("no polyclonality column found in fws_by_%s.tsv — columns: %s",
          role, paste(union(names(f), names(u)), collapse = ", "))
      next
    }
    pct_col <- pct_col[1]
    join <- full_join(
      if (is.null(f)) tibble(!!key := character()) else f %>% dplyr::select(all_of(c(key, pct_col, n_col))),
      if (is.null(u)) tibble(!!key := character()) else u %>% dplyr::select(all_of(c(key, pct_col, n_col))),
      by = key, suffix = c("_full", "_unique"))
    for (i in seq_len(nrow(join))) {
      add(sprintf("pct_polyclonal[%s=%s]", role, join[[key]][i]),
          join[[paste0(pct_col, "_full")]][i],
          join[[paste0(pct_col, "_unique")]][i], scope = role)
      if (length(n_col) > 0) {
        add(sprintf("n[%s=%s]", role, join[[key]][i]),
            join[[paste0(n_col[1], "_full")]][i],
            join[[paste0(n_col[1], "_unique")]][i], scope = role)
      }
    }
  }

} else if (STAGE == "selection") {
  # One row per model: candidate regions surviving the iHS thresholds.
  models <- union(list.dirs(FULL, recursive = FALSE, full.names = FALSE),
                  list.dirs(UNIQ, recursive = FALSE, full.names = FALSE))
  if (length(models) == 0) say("no selection models found under %s / %s", FULL, UNIQ)
  for (m in models) {
    add(sprintf("n_candidate_regions[%s]", m),
        n_rows(file.path(FULL, m, "candidate_regions_iHS.tsv")),
        n_rows(file.path(UNIQ, m, "candidate_regions_iHS.tsv")), scope = "model")
    add(sprintf("n_ihs_markers[%s]", m),
        n_rows(file.path(FULL, m, "ihs_table.tsv")),
        n_rows(file.path(UNIQ, m, "ihs_table.tsv")), scope = "model")
  }

} else {
  stop("Unknown --stage: ", STAGE,
       " (expected structure | introgression | moi | selection)")
}

# --------------------------------------------------------------------------
# Write + report
# --------------------------------------------------------------------------
out <- metrics %>%
  mutate(stage      = STAGE,
         delta      = unique - full,
         pct_change = if_else(!is.na(full) & full != 0,
                              100 * (unique - full) / full, NA_real_)) %>%
  relocate(stage) %>%
  dplyr::select(stage, metric, scope, full, unique, delta, pct_change)
write_tsv(out, OUT)

cat(sprintf("\n========= de-clonalization comparison: %s =========\n", STAGE))
cat("full = every sample | unique = one representative per clonal group\n")
cat("A metric that barely moves was not driven by clonality; one that moves a\n")
cat("lot was substantially pseudo-replication. Rising counts usually mean a\n")
cat("size-scaled threshold got easier to clear once the padding was removed.\n\n")
print(as.data.frame(out %>%
        mutate(across(c(full, unique, delta), ~ round(.x, 4)),
               pct_change = round(pct_change, 1))))
cat(sprintf("\nwrote %s\n", OUT))
cat("====================================================\n")
