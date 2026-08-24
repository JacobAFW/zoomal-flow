#!/usr/bin/env Rscript
# introgression_floor_sweep.R — how the Stage-5 result moves with the
# per-cluster support floor (spec §9.5)
# --------------------------------------------------------------------------
# WHAT: re-runs the Stage-5 aggregate filter chain (and, if a focal group is
#       configured, the headline) once per candidate
#       `introgression.per_cluster_min_samples` value, and tabulates how many
#       windows survive — overall, per cluster, and in the headline.
# WHY:  a derived floor is only trustworthy if the answer is not perched on a
#       cliff beside it. `introgression_floor_derivation.R` says which N is
#       false-discovery-safe; this says what that N costs, and whether N-1 or
#       N+1 would have said something materially different. Run them together.
#
# Detection is NOT re-run: the cached per-pair call files are the input, and
# only the cross-dataset filter chain is repeated. That makes a 20-point sweep
# minutes rather than hours.
#
# Cohort-agnostic: cluster names come from the clusters table, the focal group
# from --focal-group (omit it and the headline columns are simply absent), so
# the same script sweeps the real cohort and any benchmark fixture.
#
# CLI:
#   Rscript scripts/R/introgression_floor_sweep.R \
#     --cohort      indo \
#     --clusters    outputs/structure/admix_clusters.tsv \
#     --metadata    outputs/metadata/samples.tsv \
#     --contig-map  outputs/setup/contig_map.tsv \
#     --fai         data/reference/strain_A1_H.1.Icor.fasta.fai \
#     --gff         data/reference/PlasmoDB_version/PlasmoDB-68_PknowlesiA1H1.gff \
#     --gene-family-filters "SICA,KIR" \
#     --window-size 10000 \
#     --min-samples-per-window 2 \
#     --focal-group Aceh \
#     --n-values    "1,2,3,4,5,6,7,8,9,10,12,16,20,25,30,31,32,33,34,35" \
#     --work-dir    outputs/introgression/_floor_sweep \
#     --out         outputs/introgression/floor_stability_sweep.tsv \
#     outputs/introgression/pairs/*.tsv
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

script_args <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- dirname(sub("^--file=", "",
                          script_args[grep("^--file=", script_args)][1]))

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
required <- c("clusters", "metadata", "contig-map", "out", "work-dir")
missing  <- setdiff(required, names(args))
if (length(missing) > 0) stop("Missing args: ", paste(missing, collapse = ", "))
if (length(call_files) == 0) stop("No per-pair call files given")

COHORT   <- args[["cohort"]] %||% "cohort"
WINDOW   <- as.integer(args[["window-size"]] %||% 10000)
MIN_SAMP <- as.integer(args[["min-samples-per-window"]] %||% 2)
FOCAL    <- args[["focal-group"]]
WORK     <- args[["work-dir"]]
OUT      <- args[["out"]]
N_VALUES <- as.integer(trimws(strsplit(args[["n-values"]] %||%
                                         "1,2,3,4,5,6,7,8,9,10", ",")[[1]]))
RSCRIPT  <- file.path(R.home("bin"), "Rscript")
AGG      <- file.path(SCRIPT_DIR, "introgression_aggregate.R")
HEAD     <- file.path(SCRIPT_DIR, "introgression_headline.R")

dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
say <- function(...) message(sprintf("[floor_sweep] %s", sprintf(...)))

read_calls <- function(f) read_tsv(f, show_col_types = FALSE,
  col_types = cols(WINDOW = col_character(), SAMPLE = col_character(),
                   Cluster = col_character(), .default = col_guess()))

rows <- list()
for (N in N_VALUES) {
  vdir <- file.path(WORK, sprintf("N%02d", N))
  dir.create(vdir, recursive = TRUE, showWarnings = FALSE)
  logf <- file.path(vdir, "aggregate.log")

  st <- system2(RSCRIPT, c(
    shQuote(AGG),
    "--clusters",                shQuote(args[["clusters"]]),
    "--metadata",                shQuote(args[["metadata"]]),
    "--contig-map",              shQuote(args[["contig-map"]]),
    if (!is.null(args[["fai"]])) c("--fai", shQuote(args[["fai"]])) else NULL,
    if (!is.null(args[["gff"]])) c("--gff", shQuote(args[["gff"]])) else NULL,
    "--gene-family-filters",     shQuote(args[["gene-family-filters"]] %||% "NULL"),
    "--window-size",             WINDOW,
    "--min-samples-per-window",  MIN_SAMP,
    # The floor IS the sweep variable, so the percentage is switched off
    # entirely — otherwise the two thresholds interact and the sweep would not
    # be measuring one parameter.
    "--per-cluster-min-pct",     0,
    "--per-cluster-min-samples", N,
    "--out-dir",                 shQuote(vdir),
    shQuote(call_files)
  ), stdout = logf, stderr = logf)
  if (st != 0) stop("introgression_aggregate.R failed at N = ", N, " — see ", logf)

  filtered <- read_calls(file.path(vdir, "introgressed_windows_filtered.tsv"))
  rows[[length(rows) + 1]] <- tibble(
    cohort = COHORT, N = N, scope = "all",
    n_windows = n_distinct(filtered$WINDOW),
    n_calls   = nrow(filtered),
    n_samples = n_distinct(filtered$SAMPLE))
  if (nrow(filtered) > 0) {
    rows[[length(rows) + 1]] <- filtered %>%
      group_by(scope = Cluster) %>%
      summarise(n_windows = n_distinct(WINDOW), n_calls = dplyr::n(),
                n_samples = n_distinct(SAMPLE), .groups = "drop") %>%
      mutate(cohort = COHORT, N = N, .before = 1)
  }

  headline_n <- NA_integer_
  if (!is.null(FOCAL) && nzchar(FOCAL) && !FOCAL %in% c("NULL", "null")) {
    hlog <- file.path(vdir, "headline.log")
    uniq <- file.path(vdir, "headline_unique.tsv")
    sth <- system2(RSCRIPT, c(
      shQuote(HEAD),
      "--calls",         shQuote(file.path(vdir, "introgressed_windows_filtered.tsv")),
      "--clusters",      shQuote(args[["clusters"]]),
      "--metadata",      shQuote(args[["metadata"]]),
      "--focal-group",   shQuote(FOCAL),
      "--out-dir",       shQuote(vdir),
      "--out-unique",    shQuote(uniq),
      "--out-per-chrom", shQuote(file.path(vdir, "headline_per_chrom.tsv"))
    ), stdout = hlog, stderr = hlog)
    if (sth != 0) stop("introgression_headline.R failed at N = ", N, " — see ", hlog)
    hl <- read_tsv(uniq, show_col_types = FALSE)
    headline_n <- nrow(hl)
    rows[[length(rows) + 1]] <- tibble(
      cohort = COHORT, N = N, scope = sprintf("headline:%s", FOCAL),
      n_windows = headline_n, n_calls = NA_integer_,
      n_samples = if (nrow(hl) > 0 && "n_samples" %in% names(hl)) max(hl$n_samples) else 0L)
  }
  say("N=%2d: %4d windows / %5d calls / %3d samples%s", N,
      n_distinct(filtered$WINDOW), nrow(filtered), n_distinct(filtered$SAMPLE),
      if (is.na(headline_n)) "" else sprintf("  | headline %d windows", headline_n))
}

sweep <- bind_rows(rows) %>% arrange(cohort, N, scope)
write_tsv(sweep, OUT)

cat(sprintf("\n=========== per-cluster floor stability sweep (%s) ===========\n", COHORT))
sweep %>%
  dplyr::select(N, scope, n_windows) %>%
  pivot_wider(names_from = scope, values_from = n_windows) %>%
  as.data.frame() %>% print()
cat(sprintf("\nwrote %s\n", OUT))
cat("REPORT ONLY — this script changes no config.\n")
cat("=============================================================\n")
