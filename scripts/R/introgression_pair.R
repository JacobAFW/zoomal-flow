#!/usr/bin/env Rscript
# introgression_pair.R — per-(Kx, Ky) introgression detection (spec §2, §3)
# --------------------------------------------------------------------------
# WHAT: one cluster pair per invocation. Subsets the combined genotype table
#       to the Kx ∪ Ky samples, computes each cluster's consensus alleles and
#       every sample-window's percent-mismatch distance to both, draws the 2D
#       kernel-density cloud per cluster in (distance-to-Kx, distance-to-Ky)
#       space, and calls introgression in BOTH directions with the pluggable
#       detection rule.
# WHY:  V1's implementation was intrinsically three-cluster — two distance
#       axes with Mn/Mf/Peninsular welded in by name — so it could not scale
#       as ADMIXTURE K climbs. A clean 2-cluster 2D plot per pair keeps the
#       exact density-contour method, removes every cluster-name hardcode, and
#       works at any K. V1's three-way "ambiguous between Mf and Mn" case is
#       recovered downstream as "flagged in both the Kx-Mf and Kx-Mn pairs".
#
# The maths lives in scripts/R/introgression_core.R (distance + contour cores,
# preserved from V1) and scripts/R/introgression_detect.R (the yes/no call).
# This script is wiring: read, subset, call, write.
#
# GUARDRAIL (spec §4): all-pairs is C(K,2) and grows quadratically. The
# warning lives in the Snakefile (it knows the pair count before any pair
# runs); each pair here is a full density pass, so read that warning before
# launching a large-K run.
#
# Parallelism: none in R. Snakemake fans the {pair} wildcard across cores.
#
# CLI:
#   Rscript introgression_pair.R \
#     --genotype-table outputs/ibd/combined/hmmIBD_input.tsv \
#     --clusters       outputs/structure/admix_clusters.tsv \
#     --pair           Mf__Mn \
#     --window-size    10000 \
#     --min-snps       5 \
#     --detection-rule absolute \
#     --contour-level-other 5e-4 \
#     --contour-level-own   5e-4 \
#     --distance-margin 15 \
#     --distance-adaptive false \
#     --distance-adaptive-quantile 0.9 \
#     --out            outputs/introgression/pairs/Mf__Mn.tsv
#
# Under `--detection-rule distance` the density surface is never built: that
# rule decides from the raw per-window distances, so a contour pass would be
# wasted work AND would reintroduce the cohort dependence the rule exists to
# avoid. The LEVEL_OWN / LEVEL_OTHER output columns are NA in that mode.
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

script_args <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- dirname(sub("^--file=", "",
                          script_args[grep("^--file=", script_args)][1]))
source(file.path(SCRIPT_DIR, "introgression_core.R"))
source(file.path(SCRIPT_DIR, "introgression_detect.R"))

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
required <- c("genotype-table", "clusters", "pair", "window-size", "min-snps",
              "detection-rule", "contour-level-other", "contour-level-own", "out")
missing  <- setdiff(required, names(args))
if (length(missing) > 0) stop("Missing args: ", paste(missing, collapse = ", "))

pair_id     <- args[["pair"]]
window_size <- as.integer(args[["window-size"]])
min_snps    <- as.integer(args[["min-snps"]])
out_tsv     <- args[["out"]]
# `distance_*` are optional so an existing caller that only passes the density
# args keeps working; the defaults mirror config/config.yaml.
as_logical_arg <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  isTRUE(tolower(as.character(x)) %in% c("true", "t", "yes", "1"))
}
det_params  <- list(
  detection_rule      = args[["detection-rule"]],
  contour_level_other = as.numeric(args[["contour-level-other"]]),
  contour_level_own   = as.numeric(args[["contour-level-own"]]),
  distance_margin     = as.numeric(args[["distance-margin"]] %||% 15),
  distance_adaptive   = as_logical_arg(args[["distance-adaptive"]]),
  distance_adaptive_quantile =
    as.numeric(args[["distance-adaptive-quantile"]] %||% 0.9)
)
if (!det_params$detection_rule %in% DETECTION_RULES) {
  stop("--detection-rule must be one of: ", paste(DETECTION_RULES, collapse = ", "),
       " (got: ", det_params$detection_rule, ")")
}

# The {pair} wildcard encodes (Kx, Ky) as "Kx__Ky" — a double underscore, so
# single-underscore cluster labels survive the round trip.
parts <- strsplit(pair_id, "__", fixed = TRUE)[[1]]
if (length(parts) != 2) {
  stop("--pair must be '<Kx>__<Ky>' (double underscore separator); got: ", pair_id)
}
kx <- parts[1]; ky <- parts[2]

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

# Empty-but-valid output, so a degenerate pair never breaks the DAG.
OUT_COLS <- c("PAIR", "SAMPLE", "Cluster", "OTHER", "DIRECTION", "CHROM", "BIN",
              "WINDOW", "N_SNPS", "DIST_OWN", "DIST_OTHER",
              "LEVEL_OWN", "LEVEL_OTHER")
write_empty <- function(reason) {
  message("[introgression_pair] ", pair_id, ": ", reason,
          " — writing header-only output")
  empty <- as_tibble(setNames(
    replicate(length(OUT_COLS), character(0), simplify = FALSE), OUT_COLS))
  write_tsv(empty, out_tsv)
  quit(status = 0)
}

# -- cluster membership ----------------------------------------------------
clusters_in <- read_tsv(args[["clusters"]], show_col_types = FALSE) %>%
  dplyr::select(SAMPLE = Sample, Cluster) %>%
  distinct()
pair_members <- clusters_in %>% filter(Cluster %in% c(kx, ky))
n_x <- sum(pair_members$Cluster == kx)
n_y <- sum(pair_members$Cluster == ky)
message(sprintf("[introgression_pair] %s: %s n=%d, %s n=%d",
                pair_id, kx, n_x, ky, n_y))
if (n_x == 0 || n_y == 0) write_empty("one side of the pair has no samples")

# -- combined genotype table ------------------------------------------------
# Format (from scripts/R/genotype_table.R): chrom<TAB>pos<TAB>sample columns,
# integer CHROM (contig_map-coded upstream), hmmIBD encoding (-1 = missing),
# sample IDs already undoubled.
gt_wide <- fread(args[["genotype-table"]], header = TRUE, sep = "\t")
chrom_col <- names(gt_wide)[1]
pos_col   <- names(gt_wide)[2]
setnames(gt_wide, c(chrom_col, pos_col), c("CHROM", "POS"))

want <- intersect(pair_members$SAMPLE, names(gt_wide))
absent <- setdiff(pair_members$SAMPLE, names(gt_wide))
if (length(absent) > 0) {
  message(sprintf("[introgression_pair] %d pair sample(s) absent from the genotype table: %s",
                  length(absent), paste(head(absent, 10), collapse = ", ")))
}
if (length(want) < 2) write_empty("fewer than 2 pair samples in the genotype table")

gt_wide <- gt_wide[, c("CHROM", "POS", want), with = FALSE]

gt_long <- melt(gt_wide, id.vars = c("CHROM", "POS"),
                variable.name = "SAMPLE", value.name = "SNP",
                variable.factor = FALSE)
gt_long[, SNP := as.integer(SNP)]
gt_long <- as_tibble(gt_long) %>%
  left_join(pair_members, by = "SAMPLE")
message(sprintf("[introgression_pair] %s: %d variants x %d samples",
                pair_id, nrow(gt_wide), length(want)))

# -- core 1: consensus alleles + per-(sample, window) distances -------------
consensus <- dominant_alleles(gt_long, c(kx, ky))
message(sprintf("[introgression_pair] %s: %d positions with a tie-free consensus in both clusters",
                pair_id, nrow(consensus)))

dist_tbl <- window_distances(gt_long, consensus, kx, ky,
                             window_size_bp = window_size, min_snps = min_snps)
message(sprintf("[introgression_pair] %s: %d sample-windows (>%d SNPs, non-tied distances)",
                pair_id, nrow(dist_tbl), min_snps))
if (nrow(dist_tbl) < 10) write_empty("too few sample-windows to call on")

# -- core 2: 2D density contours per cluster -------------------------------
# Skipped entirely under `distance`: that rule reads the raw distances, so a
# contour pass would be wasted work and would put the cohort-dependent surface
# back into a decision built to avoid it.
contours <- setNames(vector("list", 2), c(kx, ky))
if (rule_needs_density(det_params$detection_rule)) {
  dens <- density_contours(dist_tbl, kx, ky)
  if (nrow(dens$points) != nrow(dist_tbl)) {
    stop("[introgression_pair] ggplot point layer lost rows (",
         nrow(dens$points), " vs ", nrow(dist_tbl), ") — cannot align calls")
  }
  lvl_counts <- dens$contours %>% count(cluster, name = "n_vertices")
  message("[introgression_pair] contour vertices: ",
          paste(sprintf("%s=%d", lvl_counts$cluster, lvl_counts$n_vertices),
                collapse = ", "))

  contours[[kx]] <- dens$contours %>% filter(cluster == kx)
  contours[[ky]] <- dens$contours %>% filter(cluster == ky)
  if (nrow(contours[[kx]]) == 0 || nrow(contours[[ky]]) == 0) {
    write_empty("one cluster produced no density contour (too few / too tied points)")
  }

  # The point layer uses the identity stat on continuous scales, so its built
  # x/y ARE the distance columns. Test against dist_tbl directly rather than
  # relying on ggplot_build preserving row order; the row-count and value
  # checks above confirm the built layer is the same point set V1 tested.
  stopifnot(isTRUE(all.equal(sort(dens$points$x), sort(dist_tbl$DIST_X))))
} else {
  message("[introgression_pair] rule '", det_params$detection_rule,
          "' needs no density surface — contour step skipped")
}

# -- detection, both directions -------------------------------------------
# Direction 1: a Kx sample whose window looks like Ky ("Kx_like_Ky").
# Direction 2: a Ky sample whose window looks like Kx.
# pair_calls() (introgression_detect.R) owns the orientation, the adaptive
# rule's per-pair calibration, and the output schema, so the pipeline and the
# Malay benchmark harness call introgression identically.
calls <- pair_calls(dist_tbl, kx, ky, contours, det_params, pair_id = pair_id)

if (is.null(calls) || nrow(calls) == 0) {
  write_empty("no introgressed sample-windows called")
}

calls <- calls %>% arrange(CHROM, BIN, Cluster, SAMPLE)
write_tsv(calls, out_tsv)

rule_desc <- if (det_params$detection_rule == "distance") {
  if (isTRUE(det_params$distance_adaptive)) {
    sprintf("rule=distance, adaptive q=%g", det_params$distance_adaptive_quantile)
  } else {
    sprintf("rule=distance, margin>=%g pp", det_params$distance_margin)
  }
} else {
  sprintf("rule=%s, other>=%g, own>=%g", det_params$detection_rule,
          det_params$contour_level_other, det_params$contour_level_own)
}
cat(sprintf(paste0("\n[introgression_pair] %s (%s): ",
                   "%d sample-window calls across %d windows / %d samples\n"),
            pair_id, rule_desc,
            nrow(calls), n_distinct(calls$WINDOW), n_distinct(calls$SAMPLE)))
calls %>% count(DIRECTION, name = "n_calls") %>%
  mutate(line = sprintf("  %-30s %d", DIRECTION, n_calls)) %>%
  pull(line) %>% cat(sep = "\n")
cat("\n")
