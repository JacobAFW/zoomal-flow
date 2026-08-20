#!/usr/bin/env Rscript
# introgression_rule_sweep.R — which detection rule / threshold is best?
# (spec §8.1; brief §6)
# --------------------------------------------------------------------------
# Runs the detection call across a grid of rules and thresholds on a fixture
# with KNOWN ground truth, and tabulates detection rate against false
# positives. This is the empirical answer to "which threshold should we use" —
# it REPORTS, it does not change any default. The shipped default stays
# `absolute` @ 5e-4 (V1) until a human decides otherwise.
#
# Deliberately NOT a pipeline rule: it needs a ground-truth table, which only
# exists for synthetic fixtures. Run it by hand against the tiny cohort:
#
#   Rscript scripts/R/introgression_rule_sweep.R \
#     --genotype-table tests/tiny_cohort/outputs/ibd/combined/hmmIBD_input.tsv \
#     --clusters       tests/tiny_cohort/outputs/structure/admix_clusters.tsv \
#     --contig-map     tests/tiny_cohort/outputs/setup/contig_map.tsv \
#     --truth          tests/tiny_cohort/data/introgression_truth.tsv \
#     --pair           RegionA1__RegionB1 \
#     --window-size    10000 \
#     --min-snps       5 \
#     --out            tests/tiny_cohort/outputs/introgression/detection_rule_sweep.tsv
#
# Scoring is on the RAW per-pair calls, before the cross-dataset filters, so
# what you see is the detection rule's own behaviour rather than the filters'.
# The truth table is (sample_id, contig, window_bin) — every other
# (sample, window) that the rule calls is a false positive.
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
required <- c("genotype-table", "clusters", "contig-map", "truth", "pair",
              "window-size", "min-snps", "out")
missing  <- setdiff(required, names(args))
if (length(missing) > 0) stop("Missing args: ", paste(missing, collapse = ", "))

window_size <- as.integer(args[["window-size"]])
min_snps    <- as.integer(args[["min-snps"]])
out_tsv     <- args[["out"]]
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

parts <- strsplit(args[["pair"]], "__", fixed = TRUE)[[1]]
kx <- parts[1]; ky <- parts[2]

# -- ground truth: (sample, window) pairs that SHOULD be called -------------
cmap <- read_tsv(args[["contig-map"]], col_names = c("CONTIG", "CHROM"),
                 show_col_types = FALSE)
truth <- read_tsv(args[["truth"]], show_col_types = FALSE) %>%
  left_join(cmap, by = c("contig" = "CONTIG"))
if (any(is.na(truth$CHROM))) {
  stop("truth table contigs not present in contig_map.tsv: ",
       paste(unique(truth$contig[is.na(truth$CHROM)]), collapse = ", "))
}
truth_keys <- truth %>%
  transmute(key = paste0(sample_id, "@", window_id(CHROM, window_bin))) %>%
  pull(key) %>% unique()
message(sprintf("[sweep] ground truth: %d (sample, window) pairs", length(truth_keys)))

# -- distances + contours, computed once ------------------------------------
clusters_in <- read_tsv(args[["clusters"]], show_col_types = FALSE) %>%
  dplyr::select(SAMPLE = Sample, Cluster) %>% distinct() %>%
  filter(Cluster %in% c(kx, ky))

gt_wide <- fread(args[["genotype-table"]], header = TRUE, sep = "\t")
setnames(gt_wide, 1:2, c("CHROM", "POS"))
want <- intersect(clusters_in$SAMPLE, names(gt_wide))
gt_long <- melt(gt_wide[, c("CHROM", "POS", want), with = FALSE],
                id.vars = c("CHROM", "POS"), variable.name = "SAMPLE",
                value.name = "SNP", variable.factor = FALSE)
gt_long[, SNP := as.integer(SNP)]
gt_long <- as_tibble(gt_long) %>% left_join(clusters_in, by = "SAMPLE")

consensus <- dominant_alleles(gt_long, c(kx, ky))
dist_tbl  <- window_distances(gt_long, consensus, kx, ky, window_size, min_snps)
dens      <- density_contours(dist_tbl, kx, ky)
cont <- list()
cont[[kx]] <- dens$contours %>% filter(cluster == kx)
cont[[ky]] <- dens$contours %>% filter(cluster == ky)
pts <- tibble(x = dist_tbl$DIST_X, y = dist_tbl$DIST_Y)
message(sprintf("[sweep] %d sample-windows; contour levels: %s=%s | %s=%s",
                nrow(dist_tbl),
                kx, paste(sort(unique(cont[[kx]]$level)), collapse = ","),
                ky, paste(sort(unique(cont[[ky]]$level)), collapse = ",")))

# The membership levels are rule-independent, so compute them once and let
# each rule read off the same two numbers. This is exactly what
# is_introgressed() does internally; doing it here keeps the sweep cheap.
lvl <- list()
lvl[[kx]] <- contour_max_level(pts, cont[[kx]])
lvl[[ky]] <- contour_max_level(pts, cont[[ky]])

own_lvl   <- ifelse(dist_tbl$Cluster == kx, lvl[[kx]], lvl[[ky]])
other_lvl <- ifelse(dist_tbl$Cluster == kx, lvl[[ky]], lvl[[kx]])
keys      <- paste0(dist_tbl$SAMPLE, "@", dist_tbl$WINDOW)

score <- function(hit, rule, thr_other, thr_own) {
  called <- unique(keys[hit])
  tp <- sum(truth_keys %in% called)
  tibble(
    detection_rule      = rule,
    contour_level_other = thr_other,
    contour_level_own   = thr_own,
    n_truth             = length(truth_keys),
    n_calls             = length(called),
    n_true_positive     = tp,
    n_false_negative    = length(truth_keys) - tp,
    n_false_positive    = length(setdiff(called, truth_keys)),
    detection_rate      = if (length(truth_keys) > 0) tp / length(truth_keys) else NA_real_,
    precision           = if (length(called) > 0) tp / length(called) else NA_real_,
    n_windows_called    = n_distinct(dist_tbl$WINDOW[hit]),
    n_fp_windows        = n_distinct(setdiff(dist_tbl$WINDOW[hit],
                                             sub("^.*@", "", truth_keys)))
  )
}

# -- the grid --------------------------------------------------------------
LEVELS_OTHER <- c(1e-4, 2e-4, 5e-4, 1e-3, 2e-3)
LEVELS_OWN   <- c(2e-4, 5e-4, 1e-3)

rows <- list()
for (lo in LEVELS_OTHER) {
  for (lw in LEVELS_OWN) {
    hit <- (other_lvl >= lo) & !(own_lvl >= lw)
    rows[[length(rows) + 1]] <- score(hit, "absolute", lo, lw)
  }
}
# `relative` has no thresholds — one row.
rows[[length(rows) + 1]] <- score(other_lvl > own_lvl, "relative", NA, NA)

sweep <- bind_rows(rows) %>%
  arrange(detection_rule, contour_level_other, contour_level_own)
write_tsv(sweep, out_tsv)

cat("\n=== detection-rule / threshold sweep ===\n")
cat(sprintf("pair=%s  fixture truth=%d (sample, window) pairs  sample-windows=%d\n\n",
            args[["pair"]], length(truth_keys), nrow(dist_tbl)))
sweep %>%
  mutate(across(c(detection_rate, precision), ~ round(.x, 3))) %>%
  dplyr::select(detection_rule, contour_level_other, contour_level_own,
                n_true_positive, n_false_negative, n_false_positive,
                detection_rate, precision, n_fp_windows) %>%
  print(n = 100)
cat(sprintf("\n[sweep] wrote %s\n", out_tsv))
cat("[sweep] REPORT ONLY — the shipped default (absolute @ 5e-4, = V1) is unchanged.\n")
