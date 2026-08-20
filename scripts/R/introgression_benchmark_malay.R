#!/usr/bin/env Rscript
# introgression_benchmark_malay.R — score every detection rule against the
# published Malaysian result (spec §"Method status"; Stage 5b)
# --------------------------------------------------------------------------
# WHAT: runs the full Stage-5 introgression detection + cross-dataset filter
#       chain on the Malay benchmark fixture ONCE PER DETECTION RULE, scores
#       each rule against the published truth, and writes a scoreboard.
# WHY:  the two density-contour rules (`absolute`, `relative`) both read a
#       fitted 2D kernel surface, and that surface moves with cohort
#       composition: the bug-fixed density method reproduces only 41 of its own
#       217 published windows, and introgression rate scales INVERSELY with
#       cluster size. The deterministic `distance` rule was added to remove
#       that dependence. This script is the evidence for choosing between them.
#       It REPORTS; it changes no default.
#
# Four measures, in the order they should be weighted (see
# data/benchmark_malay/README.md, "How to read this benchmark"):
#
#   1. SAMPLE-LEVEL CONCORDANCE (primary). Which samples show introgression is
#      the stable, trustworthy signal — the diagnosis found it near-exact and
#      robust regardless of the positional-indexing bug. Precision / recall /
#      Jaccard against the 404 Mf + 117 Mn truth samples.
#      Read the `call_everything` baseline row before reading any of these: on
#      this fixture 521 of 558 samples are introgressed in truth, so a rule
#      that calls every sample already scores 0.93 precision. Sample-level
#      separates a broken rule from a working one; it does not rank good ones.
#
#   2. WINDOW OVERLAP vs TRUTH (secondary, caveated). Compared by CHROM +
#      window START, never by WINDOW id — the diagnosis flags WINDOW ids as
#      non-comparable between runs. The truth windows come from the SAME
#      density method and may embed the same bug, so this measures
#      "reproduces the published result", not "is correct". The
#      `regen_density` row is the published density method's own ceiling.
#
#   3. CLUSTER-SIZE TEST (the discriminator). Median windows per sample
#      against cluster n, per rule. The density rules over-call small clusters
#      because a small cluster's low, broad cloud lets its own members fall
#      below their own contour and satisfy the "not in its own cloud" half of
#      the rule. A rule with no fitted surface should be flat. Reported as
#      Spearman rho over the clusters plus the per-cluster medians.
#      Medians are ZERO-FILLED over every cluster member, not just the samples
#      that got a call — otherwise a rule that calls few samples looks flat by
#      construction.
#
#   4. PERTURBATION STABILITY. Drop a random 10% and 20% of samples, re-run,
#      and take the Jaccard of the window set (by CHROM + START) against the
#      full run. This is the direct measure of the failure mode: a rule whose
#      answer depends on who else is in the cohort scores low here.
#
# The fixture is REAL, EMBARGOED data under data/benchmark_malay/ (gitignored).
# Outputs go to a gitignored path. Neither is ever committed.
#
# Fidelity: this script does not reimplement anything. Distances come from
# window_distances() (introgression_core.R), the contours from
# density_contours(), the call from pair_calls() (introgression_detect.R), and
# the filters by SHELLING OUT to introgression_aggregate.R — the same code the
# pipeline runs. What it adds is the loop over rules, the perturbation draws,
# and the scoring.
#
# Cost note: the expensive step is the per-(sample, window) distance table, and
# it is rule-independent — so it is built ONCE per (scenario, pair) and every
# rule is called against it. That is 3 scenarios x C(K,2) pairs of heavy work,
# not 3 x C(K,2) x n_rules.
#
# CLI:
#   Rscript scripts/R/introgression_benchmark_malay.R \
#     --inputs-dir data/benchmark_malay/inputs \
#     --truth-dir  data/benchmark_malay/truth \
#     --regen-dir  data/benchmark_malay/regen_density \
#     --out-dir    outputs/benchmark_malay
#
# Everything else has a default that matches the fixture's published
# parameters (window 10 kb, per-window n > 5, per-cluster n > 5, SICA/KIR mask,
# hypervariable-window filter). See --help-style listing in DEFAULTS below.
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

DEFAULTS <- list(
  "inputs-dir"               = "data/benchmark_malay/inputs",
  "truth-dir"                = "data/benchmark_malay/truth",
  "regen-dir"                = "data/benchmark_malay/regen_density",
  "out-dir"                  = "outputs/benchmark_malay",
  # Fixture parameters (data/benchmark_malay/README.md, "Parameters").
  "window-size"              = "10000",
  "min-snps"                 = "5",     # strict >: a window needs >= 6 SNPs
  "min-samples-per-window"   = "5",     # published Malay chain used n > 5
  "per-cluster-min-pct"      = "0",     # ... and an ABSOLUTE per-cluster n > 5,
  "per-cluster-min-samples"  = "6",     #     which the % alone cannot express
  "gene-family-filters"      = "SICA,KIR",
  # Detection parameters.
  "variants"                 = "absolute,relative,distance,distance_adaptive",
  "contour-level-other"      = "5e-4",
  "contour-level-own"        = "5e-4",
  "distance-margin"          = "15",
  "distance-adaptive-quantile" = "0.9",
  # Perturbation stability.
  "perturb-fractions"        = "0.1,0.2",
  "seed"                     = "20260820"
)
for (k in names(DEFAULTS)) if (is.null(args[[k]])) args[[k]] <- DEFAULTS[[k]]

INPUTS   <- args[["inputs-dir"]]
TRUTH    <- args[["truth-dir"]]
REGEN    <- args[["regen-dir"]]
OUT      <- args[["out-dir"]]
WINDOW   <- as.integer(args[["window-size"]])
MIN_SNPS <- as.integer(args[["min-snps"]])
SEED     <- as.integer(args[["seed"]])
VARIANTS <- trimws(strsplit(args[["variants"]], ",")[[1]])
PERTURB  <- as.numeric(trimws(strsplit(args[["perturb-fractions"]], ",")[[1]]))
PERTURB  <- PERTURB[!is.na(PERTURB) & PERTURB > 0]

AGGREGATE_SCRIPT <- args[["aggregate-script"]] %||%
  file.path(SCRIPT_DIR, "introgression_aggregate.R")
RSCRIPT <- file.path(R.home("bin"), "Rscript")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

say <- function(...) message(sprintf("[benchmark_malay] %s", sprintf(...)))
tick <- function() proc.time()[["elapsed"]]

# --------------------------------------------------------------------------
# Variant table — a variant is a named bundle of detection params.
# --------------------------------------------------------------------------
variant_params <- function(v) {
  base <- list(
    contour_level_other = as.numeric(args[["contour-level-other"]]),
    contour_level_own   = as.numeric(args[["contour-level-own"]]),
    distance_margin     = as.numeric(args[["distance-margin"]]),
    distance_adaptive_quantile = as.numeric(args[["distance-adaptive-quantile"]])
  )
  switch(v,
    absolute          = c(base, list(detection_rule = "absolute",  distance_adaptive = FALSE)),
    relative          = c(base, list(detection_rule = "relative",  distance_adaptive = FALSE)),
    distance          = c(base, list(detection_rule = "distance",  distance_adaptive = FALSE)),
    distance_adaptive = c(base, list(detection_rule = "distance",  distance_adaptive = TRUE)),
    stop("unknown variant: ", v, " (expected absolute | relative | distance | distance_adaptive)")
  )
}
VARIANT_PARAMS <- setNames(lapply(VARIANTS, variant_params), VARIANTS)
NEEDS_DENSITY  <- any(vapply(VARIANT_PARAMS,
                             function(p) rule_needs_density(p$detection_rule), logical(1)))

# --------------------------------------------------------------------------
# Fixture: turn the raw benchmark inputs into the four files the pipeline's
# own scripts expect. Nothing here is Malay-specific except the tess3r
# component order, which the fixture README pins.
# --------------------------------------------------------------------------
say("building fixture tables from %s", INPUTS)

fam <- read.table(file.path(INPUTS, "cleaned.fam"), header = FALSE,
                  stringsAsFactors = FALSE)
qmat <- as.matrix(read.table(file.path(INPUTS, "cleaned.3.Q"), header = FALSE))
if (nrow(fam) != nrow(qmat)) {
  stop("cleaned.fam (", nrow(fam), ") and cleaned.3.Q (", nrow(qmat),
       ") disagree on sample count")
}
# Component order is fixed by the fixture README: c("Mn","Mf","Peninsular").
Q_LABELS <- c("Mn", "Mf", "Peninsular")
clusters_all <- tibble(
  Sample  = as.character(fam[[1]]),
  Cluster = Q_LABELS[apply(qmat, 1, which.max)]
)
say("clusters: %s",
    paste(sprintf("%s n=%d", names(table(clusters_all$Cluster)),
                  as.integer(table(clusters_all$Cluster))), collapse = ", "))

# Pseudo-.fai from the coordinate bed (name, length) — the aggregate step only
# reads columns 1 and 2, and the fixture ships a bed rather than a .fai.
bed <- read.table(file.path(INPUTS, "strain_A1_H.1.Icor.fasta.bed"),
                  header = FALSE, stringsAsFactors = FALSE)
fai_tbl <- tibble(CONTIG = bed[[1]], LEN = as.numeric(bed[[3]]))
# Contig codes mirror scripts/py/contigs_from_fai.py: file order, MIT/API
# dropped, 1-based. That reproduces the "01".."14" codes in hmmIBD.tsv.
nuclear <- fai_tbl %>% filter(!grepl("MIT|API", CONTIG))
cmap_tbl <- nuclear %>% mutate(CHROM = row_number()) %>% dplyr::select(CONTIG, CHROM)

FIXTURE <- file.path(OUT, "_fixture")
dir.create(FIXTURE, recursive = TRUE, showWarnings = FALSE)
write.table(fai_tbl, file.path(FIXTURE, "ref.fasta.fai"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(cmap_tbl, file.path(FIXTURE, "contig_map.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
# Minimal metadata: sample_id only. No `geography` role, so the aggregate
# step's per-geography summary skips itself with a logged note (as designed).
write_tsv(tibble(sample_id = clusters_all$Sample),
          file.path(FIXTURE, "samples.tsv"))

GFF <- file.path(INPUTS, "strain_A1_H.1.Icor.gff3")
if (!file.exists(GFF)) stop("GFF not found: ", GFF)

# --------------------------------------------------------------------------
# Genotype table — read once, kept wide. 62k variants x ~750 samples.
# --------------------------------------------------------------------------
t0 <- tick()
gt_all <- fread(file.path(INPUTS, "hmmIBD.tsv"), header = TRUE, sep = "\t")
setnames(gt_all, 1:2, c("CHROM", "POS"))
gt_all[, CHROM := as.integer(CHROM)]
gt_all[, POS := as.integer(POS)]
in_gt <- intersect(clusters_all$Sample, names(gt_all))
say("genotype table: %d variants x %d columns; %d of %d fixture samples present (%.1fs)",
    nrow(gt_all), ncol(gt_all), length(in_gt), nrow(clusters_all), tick() - t0)
clusters_all <- clusters_all %>% filter(Sample %in% in_gt)
gt_all <- gt_all[, c("CHROM", "POS", clusters_all$Sample), with = FALSE]
stopifnot(setequal(unique(gt_all$CHROM), cmap_tbl$CHROM))

# --------------------------------------------------------------------------
# One scenario = one sample set. Runs every pair, every variant, and the
# aggregate filter chain; returns the per-variant filtered call table.
# --------------------------------------------------------------------------
run_scenario <- function(label, keep_samples) {
  clus <- clusters_all %>% filter(Sample %in% keep_samples)
  sizes <- clus %>% count(Cluster, name = "n")
  say("scenario '%s': %d samples (%s)", label, nrow(clus),
      paste(sprintf("%s=%d", sizes$Cluster, sizes$n), collapse = ", "))

  sdir <- if (label == "full") OUT else file.path(OUT, "_perturb", label)
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  clus_file <- file.path(sdir, "_clusters.tsv")
  write_tsv(clus, clus_file)

  pairs <- combn(sort(unique(clus$Cluster)), 2, simplify = FALSE)

  for (pr in pairs) {
    kx <- pr[1]; ky <- pr[2]
    pair_id <- paste0(kx, "__", ky)
    t1 <- tick()
    members <- clus %>% filter(Cluster %in% c(kx, ky))

    gt_long <- melt(gt_all[, c("CHROM", "POS", members$Sample), with = FALSE],
                    id.vars = c("CHROM", "POS"), variable.name = "SAMPLE",
                    value.name = "SNP", variable.factor = FALSE)
    gt_long[, SNP := as.integer(SNP)]
    gt_long <- gt_long[SNP >= 0]                 # drop missing calls up front
    gt_long <- as_tibble(gt_long) %>%
      left_join(members %>% dplyr::select(SAMPLE = Sample, Cluster), by = "SAMPLE")

    consensus <- dominant_alleles(gt_long, c(kx, ky))
    dist_tbl  <- window_distances(gt_long, consensus, kx, ky,
                                  window_size_bp = WINDOW, min_snps = MIN_SNPS)
    rm(gt_long, consensus); invisible(gc(FALSE))
    say("  %s / %s: %d sample-windows (%.0fs)", label, pair_id,
        nrow(dist_tbl), tick() - t1)

    contours <- setNames(vector("list", 2), c(kx, ky))
    if (NEEDS_DENSITY && nrow(dist_tbl) >= 10) {
      dens <- density_contours(dist_tbl, kx, ky)
      contours[[kx]] <- dens$contours %>% filter(cluster == kx)
      contours[[ky]] <- dens$contours %>% filter(cluster == ky)
      rm(dens); invisible(gc(FALSE))
    }

    for (v in VARIANTS) {
      pdir <- file.path(sdir, v, "pairs")
      dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
      calls <- if (nrow(dist_tbl) < 10) NULL else
        pair_calls(dist_tbl, kx, ky, contours, VARIANT_PARAMS[[v]], pair_id = pair_id)
      if (is.null(calls) || nrow(calls) == 0) {
        calls <- as_tibble(setNames(
          replicate(length(PAIR_CALL_COLS), character(0), simplify = FALSE),
          PAIR_CALL_COLS))
      } else {
        calls <- calls %>% arrange(CHROM, BIN, Cluster, SAMPLE)
      }
      write_tsv(calls, file.path(pdir, paste0(pair_id, ".tsv")))
      say("    %-18s %d calls", v, nrow(calls))
    }
    rm(dist_tbl, contours); invisible(gc(FALSE))
  }

  # Filters: the pipeline's own aggregate script, unmodified.
  out <- list()
  for (v in VARIANTS) {
    vdir <- file.path(sdir, v)
    pair_files <- list.files(file.path(vdir, "pairs"), pattern = "\\.tsv$",
                             full.names = TRUE)
    logf <- file.path(vdir, "aggregate.log")
    st <- system2(RSCRIPT, c(
      shQuote(AGGREGATE_SCRIPT),
      "--clusters",                shQuote(clus_file),
      "--metadata",                shQuote(file.path(FIXTURE, "samples.tsv")),
      "--contig-map",              shQuote(file.path(FIXTURE, "contig_map.tsv")),
      "--fai",                     shQuote(file.path(FIXTURE, "ref.fasta.fai")),
      "--gff",                     shQuote(GFF),
      "--gene-family-filters",     shQuote(args[["gene-family-filters"]]),
      "--window-size",             WINDOW,
      "--min-samples-per-window",  args[["min-samples-per-window"]],
      "--per-cluster-min-pct",     args[["per-cluster-min-pct"]],
      "--per-cluster-min-samples", args[["per-cluster-min-samples"]],
      "--out-dir",                 shQuote(vdir),
      shQuote(pair_files)
    ), stdout = logf, stderr = logf)
    if (st != 0) stop("introgression_aggregate.R failed for ", label, "/", v,
                      " — see ", logf)
    filtered <- read_tsv(file.path(vdir, "introgressed_windows_filtered.tsv"),
                         show_col_types = FALSE,
                         col_types = cols(WINDOW = col_character(),
                                          SAMPLE = col_character(),
                                          Cluster = col_character(),
                                          .default = col_guess()))
    say("  %s / %-18s filtered: %d calls, %d windows, %d samples",
        label, v, nrow(filtered),
        dplyr::n_distinct(filtered$WINDOW), dplyr::n_distinct(filtered$SAMPLE))
    out[[v]] <- filtered
  }
  list(filtered = out, clusters = clus)
}

# --------------------------------------------------------------------------
# Truth + reference sets
# --------------------------------------------------------------------------
read_window_set <- function(dir) {
  files <- list.files(dir, pattern = "_windows\\.tsv$", full.names = TRUE)
  if (length(files) == 0) stop("no *_windows.tsv under ", dir)
  files %>%
    map(~ read_tsv(.x, show_col_types = FALSE, col_types = cols(.default = col_guess()))) %>%
    bind_rows() %>%
    transmute(key = paste0(as.integer(CHROM), ":", as.integer(start))) %>%
    pull(key) %>% unique()
}
read_sample_set <- function(dir) {
  files <- list.files(dir, pattern = "_samples_window_counts\\.tsv$", full.names = TRUE)
  if (length(files) == 0) stop("no *_samples_window_counts.tsv under ", dir)
  files %>% map(~ read_tsv(.x, show_col_types = FALSE)) %>% bind_rows() %>%
    pull(SAMPLE) %>% unique()
}

truth_windows <- read_window_set(TRUTH)
truth_samples <- read_sample_set(TRUTH)
regen_windows <- read_window_set(REGEN)
say("truth: %d windows, %d samples | regen_density: %d windows",
    length(truth_windows), length(truth_samples), length(regen_windows))

# Window key from a filtered call table: CHROM + window START, NEVER the
# WINDOW id (the diagnosis flags WINDOW ids as non-comparable across runs).
window_keys <- function(filtered) {
  if (nrow(filtered) == 0) return(character(0))
  unique(paste0(as.integer(filtered$CHROM), ":", as.integer(filtered$START)))
}

set_scores <- function(called, reference, prefix) {
  inter <- length(intersect(called, reference))
  uni   <- length(union(called, reference))
  setNames(list(
    length(called), inter,
    if (length(called) > 0) inter / length(called) else NA_real_,
    if (length(reference) > 0) inter / length(reference) else NA_real_,
    if (uni > 0) inter / uni else NA_real_
  ), paste0(prefix, c("_n", "_shared", "_precision", "_recall", "_jaccard")))
}

# --------------------------------------------------------------------------
# Cluster-size test: median windows per sample vs cluster n, zero-filled over
# every cluster member.
# --------------------------------------------------------------------------
cluster_size_table <- function(filtered, clus, variant) {
  per_sample <- filtered %>% distinct(SAMPLE, WINDOW) %>%
    count(SAMPLE, name = "n_windows")
  clus %>%
    left_join(per_sample, by = c("Sample" = "SAMPLE")) %>%
    mutate(n_windows = replace_na(n_windows, 0)) %>%
    group_by(Cluster) %>%
    summarise(cluster_n = dplyr::n(),
              median_windows_per_sample = median(n_windows),
              mean_windows_per_sample   = mean(n_windows),
              pct_samples_called = 100 * mean(n_windows > 0),
              .groups = "drop") %>%
    mutate(variant = variant, .before = 1) %>%
    arrange(desc(cluster_n))
}

# Spearman rho over the clusters gives the DIRECTION of any size dependence.
# With K = 3 it has only four possible values (-1, -0.5, 0.5, 1), so a rank
# swap between two nearly-equal medians reads as rho = -0.5 — direction alone
# is not enough. The spread ratio (max median / min median) gives the
# MAGNITUDE, and the two must be read together: rho = -1 at ratio 1.05 is
# noise, rho = -0.5 at ratio 6 is a real effect.
size_bias_rho <- function(tbl) {
  if (nrow(tbl) < 3) return(NA_real_)
  suppressWarnings(cor(tbl$cluster_n, tbl$median_windows_per_sample,
                       method = "spearman"))
}
size_bias_ratio <- function(tbl) {
  m <- tbl$median_windows_per_sample
  if (length(m) < 2 || min(m) <= 0) return(NA_real_)
  max(m) / min(m)
}
size_bias_label <- function(rho, ratio) {
  if (is.na(rho)) return("undetermined")
  if (!is.na(ratio) && ratio < 1.25) return("flat (spread < 1.25x)")
  if (rho <= -0.5) "INVERSE (small clusters over-called)"
  else if (rho >= 0.5) "positive (large clusters over-called)"
  else "mixed"
}

# --------------------------------------------------------------------------
# Run: full cohort, then the perturbation draws.
# --------------------------------------------------------------------------
T_START <- tick()
full <- run_scenario("full", clusters_all$Sample)

set.seed(SEED)
perturb_results <- list()
for (frac in PERTURB) {
  # Drop uniformly at random across the whole cohort, so cluster proportions
  # move too — that is the cohort shift the stability test is about.
  n_drop <- floor(nrow(clusters_all) * frac)
  keep <- setdiff(clusters_all$Sample,
                  sample(clusters_all$Sample, n_drop))
  lbl <- sprintf("drop%02.0f", frac * 100)
  perturb_results[[lbl]] <- run_scenario(lbl, keep)
}

# --------------------------------------------------------------------------
# Scoreboard
# --------------------------------------------------------------------------
size_tbls <- list()
rows <- list()

for (v in VARIANTS) {
  filtered <- full$filtered[[v]]
  called_w <- window_keys(filtered)
  called_s <- unique(filtered$SAMPLE)

  st <- cluster_size_table(filtered, full$clusters, v)
  size_tbls[[v]] <- st
  rho   <- size_bias_rho(st)
  ratio <- size_bias_ratio(st)

  stab_labels <- names(perturb_results)
  stability <- if (length(stab_labels) == 0) list() else {
    as.list(setNames(vapply(stab_labels, function(lbl) {
      pw <- window_keys(perturb_results[[lbl]]$filtered[[v]])
      u  <- length(union(called_w, pw))
      if (u == 0) NA_real_ else length(intersect(called_w, pw)) / u
    }, numeric(1)), paste0("stability_jaccard_", stab_labels)))
  }

  row <- c(
    list(variant = v,
         detection_rule = VARIANT_PARAMS[[v]]$detection_rule,
         adaptive = isTRUE(VARIANT_PARAMS[[v]]$distance_adaptive),
         n_calls = nrow(filtered)),
    set_scores(called_s, truth_samples, "sample"),
    set_scores(called_w, truth_windows, "window"),
    set_scores(called_w, regen_windows, "vs_regen"),
    list(size_bias_rho = rho, size_bias_ratio = ratio,
         size_bias = size_bias_label(rho, ratio)),
    stability
  )
  rows[[v]] <- as_tibble(row)
}

# Reference rows, for calibration rather than competition:
#   regen_density  — the published density method's own output on this fixture
#                    (its 19% window ceiling), scored the same way.
#   call_everything— the degenerate "every sample, every window" rule, so the
#                    sample-level numbers can be read against their floor.
ref_rows <- list(
  as_tibble(c(
    list(variant = "regen_density (reference)", detection_rule = "density (published fix)",
         adaptive = NA, n_calls = NA_integer_),
    set_scores(character(0), truth_samples, "sample"),
    set_scores(regen_windows, truth_windows, "window"),
    set_scores(regen_windows, regen_windows, "vs_regen"),
    list(size_bias_rho = NA_real_, size_bias_ratio = NA_real_, size_bias = "n/a")
  )),
  as_tibble(c(
    list(variant = "call_everything (baseline)", detection_rule = "none",
         adaptive = NA, n_calls = NA_integer_),
    set_scores(clusters_all$Sample, truth_samples, "sample"),
    set_scores(character(0), truth_windows, "window"),
    set_scores(character(0), regen_windows, "vs_regen"),
    list(size_bias_rho = NA_real_, size_bias_ratio = NA_real_, size_bias = "n/a")
  ))
)

scoreboard <- bind_rows(bind_rows(rows), bind_rows(ref_rows))
write_tsv(scoreboard, file.path(OUT, "scoreboard.tsv"))
write_tsv(bind_rows(size_tbls), file.path(OUT, "cluster_size_bias.tsv"))

# Per-window detail, so a disagreement can be chased to a coordinate.
detail <- bind_rows(lapply(VARIANTS, function(v) {
  cw <- window_keys(full$filtered[[v]])
  tibble(variant = v, window_key = union(cw, truth_windows)) %>%
    mutate(called = window_key %in% cw,
           in_truth = window_key %in% truth_windows,
           in_regen_density = window_key %in% regen_windows)
}))
write_tsv(detail, file.path(OUT, "window_overlap_detail.tsv"))

# --------------------------------------------------------------------------
# Human summary
# --------------------------------------------------------------------------
cat("\n=================== Malay introgression benchmark ===================\n")
cat(sprintf("fixture: %s  |  %d samples, %d sites, %d clusters\n",
            INPUTS, nrow(clusters_all), nrow(gt_all), n_distinct(clusters_all$Cluster)))
cat(sprintf("filters: window %d bp, per-window n > %s, per-cluster n >= %s, mask [%s]\n",
            WINDOW, args[["min-samples-per-window"]],
            args[["per-cluster-min-samples"]], args[["gene-family-filters"]]))
cat(sprintf("truth: %d samples / %d windows   density ceiling: %d windows\n\n",
            length(truth_samples), length(truth_windows), length(regen_windows)))

cat("--- 1. SAMPLE-LEVEL concordance (PRIMARY) ---\n")
scoreboard %>%
  dplyr::select(variant, sample_n, sample_shared, sample_precision,
                sample_recall, sample_jaccard) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = 100)

cat("\n--- 2. WINDOW overlap vs published truth (SECONDARY, caveated) ---\n")
scoreboard %>%
  dplyr::select(variant, window_n, window_shared, window_recall, window_jaccard,
                vs_regen_shared, vs_regen_jaccard) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = 100)

cat("\n--- 3. CLUSTER-SIZE test (THE DISCRIMINATOR) ---\n")
bind_rows(size_tbls) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% print(n = 100)
cat("\n")
scoreboard %>% filter(!is.na(size_bias_rho)) %>%
  dplyr::select(variant, size_bias_rho, size_bias_ratio, size_bias) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = 100)
cat("  rho = direction (only 4 values at K=3); ratio = max/min median = magnitude.\n")
cat("  NOTE: the per-cluster filter is part of what this measures. A PERCENTAGE\n")
cat("  threshold (per_cluster_min_pct) scales with cluster n and is itself a\n")
cat("  size-dependent filter; an absolute floor (per_cluster_min_samples) is not.\n")
cat("  Re-run with --per-cluster-min-pct 0.05 --per-cluster-min-samples 0 to see\n")
cat("  how much of any apparent bias is the filter rather than the rule.\n")

if (length(perturb_results) > 0) {
  cat("\n--- 4. PERTURBATION stability (window-set Jaccard vs the full run) ---\n")
  scoreboard %>%
    dplyr::select(variant, starts_with("stability_jaccard_")) %>%
    filter(if_any(starts_with("stability_jaccard_"), ~ !is.na(.x))) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = 100)
}

cat(sprintf("\nwrote %s\n     %s\n     %s\n",
            file.path(OUT, "scoreboard.tsv"),
            file.path(OUT, "cluster_size_bias.tsv"),
            file.path(OUT, "window_overlap_detail.tsv")))
cat(sprintf("elapsed: %.1f min\n", (tick() - T_START) / 60))
cat("REPORT ONLY — the shipped default (absolute @ 5e-4) is unchanged.\n")
cat("=====================================================================\n")
