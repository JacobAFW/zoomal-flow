#!/usr/bin/env Rscript
# introgression_aggregate.R — stitch the per-pair calls, apply the cross-
# dataset filters, write the summary tables (spec §3, §6)
# --------------------------------------------------------------------------
# WHAT: reads every per-pair call file, then applies V1's four cross-dataset
#       filters IN V1's ORDER (find_introgression.R:276-350):
#
#         1. dataset low-n     window must be called in more than
#                              `min_samples_per_window` samples (V1: > 2)
#         2. per-cluster low-n window must be called in at least
#                              ceiling(cluster_n * per_cluster_min_pct)
#                              of that cluster's samples (V1: 5%)
#         3. gene family       windows overlapping any GFF feature whose
#                              attribute matches a `gene_family_filters`
#                              keyword are masked (V1: literal SICAvar + KIR)
#         4. hypervariable     windows called in more than one cluster are
#                              dropped — shared across clusters means the
#                              window is variable, not introgressed
#
# WHY:  the per-pair scripts are deliberately local (one pair, one density
#       pass); every judgement that needs to see the WHOLE dataset — how many
#       samples share a window, whether a window is hypervariable across
#       clusters, whether it sits in a hypervariable gene family — belongs
#       here, once, after the fan-out.
#
# Generalisations from V1 (all documented in docs/introgression_analysis_spec.md §7):
#   - cluster names come from admix_clusters.tsv, never a literal list;
#   - the gene-family keywords are config, not a hardcoded SICA|KIR grep;
#   - GFF contig names are reconciled to the pipeline's integer CHROM codes
#     via contig_map.tsv + the reference .fai (V1 hardcoded a PlasmoDB
#     LT727… → ordered_PKNH_NN_v2 length match);
#   - the per-cluster percentage denominator is the cluster's TRUE size from
#     admix_clusters.tsv. V1 used "samples that had any call at all" as a
#     proxy, which is always <= cluster size and so slightly more permissive.
#   - window coordinates are the window's own bounds (BIN ± window_size/2)
#     rather than V1's observed min/max SNP position inside it, so the
#     coordinates don't move with marker density.
#
# CLI (flags first, then one or more per-pair call TSVs as bare arguments):
#   Rscript introgression_aggregate.R \
#     --clusters    outputs/structure/admix_clusters.tsv \
#     --metadata    outputs/metadata/samples.tsv \
#     --contig-map  outputs/setup/contig_map.tsv \
#     --window-size 10000 \
#     --min-samples-per-window 2 \
#     --per-cluster-min-pct    0.05 \
#     --gff         data/reference/annotation.gff   # or NULL
#     --fai         data/reference/ref.fasta.fai \
#     --gene-family-filters "SICA,KIR"              # or "" to disable
#     --out-dir     outputs/introgression \
#     outputs/introgression/pairs/*.tsv
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

script_args <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- dirname(sub("^--file=", "",
                          script_args[grep("^--file=", script_args)][1]))
source(file.path(SCRIPT_DIR, "introgression_core.R"))

parse_args <- function(av) {
  flags <- list(); pos <- character()
  i <- 1
  while (i <= length(av)) {
    if (startsWith(av[i], "--")) {
      flags[[sub("^--", "", av[i])]] <- av[i + 1]; i <- i + 2
    } else {
      pos <- c(pos, av[i]); i <- i + 1
    }
  }
  list(flags = flags, pos = pos)
}
parsed    <- parse_args(commandArgs(trailingOnly = TRUE))
args      <- parsed$flags
call_files<- parsed$pos

required <- c("clusters", "metadata", "contig-map", "window-size",
              "min-samples-per-window", "per-cluster-min-pct", "out-dir")
missing  <- setdiff(required, names(args))
if (length(missing) > 0) stop("Missing args: ", paste(missing, collapse = ", "))
if (length(call_files) == 0) stop("No per-pair call files given")

window_size  <- as.integer(args[["window-size"]])
min_samp_win <- as.integer(args[["min-samples-per-window"]])
pct_min      <- as.numeric(args[["per-cluster-min-pct"]])
out_dir      <- args[["out-dir"]]
gff_path     <- args[["gff"]]
fai_path     <- args[["fai"]]
gene_kw      <- args[["gene-family-filters"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

is_null_arg <- function(x) is.null(x) || is.na(x) || !nzchar(x) || x %in% c("NULL", "None", "null")

gene_keywords <- if (is_null_arg(gene_kw)) character(0) else {
  kw <- trimws(strsplit(gene_kw, ",", fixed = TRUE)[[1]])
  kw[nzchar(kw)]
}

# -- inputs ----------------------------------------------------------------
clusters_in <- read_tsv(args[["clusters"]], show_col_types = FALSE) %>%
  dplyr::select(SAMPLE = Sample, Cluster) %>%
  distinct()
cluster_sizes <- clusters_in %>% count(Cluster, name = "cluster_n")

meta <- read_tsv(args[["metadata"]], show_col_types = FALSE)
stopifnot("sample_id" %in% names(meta))

cmap <- read_tsv(args[["contig-map"]], col_names = c("CONTIG", "CHROM"),
                 show_col_types = FALSE)

calls <- call_files %>%
  map(~ read_tsv(.x, show_col_types = FALSE,
                 col_types = cols(WINDOW = col_character(),
                                  SAMPLE = col_character(),
                                  Cluster = col_character(),
                                  CHROM = col_integer(), BIN = col_double(),
                                  .default = col_guess()))) %>%
  bind_rows()

audit <- tibble(step = character(), n_calls = integer(),
                n_windows = integer(), n_samples = integer())
log_step <- function(df, label) {
  audit <<- bind_rows(audit, tibble(
    step = label, n_calls = nrow(df),
    n_windows = n_distinct(df$WINDOW), n_samples = n_distinct(df$SAMPLE)))
  message(sprintf("[introgression_aggregate] %-26s calls=%d windows=%d samples=%d",
                  label, nrow(df), n_distinct(df$WINDOW), n_distinct(df$SAMPLE)))
  df
}

calls <- log_step(calls, "raw (all pairs)")

# Pre-filter per-window sample counts: the input to the shoulder plot, which
# exists to help the user PICK `min_samples_per_window` — so it has to be
# drawn before filter 1 removes anything (legacy
# introgression_multi_cluster.R:311-322). Compact by design (one row per
# window), so the raw call table never has to be written out.
calls %>%
  distinct(WINDOW, CHROM, BIN, SAMPLE, Cluster) %>%
  group_by(WINDOW, CHROM, BIN) %>%
  summarise(n_samples = n_distinct(SAMPLE), n_clusters = n_distinct(Cluster),
            .groups = "drop") %>%
  arrange(n_samples) %>%
  write_tsv(file.path(out_dir, "window_sample_counts_raw.tsv"))

# Everything below tolerates an empty input: a cohort with no introgression is
# a result, not a crash.
if (nrow(calls) == 0) {
  message("[introgression_aggregate] no calls from any pair — writing empty outputs")
}

# -- filter 1: dataset-level low-n ----------------------------------------
keep_windows <- calls %>%
  distinct(WINDOW, SAMPLE) %>%
  count(WINDOW, name = "n_samples") %>%
  filter(n_samples > min_samp_win) %>%
  pull(WINDOW)
calls <- calls %>% filter(WINDOW %in% keep_windows) %>%
  log_step(sprintf("1. dataset n > %d", min_samp_win))

# -- filter 2: per-cluster minimum percentage ------------------------------
per_cluster_keep <- calls %>%
  distinct(Cluster, WINDOW, SAMPLE) %>%
  count(Cluster, WINDOW, name = "n_samples") %>%
  left_join(cluster_sizes, by = "Cluster") %>%
  mutate(threshold = ceiling(cluster_n * pct_min)) %>%
  filter(n_samples >= threshold) %>%
  dplyr::select(Cluster, WINDOW)
calls <- calls %>% semi_join(per_cluster_keep, by = c("Cluster", "WINDOW")) %>%
  log_step(sprintf("2. per-cluster >= %g%%", pct_min * 100))

# -- filter 3: hypervariable gene families from the GFF --------------------
# Reconcile the GFF's contig names with the pipeline's integer CHROM codes.
# Three strategies, most trustworthy first:
#   (a) exact name match against contig_map.tsv (the easy case: the GFF and the
#       reference FASTA use the same contig names);
#   (b) exact sequence-length match via the reference .fai — V1's approach,
#       needed because PlasmoDB ships `LT727…` accessions where this project's
#       reference uses `ordered_PKNH_NN_v2`;
#   (c) MUTUAL nearest length within `tol`, which (a) and (b) both miss when the
#       two assemblies differ by a handful of corrected bases.
#
# (c) matters. On the Indo cohort's PlasmoDB-68 GFF vs the A1.H.1 Icor
# reference, only 3 of 14 chromosomes have byte-identical lengths (the Icor
# corrections shift the rest by 99-7138 bp), so exact-length matching alone
# silently masked gene families on 3 chromosomes and let them through on the
# other 11 — a bug V1 shares. `tol` is a fraction of length; the mutual-best
# requirement stops a near-collision (e.g. two chromosomes 0.3% apart) from
# producing a wrong pairing, and every inexact pairing is logged and written to
# gff_contig_map.tsv for audit.
GFF_LENGTH_TOL <- 0.01

#' Map GFF contig names -> pipeline CHROM codes. Returns tibble(GFF_NAME, CHROM,
#' method) plus (as an attribute) the GFF contigs that could not be mapped.
map_gff_contigs <- function(gff_lengths, fai_lengths, tol = GFF_LENGTH_TOL) {
  mapped <- tibble(GFF_NAME = character(), CHROM = integer(), method = character())
  # (a) exact name
  by_name <- gff_lengths %>%
    inner_join(fai_lengths, by = c("GFF_NAME" = "CONTIG"),
               suffix = c("_gff", "_fai")) %>%
    transmute(GFF_NAME, CHROM, method = "name")
  mapped <- bind_rows(mapped, by_name)

  remaining_gff <- gff_lengths %>% filter(!(GFF_NAME %in% mapped$GFF_NAME))
  remaining_fai <- fai_lengths %>% filter(!(CHROM %in% mapped$CHROM))

  # (b) exact length
  by_len <- remaining_gff %>%
    inner_join(remaining_fai, by = c("LEN" = "LEN")) %>%
    transmute(GFF_NAME, CHROM, method = "length")
  # A length shared by two contigs is ambiguous — drop those.
  by_len <- by_len %>% group_by(GFF_NAME) %>% filter(dplyr::n() == 1) %>%
    group_by(CHROM) %>% filter(dplyr::n() == 1) %>% ungroup()
  mapped <- bind_rows(mapped, by_len)

  remaining_gff <- gff_lengths %>% filter(!(GFF_NAME %in% mapped$GFF_NAME))
  remaining_fai <- fai_lengths %>% filter(!(CHROM %in% mapped$CHROM))

  # (c) mutual nearest length within tol
  if (nrow(remaining_gff) > 0 && nrow(remaining_fai) > 0) {
    grid <- tidyr::expand_grid(
        remaining_gff %>% dplyr::select(GFF_NAME, LEN_GFF = LEN),
        remaining_fai %>% dplyr::select(CHROM, LEN_FAI = LEN)
      ) %>%
      filter(!is.na(LEN_GFF), !is.na(LEN_FAI)) %>%
      mutate(rel = abs(LEN_GFF - LEN_FAI) / LEN_FAI) %>%
      filter(rel <= tol)
    best_gff <- grid %>% group_by(GFF_NAME) %>% slice_min(rel, n = 1, with_ties = FALSE) %>%
      ungroup() %>% dplyr::select(GFF_NAME, CHROM, rel)
    best_fai <- grid %>% group_by(CHROM) %>% slice_min(rel, n = 1, with_ties = FALSE) %>%
      ungroup() %>% dplyr::select(GFF_NAME, CHROM)
    mutual <- best_gff %>% semi_join(best_fai, by = c("GFF_NAME", "CHROM"))
    if (nrow(mutual) > 0) {
      for (i in seq_len(nrow(mutual))) {
        message(sprintf(
          "[introgression_aggregate] GFF contig %s -> CHROM %s by nearest length (%.3f%% apart)",
          mutual$GFF_NAME[i], mutual$CHROM[i], 100 * mutual$rel[i]))
      }
      mapped <- bind_rows(mapped, mutual %>%
                            transmute(GFF_NAME, CHROM, method = "nearest_length"))
    }
  }

  unmapped <- setdiff(gff_lengths$GFF_NAME, mapped$GFF_NAME)
  attr(mapped, "unmapped") <- unmapped
  mapped
}

masked_windows <- character(0)
if (is_null_arg(gff_path)) {
  message("[introgression_aggregate] filter 3 skipped: no `introgression.gff` configured")
} else if (length(gene_keywords) == 0) {
  message("[introgression_aggregate] filter 3 skipped: `gene_family_filters` is empty")
} else if (!file.exists(gff_path)) {
  message("[introgression_aggregate] filter 3 skipped: GFF not found at ", gff_path)
} else {
  gff_lines <- readLines(gff_path, warn = FALSE)
  seqreg <- grep("^##sequence-region", gff_lines, value = TRUE)

  fai_lengths <- if (!is_null_arg(fai_path) && file.exists(fai_path)) {
    read_tsv(fai_path, col_names = FALSE, show_col_types = FALSE) %>%
      transmute(CONTIG = X1, LEN = as.numeric(X2)) %>%
      inner_join(cmap, by = "CONTIG")
  } else {
    cmap %>% mutate(LEN = NA_real_)
  }

  gff_lengths <- if (length(seqreg) > 0) {
    tibble(line = seqreg) %>%
      separate(line, into = c("tag", "GFF_NAME", "start", "LEN"),
               sep = "\\s+", fill = "right", extra = "drop") %>%
      transmute(GFF_NAME, LEN = suppressWarnings(as.numeric(LEN)))
  } else {
    # No ##sequence-region header: fall back to the contig names in column 1.
    tibble(GFF_NAME = unique(sub("\t.*$", "", grep("^[^#]", gff_lines, value = TRUE))),
           LEN = NA_real_)
  }

  gff_name_to_chrom <- map_gff_contigs(gff_lengths, fai_lengths)
  unmapped <- attr(gff_name_to_chrom, "unmapped")
  message(sprintf("[introgression_aggregate] GFF contig mapping: %d of %d mapped (%s)",
                  nrow(gff_name_to_chrom), nrow(gff_lengths),
                  paste(sprintf("%s=%d", names(table(gff_name_to_chrom$method)),
                                as.integer(table(gff_name_to_chrom$method))),
                        collapse = ", ")))
  if (length(unmapped) > 0) {
    message(sprintf(paste0("[introgression_aggregate] %d GFF contig(s) unmapped ",
                           "(expected for contigs excluded from the nuclear set): %s"),
                    length(unmapped), paste(unmapped, collapse = ", ")))
  }
  write_tsv(gff_name_to_chrom %>% left_join(cmap, by = "CHROM"),
            file.path(out_dir, "gff_contig_map.tsv"))

  gff <- read.table(gff_path, sep = "\t", quote = "", comment.char = "#",
                    stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    setNames(c("GFF_NAME", "source", "feature", "start", "end",
               "score", "strand", "frame", "attribute")) %>%
    inner_join(gff_name_to_chrom %>% dplyr::select(GFF_NAME, CHROM), by = "GFF_NAME")

  pattern <- paste(gene_keywords, collapse = "|")
  hits <- gff %>% filter(grepl(pattern, attribute))
  message(sprintf("[introgression_aggregate] gene-family keywords [%s]: %d GFF features on %d of %d mapped contigs",
                  paste(gene_keywords, collapse = ", "), nrow(hits),
                  n_distinct(hits$CHROM), n_distinct(gff$CHROM)))

  if (nrow(hits) > 0) {
    masked_windows <- hits %>%
      mutate(start_bin = window_bin(start, window_size),
             end_bin   = window_bin(end,   window_size)) %>%
      dplyr::select(CHROM, start_bin, end_bin) %>%
      pmap(function(CHROM, start_bin, end_bin) {
        window_id(CHROM, seq(start_bin, end_bin, by = window_size))
      }) %>%
      unlist() %>% unique()
    write_tsv(tibble(WINDOW = masked_windows),
              file.path(out_dir, "gene_family_masked_windows.tsv"))
  }
  calls <- calls %>% filter(!(WINDOW %in% masked_windows)) %>%
    log_step(sprintf("3. gene-family mask (%d win)", length(masked_windows)))
}

# -- filter 4: hypervariable across clusters -------------------------------
multi_cluster <- calls %>%
  distinct(WINDOW, Cluster) %>%
  count(WINDOW, name = "n_clusters") %>%
  filter(n_clusters > 1) %>%
  pull(WINDOW)
calls <- calls %>% filter(!(WINDOW %in% multi_cluster)) %>%
  log_step(sprintf("4. hypervariable (%d win)", length(multi_cluster)))

# -- window coordinates ----------------------------------------------------
with_coords <- calls %>%
  mutate(START = as.integer(BIN - window_size %/% 2),
         END   = as.integer(BIN - window_size %/% 2 + window_size - 1)) %>%
  left_join(cmap, by = "CHROM") %>%
  relocate(CONTIG, .after = CHROM)

write_tsv(with_coords, file.path(out_dir, "introgressed_windows_filtered.tsv"))
write_tsv(audit,       file.path(out_dir, "filter_audit.tsv"))

# -- summaries -------------------------------------------------------------
# Per-cluster: distribution of per-sample window counts (V1's
# average_windows_for_clusters.tsv).
with_coords %>%
  distinct(SAMPLE, Cluster, WINDOW) %>%
  count(SAMPLE, Cluster, name = "n") %>%
  group_by(Cluster) %>%
  summarise(mean = mean(n), median = median(n), sd = sd(n),
            n_samples = dplyr::n(), .groups = "drop") %>%
  arrange(desc(median)) %>%
  write_tsv(file.path(out_dir, "average_windows_for_clusters.tsv"))

# Per-sample, with V1's tertile High/Medium/Low banding.
with_coords %>%
  distinct(SAMPLE, Cluster, WINDOW) %>%
  count(SAMPLE, Cluster, name = "n") %>%
  arrange(desc(n)) %>%
  mutate(intro_level = ntile(n, 3),
         intro_level = case_when(intro_level == 3 ~ "High",
                                 intro_level == 2 ~ "Medium",
                                 TRUE             ~ "Low")) %>%
  write_tsv(file.path(out_dir, "intro_per_sample_summary.tsv"))

# Per-window, per-cluster, with direction breakdown + coordinates.
with_coords %>%
  group_by(Cluster, WINDOW, CHROM, CONTIG, BIN, START, END) %>%
  summarise(n_samples  = n_distinct(SAMPLE),
            directions = paste(sort(unique(DIRECTION)), collapse = ";"),
            .groups = "drop") %>%
  arrange(desc(n_samples)) %>%
  write_tsv(file.path(out_dir, "windows_by_cluster.tsv"))

# Windows per chromosome per cluster.
with_coords %>%
  distinct(Cluster, CHROM, CONTIG, WINDOW) %>%
  count(Cluster, CHROM, CONTIG, name = "n_windows") %>%
  arrange(Cluster, CHROM) %>%
  write_tsv(file.path(out_dir, "windows_across_chrom.tsv"))

# Per-geography summary — role-driven, skipped (with a note) when the
# geography role isn't configured.
if ("geography" %in% names(meta)) {
  with_coords %>%
    left_join(meta %>% dplyr::select(SAMPLE = sample_id, geography),
              by = "SAMPLE") %>%
    distinct(SAMPLE, geography, WINDOW) %>%
    count(SAMPLE, geography, name = "n") %>%
    group_by(geography) %>%
    summarise(mean = mean(n), median = median(n), sd = sd(n),
              n_samples = dplyr::n(), .groups = "drop") %>%
    arrange(desc(median)) %>%
    write_tsv(file.path(out_dir, "introgression_by_geography.tsv"))
} else {
  message("[introgression_aggregate] geography role absent — per-geography summary skipped")
}

cat(sprintf(paste0("\n[introgression_aggregate] %d filtered calls | %d windows | ",
                   "%d samples | %d clusters\n"),
            nrow(with_coords), n_distinct(with_coords$WINDOW),
            n_distinct(with_coords$SAMPLE), n_distinct(with_coords$Cluster)))
