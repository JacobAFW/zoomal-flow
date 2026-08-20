#!/usr/bin/env Rscript
# introgression_core.R — the two preserved computational cores (spec §2.2)
# --------------------------------------------------------------------------
# NOT a standalone script: sourced by introgression_pair.R and by the
# rule/threshold sweep. It holds the two pieces of maths the agnostic revamp
# must NOT re-derive, lifted from V1 `scripts/R/find_introgression.R`
# (itself a port of `scripts/legacy/find_introgression_updated_framework.R`):
#
#   1. window binning + per-(sample, window) percent-mismatch-to-consensus
#      distance                                   (find_introgression.R ~182-223)
#   2. the 2D kernel-density contour extraction — build a `geom_density_2d`
#      plot, `ggplot_build()` it, pull the contour polygons + the point
#      coordinates back out                       (find_introgression.R ~226-252)
#
# Everything cohort-specific (three welded-in cluster names, two fixed
# distance axes) is gone; the maths is byte-for-byte the same operations on
# a generic (Kx, Ky) pair.
#
# Two deliberate, documented departures — structure only, not statistics:
#
#   a. V1 rebuilt the whole ggplot once per sample inside
#      `find_introgressed_regions(SAMPLENAME)`, with the focal sample as the
#      point layer. The contour layer's data is the FULL table every time and
#      `coord_cartesian()` is a zoom (it does not touch scale limits), so the
#      kde2d grid limits — `scales$x$dimension()` over the union of both
#      layers — are identical whether the point layer holds one sample or all
#      of them (layer 1's rows are a subset of layer 2's). We therefore build
#      the plot ONCE with every point in the point layer and test all points
#      against the same contours. Identical density estimate, O(1) instead of
#      O(n_samples) ggplot builds.
#   b. V1 passed every contour vertex above the level cut to a single
#      `sp::point.in.polygon()` call, concatenating nested rings and separate
#      pieces into one vertex sequence — ill-defined for nested contours. We
#      test each ring separately and take the maximum containing level
#      (`contour_max_level()` in introgression_detect.R). "Inside at level >= L"
#      is then exact, and the `relative` rule gets the density layer it needs.
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

# --------------------------------------------------------------------------
# Core 1a — window binning (verbatim from find_introgression.R:182)
# --------------------------------------------------------------------------

#' Bin a position onto its window midpoint.
#'
#' @param pos Integer vector of 1-based positions.
#' @param ws  Window size in bp (config `introgression.window_size_bp`).
#' @return    Numeric vector of window midpoints: floor(pos/ws)*ws + ws/2.
window_bin <- function(pos, ws) (floor(pos / ws) * ws) + (ws %/% 2)

#' Stable, pair-invariant window identifier.
#'
#' V1 numbered windows with a `row_number()` over the sorted (CHROM, bin)
#' pairs, which only stays consistent if every run sees exactly the same
#' variant set. Pairwise runs each see the same variants but we do not want to
#' depend on that: the id is deterministic from the coordinate alone, so
#' per-pair calls join across pairs without a shared index table.
#'
#' Format is `w<chrom>_<bin>` — the leading letter matters. A bare
#' "<chrom>:<bin>" is parsed as a clock time by readr's type guesser
#' (`1:5000` -> 01:50:00), which silently corrupts the id on round-trip.
window_id <- function(chrom, bin) paste0("w", chrom, "_", sprintf("%.0f", as.numeric(bin)))

# --------------------------------------------------------------------------
# Core 1b — cluster consensus (dominant) alleles
# (find_introgression.R:141-157, generalised from 3 named clusters to K)
# --------------------------------------------------------------------------

#' Per-(CHROM, POS, cluster) dominant allele, wide by cluster.
#'
#' Ties are dropped exactly as V1 did: V1 kept positions whose post-`max()`
#' row count was `< 4` with three clusters — i.e. one row per cluster and no
#' tie anywhere. The generic form of that condition is
#' `n_rows <= n_clusters`.
#'
#' @param gt_long Long genotype table: SAMPLE, CHROM, POS, SNP, Cluster.
#'                SNP is the hmmIBD integer encoding (-1 = missing).
#' @param clusters Character vector of the cluster labels to compute.
#' @return Tibble: CHROM, POS, one integer column per cluster label.
dominant_alleles <- function(gt_long, clusters) {
  n_clusters <- length(clusters)
  wide <- gt_long %>%
    filter(SNP >= 0) %>%
    group_by(CHROM, POS, Cluster, SNP) %>%
    summarise(Allele_count = n(), .groups = "drop") %>%
    group_by(CHROM, POS, Cluster) %>%
    filter(Allele_count == max(Allele_count)) %>%
    ungroup() %>%
    group_by(CHROM, POS) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    filter(n <= n_clusters) %>%          # V1's `filter(n < 4)` at K = 3
    dplyr::select(CHROM, POS, Cluster, SNP) %>%
    pivot_wider(names_from = Cluster, values_from = SNP)

  # A cluster with no non-missing call at a position is absent from the pivot
  # for that row; make the column exist so downstream code is uniform (NA
  # propagates through the distance sum exactly as it did in V1).
  for (cl in clusters) {
    if (!(cl %in% names(wide))) wide[[cl]] <- NA_integer_
  }
  wide
}

# --------------------------------------------------------------------------
# Core 1c — per-(sample, window) distance to each cluster consensus
# (find_introgression.R:206-223)
# --------------------------------------------------------------------------

#' Percent-mismatch distance from each sample-window to two consensuses.
#'
#' Preserved from V1: windows are kept only when the sample has MORE than
#' `min_snps` non-missing calls in them (`filter(n > min_snps)`, a strict
#' `>` — 5 in config means "at least 6"), and sample-windows whose two
#' distances are exactly equal are dropped (V1's
#' `filter(Mn_distance != Mf_distance)`) — they carry no directional signal
#' and sit on the diagonal where the two clouds meet.
#'
#' @param gt_long   Long genotypes (SAMPLE, CHROM, POS, SNP, Cluster).
#' @param consensus Output of `dominant_alleles()`.
#' @param kx,ky     The two cluster labels (column names in `consensus`).
#' @param window_size_bp,min_snps  Config values.
#' @return Tibble: SAMPLE, Cluster, CHROM, BIN, WINDOW, N_SNPS, DIST_X, DIST_Y.
window_distances <- function(gt_long, consensus, kx, ky,
                             window_size_bp, min_snps) {
  gt_long %>%
    left_join(consensus, by = c("CHROM", "POS")) %>%
    filter(SNP >= 0) %>%
    mutate(BIN = window_bin(POS, window_size_bp)) %>%
    group_by(SAMPLE, Cluster, CHROM, BIN) %>%
    mutate(n = n()) %>%
    filter(n > min_snps) %>%
    summarise(
      N_SNPS = dplyr::n(),
      DIST_X = sum(SNP != .data[[kx]], na.rm = TRUE) / dplyr::n() * 100,
      DIST_Y = sum(SNP != .data[[ky]], na.rm = TRUE) / dplyr::n() * 100,
      .groups = "drop"
    ) %>%
    filter(DIST_X != DIST_Y) %>%
    mutate(WINDOW = window_id(CHROM, BIN)) %>%
    dplyr::select(SAMPLE, Cluster, CHROM, BIN, WINDOW, N_SNPS, DIST_X, DIST_Y)
}

# --------------------------------------------------------------------------
# Core 2 — 2D kernel-density contours via ggplot_build()
# (find_introgression.R:226-252)
# --------------------------------------------------------------------------

# Fixed hexes so contour rows can be mapped back to their cluster the way V1
# did (`filter(colour == "#440154FF")`). Two entries is all a pair needs;
# the values are V1's first two viridis picks.
PAIR_COLOURS <- c("#440154FF", "#73D055FF")

#' Extract per-cluster density contours (and the built point coordinates).
#'
#' The plot is exactly V1's: a point layer over the sample-windows, a
#' `geom_density_2d(contour_var = "density")` layer coloured by cluster, a
#' manual colour scale, and `coord_cartesian(0-100, 0-100)` — the distances
#' are percentages so the panel is the full percent square.
#'
#' @param dist_tbl Output of `window_distances()` (needs DIST_X, DIST_Y, Cluster).
#' @param kx,ky    Cluster labels; kx gets PAIR_COLOURS[1], ky PAIR_COLOURS[2].
#' @return list(points = tibble(x, y), contours = tibble(cluster, level, ring, x, y))
density_contours <- function(dist_tbl, kx, ky) {
  cols <- setNames(PAIR_COLOURS, c(kx, ky))

  raster_plot <- dist_tbl %>%
    ggplot(aes(x = DIST_X, y = DIST_Y, group = Cluster)) +
    geom_point(data = dist_tbl, size = 0.75, alpha = 0.75) +
    geom_density_2d(mapping = aes(x = DIST_X, y = DIST_Y, colour = Cluster),
                    data = dist_tbl, contour_var = "density") +
    scale_colour_manual(values = cols) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100))

  pd <- ggplot_build(raster_plot)
  pts  <- pd$data[[1]]
  cont <- pd$data[[2]]

  if (!all(c("level", "x", "y", "group", "colour") %in% names(cont))) {
    stop("[introgression_core] unexpected ggplot_build contour columns: ",
         paste(names(cont), collapse = ", "))
  }

  # `group` is unique per (level, piece) — one closed ring each.
  hex_to_cluster <- setNames(names(cols), unname(cols))
  contours <- cont %>%
    as_tibble() %>%
    mutate(cluster = hex_to_cluster[as.character(colour)]) %>%
    filter(!is.na(cluster)) %>%
    transmute(cluster, level = as.numeric(level),
              ring = as.character(group), x, y)

  list(points   = tibble(x = pts$x, y = pts$y),
       contours = contours)
}
