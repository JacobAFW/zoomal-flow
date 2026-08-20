#!/usr/bin/env Rscript
# plot_introgression_shoulder.R  (agnostic port of legacy shoulder_plot.png)
# --------------------------------------------------------------------------
# Windows ranked by how many samples carry them, PRE-filter. The curve's
# shoulder is where "a handful of samples share this window" turns into "one
# or two samples have it" — i.e. where to set
# `introgression.min_samples_per_window`. The current threshold is drawn as a
# dashed line so the figure shows whether the configured cut sits on the
# shoulder or somewhere arbitrary.
#
# Source: legacy scripts/legacy/introgression_multi_cluster.R:311-322.
# Role-driven: nothing here depends on cluster names or geography, so the
# figure renders for any cohort.
#
# Usage:
#   Rscript plot_introgression_shoulder.R <window_sample_counts_raw.tsv> \
#       <min_samples_per_window> <out.png> <out.svg>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: plot_introgression_shoulder.R <counts.tsv> <threshold> <png> <svg>")
}
counts_in <- args[1]
threshold <- suppressWarnings(as.numeric(args[2]))
png_out   <- args[3]
svg_out   <- args[4]

dir.create(dirname(png_out), recursive = TRUE, showWarnings = FALSE)

theme_pub_min <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          axis.ticks = element_line())
}

counts <- read_tsv(counts_in, show_col_types = FALSE,
                   col_types = cols(WINDOW = col_character(), .default = col_guess()))

blank <- function(msg) {
  ggplot() + annotate("text", x = 0, y = 0, label = msg, size = 4.5) + theme_void()
}

if (nrow(counts) == 0) {
  p <- blank("No introgressed windows called\n(nothing to rank)")
} else {
  plot_data <- counts %>%
    arrange(n_samples) %>%
    mutate(rank = row_number())

  n_above <- sum(counts$n_samples > threshold)
  p <- ggplot(plot_data, aes(x = rank, y = n_samples)) +
    geom_col(width = 1, fill = "#3B528BFF") +
    geom_hline(yintercept = threshold, linetype = "dashed", colour = "#B22222") +
    annotate("text", x = 0, y = threshold,
             label = sprintf("  min_samples_per_window = %g  (%d window%s kept)",
                             threshold, n_above, if (n_above == 1) "" else "s"),
             hjust = 0, vjust = -0.6, size = 3.4, colour = "#B22222") +
    labs(x = "Windows, ranked by number of samples carrying them",
         y = "Samples with the window called",
         title = sprintf("Introgressed-window frequency (%d windows, pre-filter)",
                         nrow(counts))) +
    theme_pub_min()
}

ggsave(png_out, p, width = 8, height = 4.5, dpi = 300)
ggsave(svg_out, p, width = 8, height = 4.5)
cat(sprintf("[plot_introgression_shoulder] wrote %s + %s\n", png_out, svg_out))
