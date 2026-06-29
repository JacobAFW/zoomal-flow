#!/usr/bin/env Rscript
# plot_admixture_cv.R  (agnostic port — minimal cohort-neutral theme)
# --------------------------------------------------------------------------
# CV-error vs K diagnostic for ADMIXTURE. Highlights the best-K from
# best_k.txt with a vertical line so the auto-pick is visible against the
# curve. Renders PNG + SVG.
#
# Usage:
#   Rscript plot_admixture_cv.R <cv_error.tsv> <best_k.txt> <out.png> <out.svg>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript plot_admixture_cv.R <cv_error.tsv> <best_k.txt> <out.png> <out.svg>")
}
cv_in   <- args[1]
bk_in   <- args[2]
png_out <- args[3]
svg_out <- args[4]

cv     <- read_tsv(cv_in, show_col_types = FALSE)
best_k <- as.integer(readLines(bk_in, warn = FALSE)[1])

theme_pub_min <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x  = element_text(size = base),
      axis.ticks.x = element_line()
    )
}

p <- ggplot(cv, aes(K, CV_error)) +
  geom_line(linewidth = 0.5, colour = "grey40") +
  geom_point(size = 3, colour = "grey20") +
  geom_vline(xintercept = best_k, linetype = "dashed", colour = "firebrick") +
  annotate("text",
           x = best_k, y = max(cv$CV_error, na.rm = TRUE),
           label = sprintf("best K = %d", best_k),
           hjust = -0.15, vjust = 1, size = 3.5, colour = "firebrick") +
  scale_x_continuous(breaks = cv$K) +
  labs(x = "K (number of ancestry components)", y = "CV error") +
  theme_pub_min()

ggsave(png_out, p, width = 7, height = 4, dpi = 300)
ggsave(svg_out, p, width = 7, height = 4)
message(sprintf("[plot_admixture_cv] wrote %s + %s; best K = %d",
                png_out, svg_out, best_k))
