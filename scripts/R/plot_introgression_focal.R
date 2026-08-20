#!/usr/bin/env Rscript
# plot_introgression_focal.R  (agnostic port of plot_introgression_aceh_unique.R)
# --------------------------------------------------------------------------
# The headline figure: every window introgressed uniquely in the focal group,
# placed on the genome, with its height = how many focal samples carry it.
#
#   x = genomic position, one panel per contig that has a hit
#   y = number of focal-group samples carrying the window
#
# Three V1 hardcodes removed: the literal "Aceh" labelling, the
# `levels = sprintf("Chr %02d", 1:14)` 14-chromosome factor, and the
# `source("scripts/R/theme_pub.R")` V1-relative path. Contigs come from the
# headline table's own CONTIG column (which the aggregate step resolves via
# contig_map.tsv), so a cohort with 3 or 30 contigs renders the same way.
#
# Usage:
#   Rscript plot_introgression_focal.R <headline.tsv> <focal_label> \
#       <out.png> <out.svg>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: plot_introgression_focal.R <headline.tsv> <focal_label> <png> <svg>")
}
tsv_in    <- args[1]
focal_lbl <- args[2]
png_out   <- args[3]
svg_out   <- args[4]

dir.create(dirname(png_out), recursive = TRUE, showWarnings = FALSE)

theme_pub_min <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          axis.ticks = element_line(),
          strip.text = element_text(size = 11, face = "bold"),
          legend.position = "right")
}

ud <- read_tsv(tsv_in, show_col_types = FALSE,
               col_types = cols(WINDOW = col_character(), .default = col_guess()))

if (nrow(ud) == 0) {
  p <- ggplot() +
    annotate("text", x = 0, y = 0, size = 4.5,
             label = sprintf("No windows unique to '%s'", focal_lbl)) +
    theme_void()
  width <- 8; height <- 4.5
} else {
  # Contig panels ordered by the pipeline's integer CHROM code (contig_map
  # order = reference .fai order), so panels read left-to-right as the
  # reference is laid out — no hardcoded chromosome list.
  contig_order <- ud %>% distinct(CHROM, CONTIG) %>% arrange(CHROM) %>% pull(CONTIG)
  plot_data <- ud %>%
    mutate(CONTIG   = factor(CONTIG, levels = contig_order),
           midpoint = (START + END) / 2)

  n_panels <- length(contig_order)
  ncol     <- min(4, n_panels)

  # A panel holding a single window has a zero-width x range, and ggplot's
  # default expansion for that is +/- 0.5 bp — every break then rounds to the
  # same Mb label. Pad by a few window widths so the ticks are distinguishable.
  pad <- 5 * max(median(plot_data$END - plot_data$START + 1), 1)

  p <- ggplot(plot_data, aes(x = midpoint, y = n_samples)) +
    geom_segment(aes(xend = midpoint, yend = 0), colour = "grey60") +
    geom_point(aes(fill = n_samples), shape = 21, size = 3, colour = "black") +
    facet_wrap(~ CONTIG, scales = "free_x", ncol = ncol) +
    scale_fill_viridis_c(option = "D", begin = 0.2, end = 0.9) +
    scale_x_continuous(labels = function(x) sprintf("%.2f Mb", x / 1e6),
                       n.breaks = 3,
                       expand = expansion(mult = 0.05, add = pad)) +
    labs(x = "Genomic position",
         y = sprintf("'%s' samples carrying the introgressed window", focal_lbl),
         fill = "n samples",
         title = sprintf("Windows introgressed uniquely in '%s' (%d window%s)",
                         focal_lbl, nrow(ud), if (nrow(ud) == 1) "" else "s")) +
    theme_pub_min()

  width  <- min(4 + 3.2 * ncol, 20)
  height <- min(3 + 2.4 * ceiling(n_panels / ncol), 16)
}

ggsave(png_out, p, width = width, height = height, dpi = 300, limitsize = FALSE)
ggsave(svg_out, p, width = width, height = height,            limitsize = FALSE)
cat(sprintf("[plot_introgression_focal] wrote %s + %s (%d window(s))\n",
            png_out, svg_out, nrow(ud)))
