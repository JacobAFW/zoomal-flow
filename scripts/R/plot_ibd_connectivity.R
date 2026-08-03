#!/usr/bin/env Rscript
# plot_ibd_connectivity.R  (agnostic — role-driven, no STATE_PAL)
# --------------------------------------------------------------------------
# Two side-by-side igraph panels for a focal cluster (IBD >= 5% and >= 95%).
# Nodes coloured by the `geography` role via the palette generator; no
# STATE_PAL hardcode. Outputs both PDF (legacy aesthetic) and PNG so the
# report and slide handoffs share the same figure.
#
# Usage:
#   Rscript plot_ibd_connectivity.R <hmm_fract.txt> <metadata.tsv> \
#       <cluster_label> <out.pdf> <out.png>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
})

script_args <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- dirname(sub("^--file=", "",
                          script_args[grep("^--file=", script_args)][1]))
source(file.path(SCRIPT_DIR, "palettes.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: plot_ibd_connectivity.R <fract> <metadata> <cluster> <pdf> <png>")
}
fract_in    <- args[1]
meta_in     <- args[2]
cluster_lbl <- args[3]
pdf_out     <- args[4]
png_out     <- args[5]

ibd <- read.delim(fract_in, stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  select(sample1, sample2, fract_sites_IBD)

meta <- read_tsv(meta_in, show_col_types = FALSE)
stopifnot("sample_id" %in% names(meta))
has_geo <- "geography" %in% names(meta)
if (!has_geo) {
  message("[plot_ibd_connectivity] geography role absent — nodes coloured uniformly")
}

all_samples <- unique(c(ibd$sample1, ibd$sample2))
vertex_df <- tibble(name = all_samples)
if (has_geo) {
  vertex_df <- vertex_df %>%
    left_join(
      meta %>% select(name = sample_id, geography),
      by = "name"
    )
  geo_levels <- sort(unique(stats::na.omit(vertex_df$geography)))
  geo_pal    <- role_palette(geo_levels, "geography")
  vertex_df <- vertex_df %>%
    mutate(
      geography = factor(geography, levels = geo_levels),
      Colour    = geo_pal[as.character(geography)]
    )
  vertex_df$Colour[is.na(vertex_df$Colour)] <- "#CDCDCD"
  legend_levels <- geo_levels
  legend_cols   <- geo_pal
} else {
  vertex_df <- vertex_df %>%
    mutate(geography = NA_character_, Colour = "#7FA5C0")
  legend_levels <- character(0)
  legend_cols   <- character(0)
}

cutoffs <- c(0.05, 0.95)

draw_panels <- function() {
  par(mfrow = c(1, 2), mar = rep.int(1, 4) + 0.1)
  for (i in seq_along(cutoffs)) {
    cutoff <- cutoffs[i]
    edges  <- ibd %>%
      filter(fract_sites_IBD >= cutoff) %>%
      select(sample1, sample2, weight = fract_sites_IBD)
    g <- igraph::graph_from_data_frame(edges, vertices = vertex_df,
                                       directed = FALSE)
    g <- igraph::set_vertex_attr(g, "color", value = igraph::V(g)$Colour)
    set.seed(2026)
    coords <- igraph::layout_nicely(g)
    plot(
      g, layout = coords,
      vertex.size = 4, vertex.label = NA,
      main = sprintf("%s — IBD >= %d%%", cluster_lbl, round(cutoff * 100))
    )
    if (i == 1 && length(legend_levels) > 0) {
      legend("bottomleft",
             legend = legend_levels,
             fill   = legend_cols,
             cex    = 0.8, bty = "n")
    }
  }
}

pdf(pdf_out, width = 14, height = 7)
draw_panels()
dev.off()

png(png_out, width = 14, height = 7, units = "in", res = 300)
draw_panels()
dev.off()

cat("\n=== Edge counts by cutoff ===\n")
for (cutoff in cutoffs) {
  n <- sum(ibd$fract_sites_IBD >= cutoff, na.rm = TRUE)
  cat(sprintf("  IBD >= %.2f -> %d pairs\n", cutoff, n))
}
if (has_geo) {
  cat("\n=== Vertex counts by geography ===\n")
  print(vertex_df %>% count(geography, sort = TRUE), n = Inf)
}

message("[plot_ibd_connectivity] wrote ", pdf_out, " + ", png_out)
