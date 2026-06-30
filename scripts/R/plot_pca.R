#!/usr/bin/env Rscript
# plot_pca.R  (agnostic — coloured by the cluster column from admix_clusters)
# --------------------------------------------------------------------------
# Two-axis PCA scatter (PC1 vs PC2), coloured by ADMIXTURE cluster label
# from admix_clusters.tsv. Variance percent in axis labels comes from
# pca_variance.tsv.
#
# Usage:
#   Rscript plot_pca.R <Pk.eigenvec> <pca_variance.tsv> <admix_clusters.tsv> \
#       <png> <svg>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

script_args <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- dirname(sub("^--file=", "",
                          script_args[grep("^--file=", script_args)][1]))
source(file.path(SCRIPT_DIR, "palettes.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript plot_pca.R <Pk.eigenvec> <pca_variance.tsv> <admix_clusters.tsv> <png> <svg>")
}
eig_in     <- args[1]
varpct_in  <- args[2]
clusters_in<- args[3]
png_out    <- args[4]
svg_out    <- args[5]

theme_pub_min <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_rect(fill = NA, colour = "grey60"),
      legend.position  = "right"
    )
}

# PLINK 2 eigenvec has #FID + IID + PC1..PCn columns; header may or may
# not start with '#'.
eig <- read_table(eig_in, show_col_types = FALSE)
names(eig)[1] <- sub("^#", "", names(eig)[1])

clusters <- read_tsv(clusters_in, show_col_types = FALSE) %>%
  select(Sample, Cluster)

# Match sample-id column irrespective of name (IID or FID).
sample_col <- intersect(c("IID", "FID", "Sample", "sample"), names(eig))[1]
eig <- eig %>% rename(Sample = !!sym(sample_col))

joined <- eig %>% left_join(clusters, by = "Sample")
stopifnot(all(c("PC1", "PC2") %in% names(joined)))

varpct <- read_tsv(varpct_in, show_col_types = FALSE)
pc1 <- varpct$variance_percent[varpct$PC == 1]
pc2 <- varpct$variance_percent[varpct$PC == 2]

levels_cluster <- sort(unique(stats::na.omit(joined$Cluster)))
cluster_pal    <- role_palette(levels_cluster, "group")

p <- ggplot(joined, aes(PC1, PC2, colour = Cluster)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_colour_manual(values = cluster_pal, na.value = "grey70") +
  labs(
    x = sprintf("PC1 (%.1f%%)", pc1),
    y = sprintf("PC2 (%.1f%%)", pc2)
  ) +
  theme_pub_min()

ggsave(png_out, p, width = 7, height = 5.5, dpi = 300)
ggsave(svg_out, p, width = 7, height = 5.5)
message(sprintf("[plot_pca] wrote %s + %s", png_out, svg_out))
