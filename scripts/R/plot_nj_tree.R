#!/usr/bin/env Rscript
# plot_nj_tree.R  (agnostic — generic, tip-coloured by group/cluster)
# --------------------------------------------------------------------------
# Unrooted Neighbour-Joining tree from PLINK's pairwise IBS distance
# matrix. Tip points coloured by ADMIXTURE cluster (from admix_clusters.tsv);
# edges left a uniform light grey. Drops V1's Indonesia-highlight overlay
# and the hardcoded "Mf/Mn/Peninsular/Indonesia/Other" level list.
#
# Usage:
#   Rscript plot_nj_tree.R <Pk.dist> <Pk.dist.id> <admix_clusters.tsv> \
#       <png> <svg>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
  library(ggtree)
})

script_args <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- dirname(sub("^--file=", "",
                          script_args[grep("^--file=", script_args)][1]))
source(file.path(SCRIPT_DIR, "palettes.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript plot_nj_tree.R <Pk.dist> <Pk.dist.id> <admix_clusters.tsv> <png> <svg>")
}
dist_in    <- args[1]
distid_in  <- args[2]
clusters_in<- args[3]
png_out    <- args[4]
svg_out    <- args[5]

dm  <- as.matrix(read.table(dist_in, header = FALSE))
ids <- read.table(distid_in, header = FALSE, stringsAsFactors = FALSE)$V1
rownames(dm) <- colnames(dm) <- ids

tree <- ape::nj(as.dist(dm))

clusters <- read_tsv(clusters_in, show_col_types = FALSE) %>%
  select(Sample, Cluster)
tip_df <- tibble(label = tree$tip.label) %>%
  left_join(clusters, by = c("label" = "Sample")) %>%
  mutate(Cluster = if_else(is.na(Cluster), "Unassigned", Cluster))

levels_cluster <- sort(unique(tip_df$Cluster))
levels_cluster <- levels_cluster[levels_cluster != "Unassigned"]
cluster_pal    <- role_palette(levels_cluster, "group")
cluster_pal[["Unassigned"]] <- "grey75"

# Build the ggtree plot. Edges in light grey; tips coloured by Cluster
# via geom_tippoint, joined from tip_df by tip-label.
p <- ggtree(tree, layout = "daylight", colour = "grey60", size = 0.25) %<+%
  tip_df +
  geom_tippoint(aes(colour = Cluster), size = 1.6, alpha = 0.9) +
  scale_colour_manual(values = cluster_pal, na.value = "grey75",
                      name = "Cluster") +
  theme(legend.position = "right")

ggsave(png_out, p, width = 8, height = 7, dpi = 300)
ggsave(svg_out, p, width = 8, height = 7)
message(sprintf("[plot_nj_tree] wrote %s + %s", png_out, svg_out))
