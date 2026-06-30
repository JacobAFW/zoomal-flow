#!/usr/bin/env Rscript
# plot_admixture_bars.R  (agnostic — single-panel stacked bar)
# --------------------------------------------------------------------------
# Stacked ADMIXTURE bar plot at the best K, samples ordered by dominant
# component then dominant proportion. Reads the .Q, the cleaned .fam (for
# sample order), the admix_clusters.tsv (for the dominant-cluster sort),
# and renders PNG + SVG.
#
# Generic. No cohort hardcodes; ^SRR exclusion removed (those samples are
# pre-cleaned out of the input now; see data/INPUT_PREP.md).
#
# Usage:
#   Rscript plot_admixture_bars.R <Q> <fam> <admix_clusters.tsv> <png> <svg>
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
  stop("Usage: Rscript plot_admixture_bars.R <Q> <fam> <admix_clusters.tsv> <png> <svg>")
}
q_in       <- args[1]
fam_in     <- args[2]
clusters_in<- args[3]
png_out    <- args[4]
svg_out    <- args[5]

theme_pub_min <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid       = element_blank(),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      legend.position  = "right",
      plot.margin      = margin(6, 12, 6, 6)
    )
}

q   <- read.table(q_in, header = FALSE) %>% as_tibble()
fam <- read.table(fam_in, header = FALSE, stringsAsFactors = FALSE) %>%
  as_tibble() %>% select(Sample = V1)
K <- ncol(q)
component_cols <- paste0("V", seq_len(K))

admix <- bind_cols(fam, q) %>%
  rename_with(~ component_cols, all_of(seq_len(K) + 1))

clusters <- read_tsv(clusters_in, show_col_types = FALSE) %>%
  select(Sample, Cluster, Component_dom = Component, Proportion_dom = Proportion)

# Order samples by (dominant cluster) → (dominant proportion descending).
sample_order <- admix %>%
  left_join(clusters, by = "Sample") %>%
  arrange(Cluster, desc(Proportion_dom), Sample) %>%
  pull(Sample)

long <- admix %>%
  pivot_longer(all_of(component_cols),
               names_to = "Component", values_to = "Proportion") %>%
  mutate(Sample = factor(Sample, levels = sample_order))

# Component palette: use viridis qualitative across K.
comp_pal <- discrete_palette(component_cols, option = "D")

p <- ggplot(long, aes(x = Sample, y = Proportion, fill = Component)) +
  geom_col(width = 1, colour = "white", linewidth = 0.08) +
  scale_fill_manual(values = comp_pal) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.001)) +
  labs(x = NULL, y = "Ancestry proportion") +
  theme_pub_min()

ggsave(png_out, p, width = 10, height = 3.5, dpi = 300)
ggsave(svg_out, p, width = 10, height = 3.5)
message(sprintf("[plot_admixture_bars] wrote %s + %s", png_out, svg_out))
