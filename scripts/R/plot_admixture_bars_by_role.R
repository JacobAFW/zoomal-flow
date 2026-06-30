#!/usr/bin/env Rscript
# plot_admixture_bars_by_role.R  (agnostic — facet by an arbitrary role)
# --------------------------------------------------------------------------
# Faceted stacked ADMIXTURE bars; one panel per level of the given role
# column (typically "country" or "geography" from the canonical
# metadata). Samples within each panel are ordered by dominant cluster
# then by dominant proportion.
#
# Drops V1's fixed Indonesia/Malaysia / 6-province assumptions — the
# panel set comes from whatever levels exist in the role column.
#
# Usage:
#   Rscript plot_admixture_bars_by_role.R <Q> <fam> <admix_clusters.tsv> \
#       <metadata.tsv> <role_col> <png> <svg>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

script_args <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- dirname(sub("^--file=", "",
                          script_args[grep("^--file=", script_args)][1]))
source(file.path(SCRIPT_DIR, "palettes.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop(paste("Usage: Rscript plot_admixture_bars_by_role.R",
             "<Q> <fam> <admix_clusters.tsv> <metadata.tsv>",
             "<role_col> <png> <svg>"))
}
q_in        <- args[1]
fam_in      <- args[2]
clusters_in <- args[3]
meta_in     <- args[4]
role_col    <- args[5]
png_out     <- args[6]
svg_out     <- args[7]

theme_pub_min <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid      = element_blank(),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      legend.position = "right",
      strip.text      = element_text(face = "bold"),
      plot.margin     = margin(6, 12, 6, 6)
    )
}

q   <- read.table(q_in, header = FALSE) %>% as_tibble()
fam <- read.table(fam_in, header = FALSE, stringsAsFactors = FALSE) %>%
  as_tibble() %>% select(Sample = V1)
K   <- ncol(q)
component_cols <- paste0("V", seq_len(K))

admix <- bind_cols(fam, q) %>%
  rename_with(~ component_cols, all_of(seq_len(K) + 1))

clusters <- read_tsv(clusters_in, show_col_types = FALSE) %>%
  select(Sample, Cluster, Component_dom = Component, Proportion_dom = Proportion)

meta <- read_tsv(meta_in, show_col_types = FALSE)
stopifnot("sample_id" %in% names(meta))
if (!(role_col %in% names(meta))) {
  stop(sprintf("metadata has no column '%s' — was the role configured?", role_col))
}

meta_role <- meta %>%
  select(Sample = sample_id, role_value = all_of(role_col)) %>%
  filter(!is.na(role_value), nzchar(role_value))

# Panel ordering: largest cohort first (largest n samples).
panel_order <- meta_role %>%
  count(role_value) %>%
  arrange(desc(n)) %>%
  pull(role_value)

sample_order <- admix %>%
  left_join(clusters, by = "Sample") %>%
  left_join(meta_role, by = "Sample") %>%
  filter(!is.na(role_value)) %>%
  mutate(role_value = factor(role_value, levels = panel_order)) %>%
  arrange(role_value, Cluster, desc(Proportion_dom), Sample) %>%
  pull(Sample)

long <- admix %>%
  pivot_longer(all_of(component_cols),
               names_to = "Component", values_to = "Proportion") %>%
  inner_join(meta_role, by = "Sample") %>%
  mutate(
    Sample     = factor(Sample, levels = sample_order),
    role_value = factor(role_value, levels = panel_order)
  )

comp_pal <- discrete_palette(component_cols, option = "D")

# Free x scale so each panel only shows its own samples + auto-sizes width.
# Panel rows are auto-laid; choose ncol so wide panels (many samples) get a
# row of their own.
n_panels <- length(panel_order)
ncol_facet <- if (n_panels <= 2) n_panels else if (n_panels <= 6) 3 else 4

p <- ggplot(long, aes(x = Sample, y = Proportion, fill = Component)) +
  geom_col(width = 1, colour = "white", linewidth = 0.08) +
  facet_wrap(vars(role_value), scales = "free_x", ncol = ncol_facet) +
  scale_fill_manual(values = comp_pal) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.001)) +
  labs(x = NULL, y = "Ancestry proportion") +
  theme_pub_min()

# Width scales with panel count: ~3.5" per panel column.
w <- max(8, 3.2 * ncol_facet)
h <- 3 * ceiling(n_panels / ncol_facet)

ggsave(png_out, p, width = w, height = h, dpi = 300, limitsize = FALSE)
ggsave(svg_out, p, width = w, height = h, limitsize = FALSE)
message(sprintf("[plot_admixture_bars_by_role(%s)] wrote %s + %s (%d panels)",
                role_col, png_out, svg_out, n_panels))
