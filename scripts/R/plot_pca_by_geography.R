#!/usr/bin/env Rscript
# plot_pca_by_geography.R  (agnostic — colour by geography role)
# --------------------------------------------------------------------------
# Same PC1-vs-PC2 scatter as plot_pca.R, but coloured by the geography
# role (canonical column name 'geography' in outputs/metadata/samples.tsv).
#
# Usage:
#   Rscript plot_pca_by_geography.R <Pk.eigenvec> <pca_variance.tsv> \
#       <metadata.tsv> <png> <svg>
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
  stop("Usage: Rscript plot_pca_by_geography.R <Pk.eigenvec> <pca_variance.tsv> <metadata.tsv> <png> <svg>")
}
eig_in    <- args[1]
varpct_in <- args[2]
meta_in   <- args[3]
png_out   <- args[4]
svg_out   <- args[5]

theme_pub_min <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid      = element_blank(),
      panel.border    = element_rect(fill = NA, colour = "grey60"),
      legend.position = "right"
    )
}

eig <- read_table(eig_in, show_col_types = FALSE)
names(eig)[1] <- sub("^#", "", names(eig)[1])
sample_col <- intersect(c("IID", "FID", "Sample", "sample"), names(eig))[1]
eig <- eig %>% rename(Sample = !!sym(sample_col))

meta <- read_tsv(meta_in, show_col_types = FALSE)
stopifnot("sample_id" %in% names(meta))
if (!("geography" %in% names(meta))) {
  stop("metadata has no 'geography' column — was metadata.roles.geography configured?")
}

meta_geo <- meta %>%
  select(Sample = sample_id, geography) %>%
  mutate(geography = if_else(is.na(geography) | !nzchar(geography), NA_character_, geography))

joined <- eig %>% left_join(meta_geo, by = "Sample")
stopifnot(all(c("PC1", "PC2") %in% names(joined)))

varpct <- read_tsv(varpct_in, show_col_types = FALSE)
pc1 <- varpct$variance_percent[varpct$PC == 1]
pc2 <- varpct$variance_percent[varpct$PC == 2]

levels_geo <- sort(unique(stats::na.omit(joined$geography)))
geo_pal    <- role_palette(levels_geo, "geography")

p <- ggplot(joined, aes(PC1, PC2, colour = geography)) +
  geom_point(alpha = 0.85, size = 2) +
  scale_colour_manual(values = geo_pal, na.value = "grey70", name = "Geography") +
  labs(
    x = sprintf("PC1 (%.1f%%)", pc1),
    y = sprintf("PC2 (%.1f%%)", pc2)
  ) +
  theme_pub_min()

ggsave(png_out, p, width = 7.5, height = 5.5, dpi = 300)
ggsave(svg_out, p, width = 7.5, height = 5.5)
message(sprintf("[plot_pca_by_geography] wrote %s + %s", png_out, svg_out))
