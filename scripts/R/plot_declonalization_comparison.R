#!/usr/bin/env Rscript
# plot_declonalization_comparison.R — full vs unique, at a glance
# (docs/clonality.md)
# --------------------------------------------------------------------------
# WHAT: one dumbbell panel per stage. Each row is a metric; the two dots are
#       its full-set and de-clonalized values, joined by a segment whose
#       length IS the effect of collapsing clonal replicates.
# WHY:  the comparison TSVs are the record, but nobody scans four TSVs before
#       reading a result. This is the "how much did clonality drive this?"
#       glance — a row with the two dots on top of each other is a result
#       clonality did not move; a long segment is one that was substantially
#       pseudo-replication.
#
# Metrics are plotted on a per-panel free scale because they are not
# commensurable (a window count and a percentage do not share an axis), and
# the row label carries the raw pair so the numbers are readable off the
# figure without going back to the TSV.
#
# Usage:
#   Rscript plot_declonalization_comparison.R <out.png> <out.svg> <cmp1.tsv> [<cmp2.tsv> ...]
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_declonalization_comparison.R <out.png> <out.svg> <cmp.tsv> ...")
}
out_png <- args[1]
out_svg <- args[2]
in_tsvs <- args[-(1:2)]
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

cmp <- in_tsvs %>%
  keep(file.exists) %>%
  map(~ read_tsv(.x, show_col_types = FALSE)) %>%
  bind_rows()

# An empty or all-NA comparison is a legitimate state (no clonal signal, or a
# stage that produced nothing). Emit a placeholder rather than failing the
# workflow on a figure.
placeholder <- function(msg) {
  g <- ggplot() + annotate("text", x = 0, y = 0, label = msg, size = 4.5) +
    theme_void()
  ggsave(out_png, g, width = 8, height = 3, dpi = 300)
  ggsave(out_svg, g, width = 8, height = 3)
  message("[plot_declonalization_comparison] ", msg)
  quit(status = 0)
}
if (nrow(cmp) == 0) placeholder("No de-clonalization comparison data")

# Drop rows where nothing can be drawn, and the two bookkeeping rows that are
# identical in every panel (they are reported in the TSV and in the caption).
plot_df <- cmp %>%
  filter(!is.na(full), !is.na(unique)) %>%
  filter(!metric %in% c("n_samples_removed", "n_clonal_groups_collapsed")) %>%
  filter(!(full == 0 & unique == 0))
if (nrow(plot_df) == 0) placeholder("No comparable metrics between full and unique")

removed <- cmp %>% filter(metric == "n_samples_removed") %>% pull(unique) %>% first()
groups  <- cmp %>% filter(metric == "n_clonal_groups_collapsed") %>% pull(unique) %>% first()

fmt <- function(x) ifelse(x == round(x), sprintf("%d", as.integer(x)), sprintf("%.2f", x))

plot_df <- plot_df %>%
  mutate(label = sprintf("%s  (%s → %s)", metric, fmt(full), fmt(unique)),
         moved = abs(delta) > 0) %>%
  group_by(stage) %>%
  # Biggest movers at the top of each panel — that is what the figure is for.
  mutate(label = fct_reorder(label, abs(replace_na(pct_change, 0)))) %>%
  ungroup()

long <- plot_df %>%
  pivot_longer(c(full, unique), names_to = "sampleset", values_to = "value") %>%
  mutate(sampleset = factor(sampleset, levels = c("full", "unique")))

PAL <- c(full = "#4C6EF5", unique = "#F08C00")

g <- ggplot(long, aes(x = value, y = label)) +
  geom_segment(data = plot_df,
               aes(x = full, xend = unique, y = label, yend = label,
                   linetype = moved),
               inherit.aes = FALSE, colour = "grey55", linewidth = 0.6,
               show.legend = FALSE) +
  # `full` is a hollow ring and `unique` a filled dot, deliberately: when a
  # metric did not move the two coincide, and a plain overplot would hide one
  # of them entirely — "dot inside a ring" reads unambiguously as "no change",
  # which is the single most common thing this figure has to say.
  geom_point(data = ~ dplyr::filter(.x, sampleset == "full"),
             aes(colour = sampleset), size = 4.2, shape = 21, stroke = 1.2,
             fill = NA) +
  geom_point(data = ~ dplyr::filter(.x, sampleset == "unique"),
             aes(colour = sampleset), size = 2.3) +
  scale_colour_manual(values = PAL, name = "sample set",
                      labels = c(full = "full (all samples)",
                                 unique = "unique (de-clonalized)"),
                      guide = guide_legend(override.aes = list(
                        shape = c(21, 19), size = c(4.2, 2.3), fill = NA))) +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "blank")) +
  facet_wrap(~ stage, scales = "free", ncol = 2) +
  labs(
    title = "Effect of de-clonalization on each stage",
    subtitle = sprintf(
      "%s sample(s) removed by collapsing %s clonal group(s) to one representative each.\nA row whose two dots coincide is a result clonality did not drive.",
      if (is.na(removed)) "?" else fmt(removed),
      if (is.na(groups)) "?" else fmt(groups)),
    x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top",
        panel.grid.major.y = element_line(colour = "grey92"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey30"))

n_rows  <- n_distinct(plot_df$label)
n_panel <- n_distinct(plot_df$stage)
height  <- max(4, 1.6 + 0.22 * n_rows + 0.6 * ceiling(n_panel / 2))
ggsave(out_png, g, width = 11, height = height, dpi = 300, limitsize = FALSE)
ggsave(out_svg, g, width = 11, height = height, limitsize = FALSE)
message(sprintf("[plot_declonalization_comparison] wrote %s (%d metrics, %d stages)",
                out_png, n_rows, n_panel))
