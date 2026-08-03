#!/usr/bin/env Rscript
# plot_clonal_temporal.R  (agnostic — driven by clonal_clusters.tsv + date role)
# --------------------------------------------------------------------------
# Stacked bar of clonal-group sampling dates. Consumes the computed
# clonal_group column from clonal_clusters.tsv joined to the canonical
# `date` role — V1's hardcoded 8-sample tribble and Sabang-xlsx lookup
# do NOT port.
#
# Colour scheme uses the palette generator over whichever clonal_group
# values are present in the data.
#
# If a `date` column is absent (either because the role is null or every
# sample in a clonal group has NA date), the script writes an empty
# placeholder PNG/SVG with an explanatory message rather than crashing —
# so downstream FINAL_TARGETS still resolve. (The rule that calls this
# script is only wired in when the date role is set, so this is mostly a
# belt-and-braces guard.)
#
# Usage:
#   Rscript plot_clonal_temporal.R <clonal_clusters.tsv> <out.png> <out.svg>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})
try(library(ggbreak), silent = TRUE)   # optional axis break

script_args <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- dirname(sub("^--file=", "",
                          script_args[grep("^--file=", script_args)][1]))
source(file.path(SCRIPT_DIR, "palettes.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: plot_clonal_temporal.R <clonal_clusters.tsv> <png> <svg>")
}
tsv_in  <- args[1]
png_out <- args[2]
svg_out <- args[3]

blank_plot <- function(msg) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = msg, size = 5) +
    theme_void()
}

cc <- read_tsv(tsv_in, show_col_types = FALSE)

if (!("date" %in% names(cc))) {
  message("[plot_clonal_temporal] no `date` column in clonal_clusters.tsv — skipping")
  p <- blank_plot("No date role configured — temporal plot skipped.")
  ggsave(png_out, p, width = 8, height = 4, dpi = 150)
  ggsave(svg_out, p, width = 8, height = 4)
  quit(status = 0)
}

dated <- cc %>%
  filter(!is.na(clonal_group), !is.na(date), nzchar(as.character(date))) %>%
  mutate(date = suppressWarnings(as.Date(date)))

# Fallback: NA-date samples in clonal groups
n_group_no_date <- cc %>%
  filter(!is.na(clonal_group), is.na(date) | !nzchar(as.character(date))) %>%
  nrow()
if (n_group_no_date > 0) {
  message(sprintf(
    "[plot_clonal_temporal] %d clonal-group sample(s) have no date — excluded from the temporal plot",
    n_group_no_date
  ))
}

if (nrow(dated) == 0) {
  message("[plot_clonal_temporal] no dated clonal-group samples — writing placeholder")
  p <- blank_plot("No clonal-group samples have a usable date.\nTemporal plot skipped.")
  ggsave(png_out, p, width = 8, height = 4, dpi = 150)
  ggsave(svg_out, p, width = 8, height = 4)
  quit(status = 0)
}

group_levels <- sort(unique(dated$clonal_group))
group_pal    <- discrete_palette(group_levels, option = "D")

by_date <- dated %>% count(date, clonal_group)

p <- ggplot(by_date, aes(x = date, y = n, fill = clonal_group)) +
  geom_col(width = 4) +
  scale_fill_manual(values = group_pal, name = "Clonal group") +
  labs(x = "Date of collection", y = "Number of samples") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_line(colour = "grey80"),
    panel.grid.minor   = element_blank(),
    plot.title         = element_blank(),
    plot.subtitle      = element_blank()
  )

# Optional axis break: if the largest gap between consecutive dated samples
# is > 6 months, break there (matches V1's "long gap in the middle" look).
# Only applied when ggbreak is available AND a genuine gap exists.
if (requireNamespace("ggbreak", quietly = TRUE)) {
  d_sorted <- sort(unique(by_date$date))
  if (length(d_sorted) >= 2) {
    diffs <- diff(d_sorted)
    biggest <- which.max(diffs)
    if (diffs[biggest] > 180) {
      brk_lo <- d_sorted[biggest]      + 15
      brk_hi <- d_sorted[biggest + 1]  - 15
      if (brk_hi > brk_lo) {
        p <- p + ggbreak::scale_x_break(c(brk_lo, brk_hi), space = 1) +
                 scale_x_date(date_labels = "%d-%b-%Y")
        # ggbreak adds a mirrored top axis — hide it
        p <- p + theme(axis.text.x.top = element_blank())
      }
    }
  }
}

ggsave(png_out, p, width = 10, height = 5, dpi = 300, limitsize = FALSE)
ggsave(svg_out, p, width = 10, height = 5,            limitsize = FALSE)

cat("\n=== Clonal-group counts (dated) ===\n")
print(dated %>% count(clonal_group, sort = TRUE))
message("[plot_clonal_temporal] wrote ", png_out, " + ", svg_out)
