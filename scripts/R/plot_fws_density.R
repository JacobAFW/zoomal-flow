#!/usr/bin/env Rscript
# plot_fws_density.R  (agnostic port — role-driven)
# --------------------------------------------------------------------------
# Two-panel Fws density plots driven by the country / geography roles.
# Each panel pair (png + svg) is rendered only if the rule passed a real
# path for it; "NULL" skips. The V1 Indonesia-only province panel is
# dropped — Panel B is now "density by geography" across the cohort.
#
# Cohort-neutral theme (inline minimal — no V1-relative theme_pub.R source,
# consistent with Stage 1's plot_qc_distributions.R).
#
# Usage:
#   Rscript plot_fws_density.R <fws.tsv> <metadata.tsv> \
#       <country.png|NULL> <country.svg|NULL> \
#       <geography.png|NULL> <geography.svg|NULL> \
#       <cutoff>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(viridis)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop(paste(
    "Usage: Rscript plot_fws_density.R",
    "<fws.tsv> <metadata.tsv>",
    "<country.png|NULL> <country.svg|NULL>",
    "<geography.png|NULL> <geography.svg|NULL>",
    "<cutoff>"
  ))
}
fws_in    <- args[1]
meta_in   <- args[2]
cnt_png   <- args[3]
cnt_svg   <- args[4]
geo_png   <- args[5]
geo_svg   <- args[6]
cutoff    <- as.numeric(args[7])

theme_pub_min <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x  = element_text(size = base),
      axis.ticks.x = element_line(),
      legend.position = "right"
    )
}

fws  <- read_tsv(fws_in,  show_col_types = FALSE) %>% rename(Fws = Proportion)
meta <- read_tsv(meta_in, show_col_types = FALSE)

stopifnot("sample_id" %in% names(meta))

combined <- fws %>%
  rename(sample_id = sample) %>%
  left_join(meta, by = "sample_id")

density_panel <- function(df, role_col, palette_opt) {
  df_role <- df %>%
    filter(!is.na(.data[[role_col]]), .data[[role_col]] != "")
  # Order largest-n on bottom for readable legend
  level_order <- df_role %>%
    count(.data[[role_col]]) %>%
    arrange(desc(n)) %>%
    pull(.data[[role_col]])
  df_role[[role_col]] <- factor(df_role[[role_col]], levels = level_order)

  ggplot(df_role, aes(x = Fws,
                      fill = .data[[role_col]],
                      colour = .data[[role_col]])) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = cutoff, linetype = "dashed", colour = "grey30") +
    annotate(
      "text",
      x = cutoff, y = Inf, vjust = 1.8, hjust = 1.1,
      label = sprintf("polyclonal cutoff\nFws = %.2f", cutoff),
      size = 3, colour = "grey30"
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_fill_viridis_d(option = palette_opt, end = 0.85) +
    scale_colour_viridis_d(option = palette_opt, end = 0.85) +
    labs(x = "Fws", y = "Density", fill = role_col, colour = role_col) +
    theme_pub_min(12)
}

if (cnt_png != "NULL") {
  if (!"country" %in% names(combined)) {
    message("[plot_fws_density] country role configured but column 'country' absent → skipping")
  } else {
    p_country <- density_panel(combined, "country", "D")
    ggsave(cnt_png, p_country, width = 9, height = 4.5, dpi = 300)
    if (cnt_svg != "NULL") ggsave(cnt_svg, p_country, width = 9, height = 4.5)
    message(sprintf("[plot_fws_density] wrote %s + %s", cnt_png, cnt_svg))
  }
} else {
  message("[plot_fws_density] no country role → skipping country panel")
}

if (geo_png != "NULL") {
  if (!"geography" %in% names(combined)) {
    message("[plot_fws_density] geography role configured but column 'geography' absent → skipping")
  } else {
    p_geo <- density_panel(combined, "geography", "C")
    ggsave(geo_png, p_geo, width = 9, height = 4.5, dpi = 300)
    if (geo_svg != "NULL") ggsave(geo_svg, p_geo, width = 9, height = 4.5)
    message(sprintf("[plot_fws_density] wrote %s + %s", geo_png, geo_svg))
  }
} else {
  message("[plot_fws_density] no geography role → skipping geography panel")
}

# Console summary so the log captures the percentages
if ("country" %in% names(combined)) {
  cat("\n=== Polyclonality by country ===\n")
  combined %>%
    filter(!is.na(country), country != "") %>%
    group_by(country) %>%
    summarise(
      n              = n(),
      pct_polyclonal = sprintf("%.1f%%", 100 * mean(Fws < cutoff, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(desc(n)) %>%
    print(n = Inf)
}
if ("geography" %in% names(combined)) {
  cat("\n=== Polyclonality by geography ===\n")
  combined %>%
    filter(!is.na(geography), geography != "") %>%
    group_by(geography) %>%
    summarise(
      n              = n(),
      pct_polyclonal = sprintf("%.1f%%", 100 * mean(Fws < cutoff, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(desc(n)) %>%
    print(n = Inf)
}
