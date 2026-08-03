#!/usr/bin/env Rscript
# plot_ihs_scan.R  (agnostic — self-contained theme, generic CHR handling)
# --------------------------------------------------------------------------
# Genome-wide iHS scan + -log10(p) plot from the rehh output. Renders a
# placeholder if the ihs_table is empty (V1's Aceh model was a negative
# result — no candidate regions, empty scan). No cohort-specific chrom
# prefixes assumed.
#
# Usage:
#   Rscript plot_ihs_scan.R <ihs_table.tsv> \
#       <scan.png> <pvalue.png> <scan.svg> <pvalue.svg>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: plot_ihs_scan.R <ihs_table.tsv> <4 outputs>")
ihs_in   <- args[1]
scan_png <- args[2]
pval_png <- args[3]
scan_svg <- args[4]
pval_svg <- args[5]

theme_pub_min <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x  = element_text(size = base + 1),
      axis.ticks.x = element_line(),
      legend.position = "none"
    )
}

ihs <- read_tsv(ihs_in, show_col_types = FALSE)

blank <- function(msg) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = msg, size = 4.5) +
    theme_void()
}

if (nrow(ihs) == 0) {
  p <- blank("iHS scan is empty — model produced no per-SNP data.")
  for (out in c(scan_png, pval_png)) ggsave(out, p, width = 14, height = 4, dpi = 300)
  for (out in c(scan_svg, pval_svg)) ggsave(out, p, width = 14, height = 4)
  message("[plot_ihs_scan] empty input — wrote placeholders")
  quit(status = 0)
}

# CHR normalisation: attempt integer coercion for ordering; keep character
# levels for display. This works whether CHR is numeric already (Stage 3's
# plink2 --update-chr output) or a string ("ordered_PKNH_01_v2", …).
plot_data <- ihs %>%
  na.omit() %>%
  mutate(
    CHR_num = suppressWarnings(as.numeric(as.character(CHR))),
    CHR_lbl = as.character(CHR)
  ) %>%
  arrange(CHR_num, POSITION) %>%
  mutate(POS = seq_len(n())) %>%
  mutate(CHR_lbl = factor(CHR_lbl, levels = unique(CHR_lbl)))

x_axis <- plot_data %>%
  group_by(CHR_lbl) %>%
  summarise(POS = median(POS), .groups = "drop")

chrom_palette <- rep(c("#39568CFF", "#29AF7FFF"),
                     length.out = nlevels(plot_data$CHR_lbl))

hline <- plot_data %>%
  reframe(quantile = c("lower","median","upper"),
          IHS = quantile(IHS, c(0.001, 0.5, 0.999), na.rm = TRUE))

scan_p <- ggplot(plot_data, aes(x = POS, y = IHS, colour = CHR_lbl)) +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_colour_manual(values = chrom_palette) +
  scale_x_continuous(breaks = x_axis$POS, labels = x_axis$CHR_lbl) +
  geom_hline(yintercept = c(hline$IHS[1], hline$IHS[3]), linetype = "dotted") +
  labs(x = "Chromosome", y = "iHS") +
  theme_pub_min()

ggsave(scan_png, scan_p, width = 14, height = 4, dpi = 300)
ggsave(scan_svg, scan_p, width = 14, height = 4)

pval_p <- ggplot(plot_data, aes(x = POS, y = LOGPVALUE, colour = CHR_lbl)) +
  geom_jitter(size = 0.6, alpha = 0.7) +
  scale_colour_manual(values = chrom_palette) +
  scale_x_continuous(breaks = x_axis$POS, labels = x_axis$CHR_lbl) +
  labs(x = "Chromosome", y = expression(iHS ~ -log[10](italic(p)))) +
  theme_pub_min()

ggsave(pval_png, pval_p, width = 14, height = 4, dpi = 300)
ggsave(pval_svg, pval_p, width = 14, height = 4)

message(sprintf("[plot_ihs_scan] wrote %s + %s + %s + %s",
                scan_png, pval_png, scan_svg, pval_svg))
