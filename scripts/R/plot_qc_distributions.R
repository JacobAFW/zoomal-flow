#!/usr/bin/env Rscript
# plot_qc_distributions.R  (agnostic port)
# --------------------------------------------------------------------------
# Density / boxplot panels of MAF and per-sample/per-variant missingness from
# PLINK --freq / --missing output. Identical visual output to the V1 script,
# minus the cohort-specific contig-label sub() (V1 stripped 'ordered_PKNH_0?'
# and '_v2' to abbreviate contig names on the boxplot x-axis). Here the
# contig labels stay raw — generic across cohorts.
#
# Ported from: ../../scripts/R/plot_qc_distributions.R (V1).
#
# Usage:
#   Rscript plot_qc_distributions.R <Pk.frq> <Pk.imiss> <Pk.lmiss> \
#       <maf_overall.png> <maf_bychr.png> \
#       <miss_imiss.png> <miss_lmiss.png>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(viridis)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop(paste(
    "Usage: Rscript plot_qc_distributions.R",
    "<Pk.frq> <Pk.imiss> <Pk.lmiss>",
    "<maf_overall.png> <maf_bychr.png>",
    "<miss_imiss.png> <miss_lmiss.png>"
  ))
}
frq_in            <- args[1]
imiss_in          <- args[2]
lmiss_in          <- args[3]
maf_overall_out   <- args[4]
maf_bychr_out     <- args[5]
miss_imiss_out    <- args[6]
miss_lmiss_out    <- args[7]

# A minimal "publication" theme: cohort-neutral. V1 sources theme_pub.R; the
# agnostic tree carries no cohort palettes, so we inline a small theme here
# and keep the figures readable rather than styled.
theme_pub_min <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x  = element_text(size = base + 1),
      axis.ticks.x = element_line()
    )
}

# MAF: tolerate PLINK 1.9 (MAF, CHR) and PLINK 2 (#CHROM, ALT_FREQS) schemas.
frq <- read_table(frq_in, show_col_types = FALSE)
if ("ALT_FREQS" %in% names(frq)) {
  frq <- frq %>%
    mutate(
      MAF = pmin(as.numeric(ALT_FREQS), 1 - as.numeric(ALT_FREQS)),
      CHR = `#CHROM`
    )
} else if ("MAF" %in% names(frq)) {
  frq <- frq %>% mutate(MAF = as.numeric(MAF))
  if (!"CHR" %in% names(frq) && "#CHROM" %in% names(frq)) frq$CHR <- frq$`#CHROM`
} else {
  stop("Could not find MAF or ALT_FREQS column in frq file")
}

p_maf <- ggplot(frq, aes(x = MAF)) +
  geom_density(fill = viridis::viridis(1, option = "D", begin = 0.45), alpha = 0.6) +
  labs(x = "MAF", y = "Density") +
  theme_pub_min()

frq_box <- frq %>%
  mutate(Contig = factor(as.character(CHR), levels = sort(unique(as.character(CHR)))))

p_maf_by_chr <- ggplot(frq_box, aes(x = Contig, y = MAF, fill = Contig)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
  scale_fill_viridis_d(option = "H") +
  labs(x = "Contig", y = "MAF") +
  theme_pub_min() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(maf_overall_out, p_maf,        width = 8, height = 4, dpi = 150)
ggsave(maf_bychr_out,   p_maf_by_chr, width = 8, height = 4, dpi = 150)
message(sprintf("[plot_qc] wrote %s + %s", maf_overall_out, maf_bychr_out))

# Missingness — PLINK 1.9 (F_MISS) and PLINK 2 (MISSING_CT + OBS_CT) schemas.
imiss <- read_table(imiss_in, show_col_types = FALSE)
lmiss <- read_table(lmiss_in, show_col_types = FALSE)

normalise_missing <- function(df) {
  if ("F_MISS" %in% names(df)) {
    df %>% mutate(prop_missing = as.numeric(F_MISS))
  } else if (all(c("MISSING_CT", "OBS_CT") %in% names(df))) {
    df %>% mutate(
      prop_missing = as.numeric(MISSING_CT) /
        (as.numeric(MISSING_CT) + as.numeric(OBS_CT))
    )
  } else {
    stop("Could not find F_MISS or MISSING_CT/OBS_CT in missingness file")
  }
}
imiss <- normalise_missing(imiss)
lmiss <- normalise_missing(lmiss)

p_imiss <- ggplot(imiss, aes(x = prop_missing)) +
  geom_histogram(bins = 50,
                 fill = viridis::viridis(1, option = "D", begin = 0.55),
                 colour = "white") +
  labs(x = "Proportion missing", y = "Count") +
  theme_pub_min()

p_lmiss <- ggplot(lmiss, aes(x = prop_missing)) +
  geom_histogram(bins = 50,
                 fill = viridis::viridis(1, option = "D", begin = 0.95),
                 colour = "white") +
  labs(x = "Proportion missing", y = "Count") +
  theme_pub_min()

ggsave(miss_imiss_out, p_imiss, width = 8, height = 4, dpi = 150)
ggsave(miss_lmiss_out, p_lmiss, width = 8, height = 4, dpi = 150)
message(sprintf("[plot_qc] wrote %s + %s", miss_imiss_out, miss_lmiss_out))

cat("\n=== QC summary ===\n")
cat(sprintf("Variants:        %s\n", format(nrow(frq), big.mark = ",")))
cat(sprintf("Samples:         %s\n", format(nrow(imiss), big.mark = ",")))
cat(sprintf("MAF: median %.3f, mean %.3f\n",
            median(frq$MAF, na.rm = TRUE), mean(frq$MAF, na.rm = TRUE)))
cat(sprintf("Per-sample missing: median %.3f, max %.3f\n",
            median(imiss$prop_missing, na.rm = TRUE),
            max(imiss$prop_missing, na.rm = TRUE)))
cat(sprintf("Per-variant missing: median %.3f, max %.3f\n",
            median(lmiss$prop_missing, na.rm = TRUE),
            max(lmiss$prop_missing, na.rm = TRUE)))
