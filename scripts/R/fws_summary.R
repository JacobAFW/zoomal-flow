#!/usr/bin/env Rscript
# fws_summary.R  (agnostic port — role-driven)
# --------------------------------------------------------------------------
# Per-country and per-geography Fws summaries. Reads the canonical
# role-renamed metadata table (columns: sample_id, country, geography, ...)
# and groups by whichever roles are present. Either output path can be
# the string "NULL" — that role's summary is then skipped (the validator
# at Stage 0 will already have noted the absent role).
#
# Usage:
#   Rscript fws_summary.R <fws.tsv> <metadata.tsv> \
#       <by_country.tsv|NULL> <by_geography.tsv|NULL> <cutoff> [exclude_list|NULL]
#
# DE-CLONALIZATION (optional 6th argument). An exclusion list drops samples
# COUNTED here (docs/clonality.md). The per-sample Fws itself is untouched —
# it does not depend on which other samples are present — but the by-group
# rates do: a clonal block of three polyclonal samples from one province
# counts that province's polyclonality three times. `full` passes every
# sample (a no-op); `unique` passes one representative per clonal group.
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(paste(
    "Usage: Rscript fws_summary.R",
    "<fws.tsv> <metadata.tsv>",
    "<by_country.tsv|NULL> <by_geography.tsv|NULL>",
    "<cutoff>"
  ))
}
fws_in    <- args[1]
meta_in   <- args[2]
out_cnt   <- args[3]   # "NULL" → skip
out_geo   <- args[4]   # "NULL" → skip
cutoff    <- as.numeric(args[5])
drop_in   <- if (length(args) >= 6) args[6] else "NULL"

fws  <- read_tsv(fws_in,  show_col_types = FALSE) %>% rename(Fws = Proportion)
# A DROP-list, not a keep-list: this stage scores every QC-surviving sample,
# a larger universe than the PLINK/clustered sets downstream. Intersecting
# against a keep-list drawn from one of those would silently shrink the `full`
# arm; removing named samples is a true no-op when the list is empty.
if (!(drop_in %in% c("NULL", "None", "")) && file.exists(drop_in)) {
  drop <- readLines(drop_in, warn = FALSE); drop <- drop[nzchar(drop)]
  if (length(drop) > 0) {
    n_before <- nrow(fws)
    fws <- fws %>% filter(!(sample %in% drop))
    message(sprintf("[fws_summary] exclusion list %s: dropped %d, %d samples counted",
                    basename(drop_in), n_before - nrow(fws), nrow(fws)))
  }
}
meta <- read_tsv(meta_in, show_col_types = FALSE)

# Canonical roles live as lowercase column names in samples.tsv. We only
# guarantee `sample_id`; country / geography may or may not be present.
stopifnot("sample_id" %in% names(meta))

combined <- fws %>%
  rename(sample_id = sample) %>%
  left_join(meta, by = "sample_id")

summarise_by <- function(df, group_cols, cutoff) {
  df %>%
    filter(if_all(all_of(group_cols), ~ !is.na(.x))) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n              = n(),
      fws_median     = median(Fws, na.rm = TRUE),
      fws_min        = min(Fws, na.rm = TRUE),
      fws_max        = max(Fws, na.rm = TRUE),
      n_polyclonal   = sum(Fws < cutoff, na.rm = TRUE),
      pct_polyclonal = 100 * mean(Fws < cutoff, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n))
}

if (out_cnt != "NULL") {
  if (!"country" %in% names(combined)) {
    message("[fws_summary] country role configured but column 'country' absent → skipping")
  } else {
    by_country <- summarise_by(combined, "country", cutoff)
    write_tsv(by_country, out_cnt)
    message(sprintf("[fws_summary] wrote %s", out_cnt))
    cat("\n=== Fws by country ===\n")
    print(by_country, n = Inf)
  }
} else {
  message("[fws_summary] no country role → skipping country summary")
}

if (out_geo != "NULL") {
  if (!"geography" %in% names(combined)) {
    message("[fws_summary] geography role configured but column 'geography' absent → skipping")
  } else {
    # Composite (country, geography) when country role is present — collision-
    # safe across cohorts that share a region name (and reproduces V1's
    # Sabah-by-country split: Malaysia 819 + Macaque 2). Fall back to
    # geography-only when country is absent.
    geo_keys <- if ("country" %in% names(combined)) {
      c("country", "geography")
    } else {
      "geography"
    }
    by_geo <- summarise_by(combined, geo_keys, cutoff)
    write_tsv(by_geo, out_geo)
    message(sprintf("[fws_summary] wrote %s (grouped by %s)",
                    out_geo, paste(geo_keys, collapse = " + ")))
    cat("\n=== Fws by ", paste(geo_keys, collapse = " + "), " ===\n", sep = "")
    print(by_geo, n = Inf)
  }
} else {
  message("[fws_summary] no geography role → skipping geography summary")
}
