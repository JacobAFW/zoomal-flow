#!/usr/bin/env Rscript
# gis_join.R  (agnostic port — role-driven, no _DK.* hardcode)
# --------------------------------------------------------------------------
# Join cluster assignments to a per-sample GIS table (sample_id + lat/long).
# The GIS TSV should have a sample-ID column (any case, may be named
# `Sample` / `sample` / `sample_id`) plus lat/long under one of the common
# spellings: lat/long, Lat/Long, Latitude/Longitude.
#
# Sample IDs are matched verbatim — no _DK.*-style strip. If your GIS table
# needs a strip, do it before saving the file, or use the same canonical
# IDs in both. (V1 stripped Illumina lane suffix; the agnostic version is
# stricter on IDs so cohorts that don't have that suffix don't get bitten.)
#
# Usage:
#   Rscript gis_join.R <admix_clusters.tsv> <gis_ref.tsv> <metadata.tsv> <out.tsv>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: see header")
clusters_in <- args[1]
gis_ref_in  <- args[2]
meta_in     <- args[3]
out_tsv     <- args[4]

clusters <- read_tsv(clusters_in, show_col_types = FALSE) %>%
  distinct(Sample, Cluster)

meta <- read_tsv(meta_in, show_col_types = FALSE)
stopifnot("sample_id" %in% names(meta))

# Optional country / geography columns; pass through when present.
meta_keep <- meta %>%
  select(any_of(c("sample_id", "country", "geography"))) %>%
  rename(Sample = sample_id) %>%
  distinct(Sample, .keep_all = TRUE)

gis_raw <- read_tsv(gis_ref_in, show_col_types = FALSE)

# Harmonise sample-ID column name → "Sample"
id_candidates <- intersect(c("Sample", "sample", "sample_id", "ID", "id"),
                           names(gis_raw))
if (length(id_candidates) == 0) {
  stop("GIS file has no sample-ID column (expected one of: Sample, sample, sample_id, ID, id)")
}
gis_raw <- gis_raw %>% rename(Sample = !!sym(id_candidates[1]))

# Harmonise lat/long column names
lat_candidates  <- intersect(c("lat", "Lat", "Latitude"),  names(gis_raw))
long_candidates <- intersect(c("long", "Long", "Longitude"), names(gis_raw))
if (length(lat_candidates) == 0 || length(long_candidates) == 0) {
  stop("GIS file is missing a lat/long column pair")
}
gis <- gis_raw %>%
  mutate(
    lat  = suppressWarnings(as.numeric(.data[[lat_candidates[1]]])),
    long = suppressWarnings(as.numeric(.data[[long_candidates[1]]]))
  ) %>%
  distinct(Sample, .keep_all = TRUE) %>%
  select(Sample, lat, long)

# Sanity filter: drop rows with physically-impossible coords (lat/long
# swap typos in the source file). Lat in [-90, 90], long in [-180, 180].
bad_coords <- gis %>%
  filter(!is.na(lat) & !is.na(long) &
         (abs(lat) > 90 | abs(long) > 180))
if (nrow(bad_coords) > 0) {
  message(sprintf("[gis_join] %d sample(s) have invalid lat/long (likely swapped) — dropping:",
                  nrow(bad_coords)))
  print(bad_coords)
  gis <- gis %>%
    mutate(
      lat  = if_else(abs(lat) > 90 | abs(long) > 180, NA_real_, lat),
      long = if_else(abs(lat) > 90 | abs(long) > 180, NA_real_, long)
    )
}

joined <- clusters %>%
  left_join(meta_keep, by = "Sample") %>%
  left_join(gis,       by = "Sample")

write_tsv(joined, out_tsv)

cat("\n=== GIS join coverage ===\n")
print(joined %>%
        summarise(
          total            = n(),
          with_lat_long    = sum(!is.na(lat) & !is.na(long)),
          missing_lat_long = sum(is.na(lat) | is.na(long))
        ))
if ("country" %in% names(joined)) {
  cat("\nBy country (samples with lat/long):\n")
  print(joined %>%
          filter(!is.na(lat) & !is.na(long)) %>%
          count(country, sort = TRUE))
}
