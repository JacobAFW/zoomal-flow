#!/usr/bin/env Rscript
# plot_province_map.R  (agnostic — countries + extent derived from data)
# --------------------------------------------------------------------------
# Admin-1 (province / state) choropleth + jittered sample points coloured
# by cluster. The countries to draw and the map extent are both derived
# from the data — no hardcoded Malaysia / Indonesia / bbox values.
#
# Polygons: rnaturalearth::ne_states(country = <countries seen in data>).
# Region fill: by geography role via the palette generator.
# Points:    by Cluster from admix_clusters.tsv (via gis_join.R input).
#
# Inputs: admix_clusters_gis.tsv (output of gis_join.R) with columns
#         Sample, Cluster, country, geography, lat, long
#
# Usage:
#   Rscript plot_province_map.R <admix_clusters_gis.tsv> <out.png> <out.svg>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(rnaturalearth)
})

script_args <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- dirname(sub("^--file=", "",
                          script_args[grep("^--file=", script_args)][1]))
source(file.path(SCRIPT_DIR, "palettes.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript plot_province_map.R <gis.tsv> <png> <svg>")
gis_in  <- args[1]
png_out <- args[2]
svg_out <- args[3]

map_data <- read_tsv(gis_in, show_col_types = FALSE) %>%
  filter(!is.na(lat), !is.na(long)) %>%
  mutate(lat = as.numeric(lat), long = as.numeric(long))

if (nrow(map_data) == 0) {
  stop("[plot_province_map] no rows with lat/long; check gis_join output and metadata.gis")
}
if (!("country" %in% names(map_data))) {
  stop("[plot_province_map] gis_join output has no 'country' column — country role is required for the country-list lookup")
}

# Map extent from data bbox + 1° padding on each side.
bbox <- with(map_data, list(
  xmin = min(long, na.rm = TRUE) - 1.0,
  xmax = max(long, na.rm = TRUE) + 1.0,
  ymin = min(lat,  na.rm = TRUE) - 1.0,
  ymax = max(lat,  na.rm = TRUE) + 1.0
))

# Countries to fetch polygons for: derived from the country role.
countries <- map_data %>%
  distinct(country) %>%
  filter(!is.na(country), nzchar(country)) %>%
  pull(country)
if (length(countries) == 0) {
  stop("[plot_province_map] no countries in gis data — cannot determine polygon set")
}
message(sprintf("[plot_province_map] fetching admin-1 polygons for: %s",
                paste(countries, collapse = ", ")))

sf::sf_use_s2(FALSE)
provinces <- ne_states(country = tolower(countries), returnclass = "sf")

# Geography palette: levels = distinct geography values in the data.
geo_levels <- sort(unique(stats::na.omit(map_data$geography)))
geo_pal    <- role_palette(geo_levels, "geography")

# Cluster palette: levels = distinct Cluster values in the data.
cluster_levels <- sort(unique(stats::na.omit(map_data$Cluster)))
cluster_pal    <- role_palette(cluster_levels, "group")

# Join province polygons to a "geography" attribute via name match against
# the canonical geography column. rnaturalearth's `name` field is the admin-1
# name; we match on it. Provinces that don't appear in the data get NA fill
# (rendered with na.value below).
provinces_keyed <- provinces %>%
  mutate(geography = if_else(name %in% geo_levels, name, NA_character_))

p <- ggplot() +
  geom_sf(data = provinces_keyed,
          aes(fill = geography),
          colour = "grey25", linewidth = 0.15) +
  scale_fill_manual(values = geo_pal, na.value = "grey85",
                    name = "Geography", na.translate = FALSE) +
  ggnewscale::new_scale_colour() +
  geom_jitter(data = map_data,
              aes(x = long, y = lat, colour = Cluster),
              width = 0.25, height = 0.25, alpha = 0.85, size = 2) +
  scale_colour_manual(values = cluster_pal, na.value = "grey50",
                      name = "Cluster") +
  coord_sf(xlim = c(bbox$xmin, bbox$xmax),
           ylim = c(bbox$ymin, bbox$ymax),
           expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_linedraw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_blank()
  )

ggsave(png_out, p, width = 10, height = 8, dpi = 300, limitsize = FALSE)
ggsave(svg_out, p, width = 10, height = 8,            limitsize = FALSE)
message(sprintf("[plot_province_map] wrote %s + %s (%d countries, %d samples)",
                png_out, svg_out, length(countries), nrow(map_data)))

cat("\n=== Map sample counts by (country, geography) ===\n")
print(map_data %>% count(country, geography, sort = TRUE), n = Inf)
