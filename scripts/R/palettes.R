#!/usr/bin/env Rscript
# palettes.R  (agnostic — palette GENERATOR, not hardcoded vectors)
# --------------------------------------------------------------------------
# Given the actual levels present in the data (clusters, geographies, etc.),
# return a deterministic level → colour mapping. Same levels in two runs →
# same colours, so cross-figure consistency is automatic.
#
# Two palette families:
#   discrete_palette(levels, option = "D")  viridis qualitative
#   highlight_palette(levels, highlight)    grey + accent for the highlight
#
# Source from each plot script:
#   source(file.path(SCRIPT_DIR, "palettes.R"))
# where SCRIPT_DIR is the dirname of the plot script (set at top of script
# from commandArgs(trailingOnly = FALSE)).
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(viridisLite)
})

#' Return a named character vector mapping `levels` → hex colour.
#'
#' @param levels      Character vector of factor levels (or any sortable set).
#' @param option      Viridis option ("D", "C", "B", "A", "E", "F", "G", "H").
#' @param begin,end   Trim the palette range (0..1).
#' @param na_colour   Optional grey for any NA / missing-level entries; if
#'                    supplied, an entry named "" or "NA" is added.
#' @return Named character vector, length = length(levels), names = levels.
discrete_palette <- function(levels,
                             option   = "D",
                             begin    = 0.05,
                             end      = 0.90,
                             na_colour = NULL) {
  levels <- unique(as.character(levels))
  levels <- levels[!is.na(levels) & nzchar(levels)]
  levels <- sort(levels)                # deterministic ordering
  n <- length(levels)
  if (n == 0) return(character(0))
  cols <- viridisLite::viridis(n, option = option, begin = begin, end = end)
  pal <- setNames(cols, levels)
  if (!is.null(na_colour)) {
    pal <- c(pal, setNames(na_colour, "NA"))
  }
  pal
}

#' Grey-and-highlight palette: every level in `levels` gets the same grey
#' EXCEPT `highlight`, which gets a distinguished colour. Useful for
#' "this group versus everyone else" emphasis figures.
#'
#' @param levels      Character vector of levels.
#' @param highlight   The one level to call out (must be in `levels`).
#' @param accent      Hex colour for the highlight.
#' @param base        Hex colour for every other level.
highlight_palette <- function(levels, highlight,
                              accent = "#E69F00", base = "grey75") {
  levels <- unique(as.character(levels))
  levels <- levels[!is.na(levels) & nzchar(levels)]
  pal <- setNames(rep(base, length(levels)), levels)
  if (highlight %in% levels) pal[[highlight]] <- accent
  pal
}

#' Resolve a level vector → colour via the role-default palette. Keeps
#' the role↔palette-option mapping in one place so different scripts can
#' agree on (e.g.) "geography uses option C".
role_palette <- function(levels, role = c("group", "geography", "country")) {
  role <- match.arg(role)
  opt  <- switch(role, group = "D", geography = "C", country = "H")
  discrete_palette(levels, option = opt)
}
