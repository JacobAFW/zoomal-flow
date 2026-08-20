#!/usr/bin/env Rscript
# introgression_detect.R — the isolated, pluggable introgression call (spec §2.3)
# --------------------------------------------------------------------------
# NOT a standalone script: sourced by introgression_pair.R and by the
# rule/threshold sweep (scripts/R/introgression_rule_sweep.R). It exposes ONE
# decision function plus the point-in-polygon primitive it is built on, and
# nothing else — swapping or adding a rule never touches the distance or
# density code in introgression_core.R.
#
#   is_introgressed(points, contours_own, contours_other, params) -> logical
#
# Two rules behind the `detection_rule` switch:
#
#   absolute (default, = V1)   inside the OTHER cluster's contour at density
#                              level >= contour_level_other, AND NOT inside
#                              its OWN cluster's contour at level >=
#                              contour_level_own. Both default to 5e-4,
#                              which is V1 (find_introgression.R:241, 257-259).
#   relative                   the point reaches a HIGHER density layer in the
#                              other cluster than in its own — "more core to
#                              the other than to itself". No absolute cutoff,
#                              so it does not drift when one cluster is more
#                              diffuse than another.
#
# The "point sits in both clusters' peripheries" definition was considered and
# rejected in the spec (§2.3): that is ambiguity, not introgression, and it is
# where hypervariable noise lands.
#
# V1-comparability caveat: V1 wrote `filter(level > 5e-4)` (strict). The
# default contour breaks are `pretty(range(density), 10)`, whose lowest break
# frequently IS 5e-4 exactly — so V1's strict `>` silently dropped the
# outermost ring. The spec asks for `>=`, which is what is implemented here.
# On the same data, `absolute` @ 5e-4 is therefore very slightly more
# permissive than V1 at the same nominal number. Raise contour_level_other to
# 1e-3 to recover V1's effective region when the breaks land on 5e-4.
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sp)
})

#' Deepest density layer of `contours` that contains each point.
#'
#' Contour rings from `stat_contour` are nested (higher level inside lower),
#' so "inside the cluster at level >= L" is equivalent to "inside at least one
#' ring whose level is >= L", and the deepest containing level fully
#' characterises a point's position in the cloud. Each ring is tested
#' separately with `sp::point.in.polygon()` — V1 concatenated every ring's
#' vertices into one call, which is undefined for nested rings.
#'
#' @param points   data.frame/tibble with numeric `x`, `y`.
#' @param contours tibble with `level`, `ring`, `x`, `y` (one closed ring per
#'                 distinct `ring` value). Zero rows = no cloud (returns -Inf).
#' @return Numeric vector, length nrow(points): deepest containing level, or
#'         -Inf for points outside every ring.
contour_max_level <- function(points, contours) {
  out <- rep(-Inf, nrow(points))
  if (is.null(contours) || nrow(contours) == 0 || nrow(points) == 0) return(out)
  for (rg in unique(contours$ring)) {
    ring <- contours[contours$ring == rg, , drop = FALSE]
    lvl  <- ring$level[1]
    # Rings already deeper than the running max cannot improve it.
    if (all(out >= lvl)) next
    pip <- sp::point.in.polygon(
      point.x = points$x, point.y = points$y,
      pol.x   = ring$x,   pol.y   = ring$y
    )
    inside <- pip > 0            # 1 = interior, 2 = edge, 3 = vertex
    out[inside & lvl > out] <- lvl
  }
  out
}

#' Is each sample-window introgressed?
#'
#' Pure: no file I/O, no globals, no cluster names. Given the point(s), the
#' own-cluster contours, the other-cluster contours, and the rule parameters,
#' return one boolean per point.
#'
#' @param points        data.frame with numeric `x`, `y` (the sample-window's
#'                      (distance-to-Kx, distance-to-Ky) coordinates).
#' @param contours_own  Contours of the cluster the sample is ASSIGNED to.
#' @param contours_other Contours of the cluster it might be introgressed from.
#' @param params        list(detection_rule = "absolute"|"relative",
#'                           contour_level_other = 5e-4,
#'                           contour_level_own   = 5e-4)
#' @return list(introgressed = logical, level_own = numeric, level_other = numeric)
is_introgressed <- function(points, contours_own, contours_other, params) {
  rule <- params$detection_rule %||% "absolute"
  if (!rule %in% c("absolute", "relative")) {
    stop("[introgression_detect] unknown detection_rule: ", rule,
         " (expected 'absolute' or 'relative')")
  }

  level_own   <- contour_max_level(points, contours_own)
  level_other <- contour_max_level(points, contours_other)

  introgressed <- if (rule == "absolute") {
    thr_other <- as.numeric(params$contour_level_other)
    thr_own   <- as.numeric(params$contour_level_own)
    (level_other >= thr_other) & !(level_own >= thr_own)
  } else {
    # relative: strictly more core to the other cluster than to its own.
    # -Inf > -Inf is FALSE, so points outside both clouds are never called.
    level_other > level_own
  }

  list(introgressed = introgressed,
       level_own    = level_own,
       level_other  = level_other)
}

# Small null-coalesce so `params` may omit keys the chosen rule doesn't use.
`%||%` <- function(a, b) if (is.null(a)) b else a
