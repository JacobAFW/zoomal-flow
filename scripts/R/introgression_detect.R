#!/usr/bin/env Rscript
# introgression_detect.R — the isolated, pluggable introgression call (spec §2.3)
# --------------------------------------------------------------------------
# NOT a standalone script: sourced by introgression_pair.R, by the
# rule/threshold sweep (scripts/R/introgression_rule_sweep.R) and by the Malay
# benchmark (scripts/R/introgression_benchmark_malay.R). It exposes the
# detection layer and nothing else — swapping or adding a rule never touches
# the distance or density code in introgression_core.R:
#
#   contour_max_level(points, contours)                  the PIP primitive
#   is_introgressed(points, contours_own, contours_other, params) -> logical
#   distance_calibration(dist_tbl, own, other, kx, ky)   adaptive-rule calibration
#   pair_calls(dist_tbl, kx, ky, contours, params)       both directions -> call table
#
# Three rules behind the `detection_rule` switch:
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
#   distance                   DETERMINISTIC. Decides purely from the two raw
#                              per-window match-fraction distances (d_own,
#                              d_other). Never builds a density surface and
#                              never calls contour_max_level(), so its answer
#                              cannot move when cohort composition changes.
#
# Why `distance` exists (spec §"Method status"). Both density rules inherit
# the density surface's cohort dependence, and two findings show that surface
# is not stable: the bug-fixed density method reproduces only 41/217 (19%) of
# its own published Malay windows once cohort composition shifts, and
# introgression rate scales INVERSELY with cluster size (a small cluster's low,
# broad cloud lets its own members fall below their own contour and satisfy the
# "not in its own cloud" half of the rule). Sample-level answers survive both;
# window-level answers do not. A rule that thresholds the distances directly
# has no fitted surface to drift.
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
  library(dplyr)
  library(tibble)
})

DETECTION_RULES <- c("absolute", "relative", "distance")

#' Rules that need a fitted 2D density surface.
#'
#' Callers use this to skip the (expensive, and for `distance` meaningless)
#' contour build entirely rather than building a surface nothing reads.
rule_needs_density <- function(rule) rule %in% c("absolute", "relative")

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

#' Per-pair calibration for the ADAPTIVE distance rule.
#'
#' The adaptive rule asks two questions that only make sense relative to how
#' far a cluster's members normally sit from their OWN consensus: "is this
#' window atypically far from its own cluster?" and "is it as close to the
#' other cluster as that cluster's own members are?". Both need the two
#' clusters' self-distance distributions, which the per-direction point subset
#' handed to `is_introgressed()` does not contain — so they are computed once
#' per pair from the full distance table and passed in through `params`.
#'
#' Deterministic: same distance table in, same numbers out. No fitting.
#'
#' @param dist_tbl Output of `window_distances()` (SAMPLE, Cluster, DIST_X, DIST_Y).
#' @param own,other The two cluster labels, oriented for one direction.
#' @param kx,ky    The pair's axis labels (DIST_X is distance to kx).
#' @return list(own_self = numeric, other_self = numeric)
distance_calibration <- function(dist_tbl, own, other, kx, ky) {
  self_of <- function(cl) {
    idx <- which(dist_tbl$Cluster == cl)
    if (identical(cl, kx)) dist_tbl$DIST_X[idx] else dist_tbl$DIST_Y[idx]
  }
  list(own_self = self_of(own), other_self = self_of(other))
}

#' The deterministic distance rule (spec §2.3, `detection_rule: "distance"`).
#'
#' Two margin modes, both same-input-same-output:
#'
#'   fixed (default)  introgressed when the window matches the other cluster's
#'                    consensus at least `distance_margin` percentage points
#'                    better than its own: `d_own - d_other >= margin`.
#'   adaptive         `distance_adaptive: true`. The margin comes from the two
#'                    clusters' own within-cluster spread instead of a config
#'                    number: introgressed when the window is atypically far
#'                    from its own consensus (`d_own >=` the q-quantile of the
#'                    own cluster's self-distances) AND no further from the
#'                    other consensus than that cluster's own members typically
#'                    are (`d_other <=` the q-quantile of the other cluster's
#'                    self-distances) AND still strictly closer to the other
#'                    than to its own. `q` = `distance_adaptive_quantile`
#'                    (default 0.9). The third clause costs nothing and rules
#'                    out the degenerate case where both quantile tests pass on
#'                    a window that is still marginally closer to its own.
#'
#' @param points data.frame with numeric `d_own`, `d_other` (percentages).
#' @param params see is_introgressed().
#' @return logical vector, length nrow(points).
detect_by_distance <- function(points, params) {
  if (!all(c("d_own", "d_other") %in% names(points))) {
    stop("[introgression_detect] detection_rule 'distance' needs `d_own` and ",
         "`d_other` columns on `points` (the per-window distances oriented by ",
         "the sample's cluster assignment); got: ",
         paste(names(points), collapse = ", "))
  }
  d_own   <- as.numeric(points$d_own)
  d_other <- as.numeric(points$d_other)

  adaptive <- isTRUE(as.logical(params$distance_adaptive %||% FALSE))
  hit <- if (!adaptive) {
    margin <- as.numeric(params$distance_margin %||% 15)
    if (is.na(margin)) stop("[introgression_detect] distance_margin is not numeric")
    (d_own - d_other) >= margin
  } else {
    cal <- params$distance_calibration
    if (is.null(cal) || is.null(cal$own_self) || is.null(cal$other_self)) {
      stop("[introgression_detect] distance_adaptive needs `distance_calibration` ",
           "in params (list(own_self=, other_self=)) — build it with ",
           "distance_calibration(); pair_calls() does this for you")
    }
    q <- as.numeric(params$distance_adaptive_quantile %||% 0.9)
    if (is.na(q) || q < 0 || q > 1) {
      stop("[introgression_detect] distance_adaptive_quantile must be in [0, 1]")
    }
    own_hi   <- stats::quantile(cal$own_self,   probs = q, na.rm = TRUE, names = FALSE)
    other_hi <- stats::quantile(cal$other_self, probs = q, na.rm = TRUE, names = FALSE)
    (d_own >= own_hi) & (d_other <= other_hi) & (d_own > d_other)
  }

  hit[is.na(hit)] <- FALSE
  hit
}

#' Is each sample-window introgressed?
#'
#' Pure: no file I/O, no globals, no cluster names. Given the point(s), the
#' own-cluster contours, the other-cluster contours, and the rule parameters,
#' return one boolean per point.
#'
#' @param points        data.frame with numeric `x`, `y` (the sample-window's
#'                      (distance-to-Kx, distance-to-Ky) coordinates) and, for
#'                      the `distance` rule, `d_own` / `d_other` — the SAME two
#'                      numbers relabelled by the sample's cluster assignment.
#'                      Callers build all four columns; density rules read
#'                      x/y, the distance rule reads d_own/d_other.
#' @param contours_own  Contours of the cluster the sample is ASSIGNED to.
#'                      Ignored (and may be NULL) under `distance`.
#' @param contours_other Contours of the cluster it might be introgressed from.
#'                      Ignored (and may be NULL) under `distance`.
#' @param params        list(detection_rule = "absolute"|"relative"|"distance",
#'                           contour_level_other = 5e-4,   # absolute
#'                           contour_level_own   = 5e-4,   # absolute
#'                           distance_margin     = 15,     # distance, fixed
#'                           distance_adaptive   = FALSE,  # distance
#'                           distance_adaptive_quantile = 0.9,
#'                           distance_calibration = list(own_self=, other_self=))
#' @return list(introgressed = logical, level_own = numeric, level_other = numeric).
#'         Under `distance` the two level vectors are NA — no density surface
#'         is built, which is the entire point of the rule.
is_introgressed <- function(points, contours_own, contours_other, params) {
  rule <- params$detection_rule %||% "absolute"
  if (!rule %in% DETECTION_RULES) {
    stop("[introgression_detect] unknown detection_rule: ", rule,
         " (expected one of ", paste(DETECTION_RULES, collapse = ", "), ")")
  }

  if (rule == "distance") {
    return(list(introgressed = detect_by_distance(points, params),
                level_own    = rep(NA_real_, nrow(points)),
                level_other  = rep(NA_real_, nrow(points))))
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

# Output schema of a per-pair call file (introgression_aggregate.R reads it).
PAIR_CALL_COLS <- c("PAIR", "SAMPLE", "Cluster", "OTHER", "DIRECTION",
                    "CHROM", "BIN", "WINDOW", "N_SNPS",
                    "DIST_OWN", "DIST_OTHER", "LEVEL_OWN", "LEVEL_OTHER")

#' Call introgression in BOTH directions for one cluster pair.
#'
#' Direction 1: a Kx sample whose window looks like Ky ("Kx_like_Ky").
#' Direction 2: a Ky sample whose window looks like Kx.
#'
#' Shared by introgression_pair.R (the pipeline) and
#' introgression_benchmark_malay.R (the method scoreboard), so both exercise
#' the identical orientation, calibration and output schema.
#'
#' `DIST_OWN` / `DIST_OTHER` are the pair's two axes relabelled relative to the
#' sample's own cluster, so both directions share one schema — and they are the
#' same two numbers the `distance` rule decides on.
#'
#' @param dist_tbl Output of `window_distances()`.
#' @param kx,ky    Cluster labels; DIST_X is the distance to kx.
#' @param contours Named list with `contours[[kx]]` / `contours[[ky]]`, or NULL
#'                 for rules that need no density surface.
#' @param params   Detection params (see is_introgressed). The adaptive
#'                 distance rule's per-pair calibration is derived here.
#' @param pair_id  Value written into the PAIR column.
#' @return Tibble with PAIR_CALL_COLS, hits only (0 rows is a valid answer).
pair_calls <- function(dist_tbl, kx, ky, contours, params,
                       pair_id = paste0(kx, "__", ky)) {
  one_direction <- function(own, other) {
    idx <- which(dist_tbl$Cluster == own)
    if (length(idx) == 0) return(NULL)
    sub <- dist_tbl[idx, , drop = FALSE]
    d_own   <- if (identical(own, kx)) sub$DIST_X else sub$DIST_Y
    d_other <- if (identical(own, kx)) sub$DIST_Y else sub$DIST_X
    pts <- tibble::tibble(x = sub$DIST_X, y = sub$DIST_Y,
                          d_own = d_own, d_other = d_other)

    p <- params
    if (identical(p$detection_rule, "distance") &&
        isTRUE(as.logical(p$distance_adaptive %||% FALSE))) {
      p$distance_calibration <- distance_calibration(dist_tbl, own, other, kx, ky)
    }
    res <- is_introgressed(pts, contours[[own]], contours[[other]], p)

    sub %>%
      dplyr::mutate(
        PAIR        = pair_id,
        OTHER       = other,
        DIRECTION   = paste0(own, "_like_", other),
        DIST_OWN    = d_own,
        DIST_OTHER  = d_other,
        LEVEL_OWN   = res$level_own,
        LEVEL_OTHER = res$level_other,
        HIT         = res$introgressed
      ) %>%
      dplyr::filter(HIT) %>%
      dplyr::select(dplyr::all_of(PAIR_CALL_COLS))
  }
  dplyr::bind_rows(one_direction(kx, ky), one_direction(ky, kx))
}

# Small null-coalesce so `params` may omit keys the chosen rule doesn't use.
`%||%` <- function(a, b) if (is.null(a)) b else a
