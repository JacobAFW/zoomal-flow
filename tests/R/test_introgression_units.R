#!/usr/bin/env Rscript
# test_introgression_units.R — component unit tests for the Stage 5 cores
# --------------------------------------------------------------------------
# Covers, each in isolation (spec §8.2):
#   - window_bin()                 window binning arithmetic
#   - window_id()                  stable, readr-safe window identifier
#   - dominant_alleles()           per-cluster consensus + tie handling
#   - window_distances()           percent-mismatch distance + min-SNP cut
#   - contour_max_level()          point-in-polygon membership on toy contours
#   - is_introgressed()            ALL THREE detection rules
#   - distance_calibration()       the adaptive rule's per-pair quantile inputs
#   - pair_calls()                 both-directions orientation + output schema
#
# Deliberately dependency-free (no testthat in the vvg-box env): a tiny
# ok()/expect_equal() harness, non-zero exit on any failure. Run directly or
# via tests/test_introgression_r_units.py.
#
# Usage:
#   Rscript tests/R/test_introgression_units.R
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
})

script_args <- commandArgs(trailingOnly = FALSE)
SCRIPT_DIR <- dirname(sub("^--file=", "",
                          script_args[grep("^--file=", script_args)][1]))
R_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", "..", "scripts", "R"))
source(file.path(R_DIR, "introgression_core.R"))
source(file.path(R_DIR, "introgression_detect.R"))

# -- micro test harness ----------------------------------------------------
.n_pass <- 0L
.n_fail <- 0L
ok <- function(cond, label) {
  if (isTRUE(cond)) {
    .n_pass <<- .n_pass + 1L
    cat(sprintf("  PASS  %s\n", label))
  } else {
    .n_fail <<- .n_fail + 1L
    cat(sprintf("  FAIL  %s\n", label))
  }
}
eq <- function(actual, expected, label, tol = 1e-9) {
  same <- tryCatch(isTRUE(all.equal(actual, expected, tolerance = tol)),
                   error = function(e) FALSE)
  if (!same) {
    cat(sprintf("        expected: %s\n        actual:   %s\n",
                paste(utils::capture.output(print(expected)), collapse = " | "),
                paste(utils::capture.output(print(actual)),   collapse = " | ")))
  }
  ok(same, label)
}
section <- function(x) cat(sprintf("\n[%s]\n", x))

# --------------------------------------------------------------------------
section("window_bin")
# Midpoint of the containing window: floor(pos/ws)*ws + ws/2.
eq(window_bin(1, 10000),      5000, "position 1 -> first window midpoint")
eq(window_bin(9999, 10000),   5000, "last position of window 1")
eq(window_bin(10000, 10000), 15000, "first position of window 2")
eq(window_bin(45123, 10000), 45000, "mid-genome position")
eq(window_bin(c(1, 10000, 20001), 10000), c(5000, 15000, 25000),
   "vectorised over positions")
eq(window_bin(1500, 1000), 1500, "non-default window size (1 kb)")

section("window_id")
eq(window_id(2, 45000), "w2_45000", "id is w<chrom>_<bin>")
eq(window_id(c(1, 2), c(5000, 45000)), c("w1_5000", "w2_45000"), "vectorised")
# Regression: a bare "1:5000" is parsed as a clock time by readr's guesser,
# which silently rewrote the id to 01:50:00 on round-trip.
ok(!grepl(":", window_id(1, 5000)),
   "id has no colon (readr would parse it as a time)")
eq(window_id(2, 5e4), "w2_50000", "large bins never render in scientific notation")

# --------------------------------------------------------------------------
section("dominant_alleles")
# Two clusters, three positions:
#   pos 1: A is 3x allele 0 vs 1x allele 1 -> 0;  B is 3x allele 1        -> 1
#   pos 2: A ties 2-2                      -> position dropped entirely
#   pos 3: A all 1, B all 1                -> 1 / 1 (no divergence)
gt <- tribble(
  ~SAMPLE, ~CHROM, ~POS, ~SNP, ~Cluster,
  "a1", 1L, 1L, 0L, "A",  "a2", 1L, 1L, 0L, "A",
  "a3", 1L, 1L, 0L, "A",  "a4", 1L, 1L, 1L, "A",
  "b1", 1L, 1L, 1L, "B",  "b2", 1L, 1L, 1L, "B",  "b3", 1L, 1L, 1L, "B",
  "a1", 1L, 2L, 0L, "A",  "a2", 1L, 2L, 0L, "A",
  "a3", 1L, 2L, 1L, "A",  "a4", 1L, 2L, 1L, "A",
  "b1", 1L, 2L, 1L, "B",  "b2", 1L, 2L, 1L, "B",  "b3", 1L, 2L, 1L, "B",
  "a1", 1L, 3L, 1L, "A",  "a2", 1L, 3L, 1L, "A",
  "a3", 1L, 3L, 1L, "A",  "a4", 1L, 3L, 1L, "A",
  "b1", 1L, 3L, 1L, "B",  "b2", 1L, 3L, 1L, "B",  "b3", 1L, 3L, 1L, "B"
)
cons <- dominant_alleles(gt, c("A", "B"))
eq(sort(cons$POS), c(1L, 3L), "tied position is dropped (V1's n < K+1 rule)")
eq(cons$A[cons$POS == 1], 0L, "cluster A dominant allele at pos 1")
eq(cons$B[cons$POS == 1], 1L, "cluster B dominant allele at pos 1")
eq(cons$A[cons$POS == 3], 1L, "no-divergence position keeps both consensuses")

# Missing calls (-1) are ignored, not counted as an allele.
gt_missing <- tribble(
  ~SAMPLE, ~CHROM, ~POS, ~SNP, ~Cluster,
  "a1", 1L, 1L, -1L, "A",  "a2", 1L, 1L, 1L, "A",
  "b1", 1L, 1L,  0L, "B",  "b2", 1L, 1L, 0L, "B"
)
cons_missing <- dominant_alleles(gt_missing, c("A", "B"))
eq(cons_missing$A, 1L, "missing (-1) calls excluded from the consensus vote")

# A cluster with NO non-missing call still gets a column (all NA).
gt_absent <- tribble(
  ~SAMPLE, ~CHROM, ~POS, ~SNP, ~Cluster,
  "a1", 1L, 1L,  0L, "A",  "a2", 1L, 1L, 0L, "A",
  "b1", 1L, 1L, -1L, "B",  "b2", 1L, 1L, -1L, "B"
)
cons_absent <- dominant_alleles(gt_absent, c("A", "B"))
ok("B" %in% names(cons_absent) && all(is.na(cons_absent$B)),
   "cluster with no calls yields an all-NA consensus column")

# --------------------------------------------------------------------------
section("window_distances")
# One window (positions 1..10 -> bin 5000 at ws = 10000), two clusters.
# Consensus: A = 0 everywhere, B = 1 everywhere.
mk_gt <- function(sample, cluster, snps) {
  tibble(SAMPLE = sample, CHROM = 1L, POS = seq_along(snps),
         SNP = as.integer(snps), Cluster = cluster)
}
cons10 <- tibble(CHROM = 1L, POS = 1:10, A = 0L, B = 1L)
# a_pure matches A at all 10 -> DIST_X 0, DIST_Y 100
# a_intro matches B at all 10 -> DIST_X 100, DIST_Y 0  (the introgressed shape)
# a_half  matches A at 7 of 10 -> DIST_X 30, DIST_Y 70
gt10 <- bind_rows(
  mk_gt("a_pure",  "A", rep(0, 10)),
  mk_gt("a_intro", "A", rep(1, 10)),
  mk_gt("a_half",  "A", c(rep(0, 7), rep(1, 3))),
  mk_gt("b_pure",  "B", rep(1, 10))
)
d <- window_distances(gt10, cons10, "A", "B", window_size_bp = 10000, min_snps = 5)
eq(nrow(d), 4L, "one row per (sample, window) clearing min_snps")
eq(d$DIST_X[d$SAMPLE == "a_pure"],    0, "sample matching its own consensus: distance 0")
eq(d$DIST_Y[d$SAMPLE == "a_pure"],  100, "...and 100% from the other consensus")
eq(d$DIST_X[d$SAMPLE == "a_intro"], 100, "consensus-swapped sample: 100% from own")
eq(d$DIST_Y[d$SAMPLE == "a_intro"],   0, "...and 0% from the other")
eq(d$DIST_X[d$SAMPLE == "a_half"],   30, "percent mismatch is a straight proportion")
eq(d$N_SNPS[d$SAMPLE == "a_half"],   10L, "N_SNPS counts non-missing calls used")
eq(unique(d$WINDOW), "w1_5000", "window id derived from CHROM + bin")

# min_snps is a STRICT >: 5 means "at least 6 calls in the window".
d_thin <- window_distances(
  bind_rows(mk_gt("a_thin", "A", rep(0, 5)), mk_gt("b_thin", "B", rep(1, 5))),
  tibble(CHROM = 1L, POS = 1:5, A = 0L, B = 1L),
  "A", "B", window_size_bp = 10000, min_snps = 5)
eq(nrow(d_thin), 0L, "window with exactly min_snps calls is dropped (strict >)")

# Equidistant sample-windows are dropped (V1's DIST_X != DIST_Y filter).
cons_tie <- tibble(CHROM = 1L, POS = 1:10, A = 0L, B = 0L)
d_tie <- window_distances(
  bind_rows(mk_gt("a_tie", "A", c(rep(0, 5), rep(1, 5))),
            mk_gt("b_tie", "B", c(rep(0, 5), rep(1, 5)))),
  cons_tie, "A", "B", window_size_bp = 10000, min_snps = 5)
eq(nrow(d_tie), 0L, "equidistant sample-windows are dropped")

# Missing calls are excluded from both numerator and denominator.
d_miss <- window_distances(
  bind_rows(mk_gt("a_miss", "A", c(rep(-1, 2), rep(0, 8))),
            mk_gt("b_miss", "B", rep(1, 10))),
  cons10, "A", "B", window_size_bp = 10000, min_snps = 5)
eq(d_miss$N_SNPS[d_miss$SAMPLE == "a_miss"], 8L, "missing calls drop out of N_SNPS")
eq(d_miss$DIST_X[d_miss$SAMPLE == "a_miss"], 0, "...and out of the distance denominator")

# --------------------------------------------------------------------------
section("contour_max_level")
# Toy contours: two concentric squares per cluster.
#   cluster A: level 1e-4 ring spans [0,40]^2, level 1e-3 ring spans [10,30]^2
#   cluster B: level 1e-4 ring spans [60,100]^2, level 1e-3 ring [70,90]^2
square <- function(cluster, level, ring, lo, hi) {
  tibble(cluster = cluster, level = level, ring = ring,
         x = c(lo, hi, hi, lo, lo), y = c(lo, lo, hi, hi, lo))
}
cont_a <- bind_rows(square("A", 1e-4, "a-outer",  0, 40),
                    square("A", 1e-3, "a-inner", 10, 30))
cont_b <- bind_rows(square("B", 1e-4, "b-outer", 60, 100),
                    square("B", 1e-3, "b-inner", 70,  90))

pts <- tibble(x = c(20,  5, 50, 80,  65),
              y = c(20,  5, 50, 80,  65))
lv_a <- contour_max_level(pts, cont_a)
lv_b <- contour_max_level(pts, cont_b)
eq(lv_a, c(1e-3, 1e-4, -Inf, -Inf, -Inf), "deepest containing level, cluster A")
eq(lv_b, c(-Inf, -Inf, -Inf, 1e-3, 1e-4), "deepest containing level, cluster B")
eq(contour_max_level(pts, cont_a[0, ]), rep(-Inf, 5),
   "no contours -> every point outside")
eq(contour_max_level(pts[0, ], cont_a), numeric(0), "no points -> empty result")
ok(contour_max_level(tibble(x = 0, y = 20), cont_a) == 1e-4,
   "a point on the ring boundary counts as inside")

# --------------------------------------------------------------------------
section("is_introgressed — absolute rule")
abs_params <- list(detection_rule = "absolute",
                   contour_level_other = 5e-4, contour_level_own = 5e-4)
# An A-assigned sample-window sitting deep in B's core, outside A's cloud.
r <- is_introgressed(tibble(x = 80, y = 80), cont_a, cont_b, abs_params)
ok(r$introgressed, "deep in other + outside own -> introgressed")
eq(r$level_other, 1e-3, "reports the other-cluster level it reached")
eq(r$level_own,  -Inf,  "reports the own-cluster level it reached")

# Only in B's shallow layer (1e-4 < 5e-4): not deep enough.
ok(!is_introgressed(tibble(x = 65, y = 65), cont_a, cont_b, abs_params)$introgressed,
   "other-cluster level below contour_level_other -> not called")
# In its own core: rejected even though it is nowhere near the other cluster.
ok(!is_introgressed(tibble(x = 20, y = 20), cont_a, cont_b, abs_params)$introgressed,
   "inside its own core -> not called")
# Outside everything.
ok(!is_introgressed(tibble(x = 50, y = 50), cont_a, cont_b, abs_params)$introgressed,
   "outside both clouds -> not called")
# Vectorised over several points at once.
eq(is_introgressed(tibble(x = c(80, 20, 50), y = c(80, 20, 50)),
                   cont_a, cont_b, abs_params)$introgressed,
   c(TRUE, FALSE, FALSE), "vectorised over points")

# Lowering contour_level_other admits the shallow layer.
loose <- modifyList(abs_params, list(contour_level_other = 1e-4))
ok(is_introgressed(tibble(x = 65, y = 65), cont_a, cont_b, loose)$introgressed,
   "lowering contour_level_other admits a shallower hit")
# Threshold comparison is >= (spec), not V1's strict >.
edge <- modifyList(abs_params, list(contour_level_other = 1e-3))
ok(is_introgressed(tibble(x = 80, y = 80), cont_a, cont_b, edge)$introgressed,
   "level exactly equal to the threshold counts (>=, per spec §2.3)")

section("is_introgressed — relative rule")
rel_params <- list(detection_rule = "relative")
ok(is_introgressed(tibble(x = 65, y = 65), cont_a, cont_b, rel_params)$introgressed,
   "shallow in other, outside own -> introgressed (no absolute cutoff)")
ok(!is_introgressed(tibble(x = 20, y = 20), cont_a, cont_b, rel_params)$introgressed,
   "deep in own, outside other -> not called")
ok(!is_introgressed(tibble(x = 50, y = 50), cont_a, cont_b, rel_params)$introgressed,
   "outside both -> not called (-Inf is not > -Inf)")
# Equal depth in both clouds is ambiguity, not introgression.
cont_overlap <- bind_rows(square("B", 1e-4,  "b-ov",  0, 40),
                          square("B", 1e-3, "b-ov2", 10, 30))
ok(!is_introgressed(tibble(x = 20, y = 20), cont_a, cont_overlap, rel_params)$introgressed,
   "equally core to both clusters -> not called (strict >)")


# --------------------------------------------------------------------------
section("is_introgressed — distance rule (deterministic, fixed margin)")
# The distance rule reads d_own / d_other only; the contours are ignored
# entirely (passed as NULL below to prove no density surface is consulted).
dist_params <- list(detection_rule = "distance", distance_margin = 15)

# Clearly closer to the OTHER consensus: 30% from own, 8% from other.
ok(is_introgressed(tibble(x = 0, y = 0, d_own = 30, d_other = 8),
                   NULL, NULL, dist_params)$introgressed,
   "clearly closer to the other consensus -> introgressed")
# A core-own point: 5% from own, 25% from other.
ok(!is_introgressed(tibble(x = 0, y = 0, d_own = 5, d_other = 25),
                    NULL, NULL, dist_params)$introgressed,
   "core-own point -> not called")
# Margin boundary: difference exactly the margin counts (>=), one point under
# does not.
eq(is_introgressed(tibble(x = 0, y = 0,
                          d_own   = c(30,   29.9, 30.1),
                          d_other = c(15,   15,   15)),
                   NULL, NULL, dist_params)$introgressed,
   c(TRUE, FALSE, TRUE), "margin boundary is inclusive (>=), one pp under is not")
# Raising the margin retracts a marginal call; lowering it admits one.
ok(!is_introgressed(tibble(x = 0, y = 0, d_own = 30, d_other = 8), NULL, NULL,
                    modifyList(dist_params, list(distance_margin = 25)))$introgressed,
   "raising distance_margin retracts a marginal call")
ok(is_introgressed(tibble(x = 0, y = 0, d_own = 12, d_other = 8), NULL, NULL,
                   modifyList(dist_params, list(distance_margin = 3)))$introgressed,
   "lowering distance_margin admits a shallower call")
# No density surface is built, so the level columns are NA by contract.
r_dist <- is_introgressed(tibble(x = 0, y = 0, d_own = 30, d_other = 8),
                          cont_a, cont_b, dist_params)
ok(is.na(r_dist$level_own) && is.na(r_dist$level_other),
   "distance rule reports NA density levels (it fits no surface)")
# Missing distances are a caller error, not a silent FALSE.
ok(inherits(tryCatch(is_introgressed(tibble(x = 1, y = 1), NULL, NULL, dist_params),
                     error = function(e) e), "error"),
   "distance rule without d_own/d_other raises an error")
# Determinism: identical input, identical output, no state carried between calls.
eq(is_introgressed(tibble(x = 0, y = 0, d_own = 30, d_other = 8),
                   NULL, NULL, dist_params)$introgressed,
   is_introgressed(tibble(x = 0, y = 0, d_own = 30, d_other = 8),
                   NULL, NULL, dist_params)$introgressed,
   "same input -> same output")

section("is_introgressed — distance rule (adaptive margin)")
# Cluster A sits 2-10% from its own consensus, cluster B 3-11% from its own.
cal <- list(own_self = seq(2, 10, length.out = 50),
            other_self = seq(3, 11, length.out = 50))
ad_params <- list(detection_rule = "distance", distance_adaptive = TRUE,
                  distance_adaptive_quantile = 0.9, distance_calibration = cal)
# q(0.9) is 9.2 for own and 10.2 for other.
ok(is_introgressed(tibble(x = 0, y = 0, d_own = 20, d_other = 6),
                   NULL, NULL, ad_params)$introgressed,
   "atypically far from own + typical for other -> introgressed")
ok(!is_introgressed(tibble(x = 0, y = 0, d_own = 4, d_other = 6),
                    NULL, NULL, ad_params)$introgressed,
   "typical distance from its own cluster -> not called")
ok(!is_introgressed(tibble(x = 0, y = 0, d_own = 20, d_other = 18),
                    NULL, NULL, ad_params)$introgressed,
   "far from own but also atypical for the other -> not called")
ok(!is_introgressed(tibble(x = 0, y = 0, d_own = 5, d_other = 6),
                    NULL, NULL, ad_params)$introgressed,
   "closer to its own than to the other -> not called")
# A tighter own cluster lowers the bar for "atypical".
tight <- modifyList(ad_params, list(
  distance_calibration = list(own_self = seq(1, 2, length.out = 50),
                              other_self = cal$other_self)))
ok(is_introgressed(tibble(x = 0, y = 0, d_own = 6, d_other = 5),
                   NULL, NULL, tight)$introgressed,
   "a tighter own cluster makes a smaller departure atypical (self-calibrating)")
ok(inherits(tryCatch(is_introgressed(tibble(x = 0, y = 0, d_own = 20, d_other = 6),
                                     NULL, NULL,
                                     list(detection_rule = "distance",
                                          distance_adaptive = TRUE)),
                     error = function(e) e), "error"),
   "adaptive mode without a calibration raises an error")

section("distance_calibration")
cal_tbl <- tibble(SAMPLE  = c("a1", "a2", "b1", "b2"),
                  Cluster = c("A", "A", "B", "B"),
                  DIST_X  = c(2, 4, 20, 22),
                  DIST_Y  = c(20, 22, 3, 5))
cc <- distance_calibration(cal_tbl, own = "A", other = "B", kx = "A", ky = "B")
eq(cc$own_self,   c(2, 4), "own_self = A's distance to A (DIST_X)")
eq(cc$other_self, c(3, 5), "other_self = B's distance to B (DIST_Y)")
cc2 <- distance_calibration(cal_tbl, own = "B", other = "A", kx = "A", ky = "B")
eq(cc2$own_self,   c(3, 5), "orientation flips with the direction")
eq(cc2$other_self, c(2, 4), "orientation flips with the direction (other)")

section("pair_calls — orientation and schema")
pc_tbl <- tibble(
  SAMPLE  = c("a1", "a2", "b1"),
  Cluster = c("A",  "A",  "B"),
  CHROM   = c(1L, 1L, 1L),
  BIN     = c(5000, 15000, 5000),
  WINDOW  = c("w1_5000", "w1_15000", "w1_5000"),
  N_SNPS  = c(10L, 10L, 10L),
  DIST_X  = c(2,  30,  25),      # a2 looks like B; b1 looks like A
  DIST_Y  = c(28,  4,   3)
)
pc <- pair_calls(pc_tbl, "A", "B", contours = NULL,
                 params = list(detection_rule = "distance", distance_margin = 15))
eq(sort(pc$SAMPLE), c("a2"), "only the A sample that looks like B is called")
eq(pc$DIRECTION, "A_like_B", "direction is labelled own_like_other")
eq(pc$DIST_OWN,   30, "DIST_OWN is the distance to the sample's OWN cluster")
eq(pc$DIST_OTHER,  4, "DIST_OTHER is the distance to the other cluster")
eq(names(pc), PAIR_CALL_COLS, "output schema matches PAIR_CALL_COLS")
# Flip the injected signal onto the B sample and confirm the reverse direction.
pc_rev <- pc_tbl; pc_rev$DIST_X[3] <- 3; pc_rev$DIST_Y[3] <- 25
pc2 <- pair_calls(pc_rev, "A", "B", contours = NULL,
                  params = list(detection_rule = "distance", distance_margin = 15))
eq(sort(pc2$DIRECTION), c("A_like_B", "B_like_A"),
   "both directions are called from one table")
eq(pc2$DIST_OWN[pc2$SAMPLE == "b1"], 25,
   "the reverse direction reads DIST_Y as its own distance")

# --------------------------------------------------------------------------
section("is_introgressed — errors")
ok(inherits(tryCatch(is_introgressed(tibble(x = 1, y = 1), cont_a, cont_b,
                                     list(detection_rule = "nonsense")),
                     error = function(e) e), "error"),
   "unknown detection_rule raises an error")

# --------------------------------------------------------------------------
cat(sprintf("\n%d passed, %d failed\n", .n_pass, .n_fail))
if (.n_fail > 0) quit(status = 1)
