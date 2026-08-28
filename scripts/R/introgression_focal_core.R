#!/usr/bin/env Rscript
# introgression_focal_core.R — the focal-enrichment statistic (spec §9.6)
# --------------------------------------------------------------------------
# NOT a standalone script: sourced by introgression_focal_test.R and by
# tests/R/test_introgression_units.R. It holds one function, isolated from all
# I/O so the statistic can be unit-tested against hand-computable cases.
#
# Unlike introgression_core.R, nothing here is ported from V1 — V1 had no
# focal-level test at all, only a set difference ("windows unique to Aceh")
# with no null. This is new machinery, so it gets its own file rather than
# being mixed into the preserved cores.
# --------------------------------------------------------------------------

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' Per-window enrichment of a focal subgroup inside its own cluster.
#'
#' THE NULL. Hold each window's set of introgressed samples fixed, and randomly
#' reassign which `n_focal` of the cluster's members carry the focal label.
#' Size-preserving: every replicate has exactly `n_focal` focal members, so the
#' comparison is always "this subgroup vs a random subgroup of the same size".
#'
#' TWO FORMS OF THE SAME NULL. Permuting the label with the called set held
#' fixed makes each window's marginal null exactly hypergeometric:
#'
#'   focal_support ~ Hypergeometric(N = cluster_n, K = n_called, k = n_focal)
#'
#' so `p_exact` (the phyper tail) and `p_perm` (the Monte-Carlo tail) describe
#' the SAME distribution — the first evaluated, the second sampled. Both are
#' returned. `p_exact` is the one to adjust and call on: `p_perm` has a
#' resolution floor of 1/(nperm + 1), which after BH across several hundred
#' windows can make significance unreachable no matter how strong the signal.
#' `p_perm` earns its keep as a live check that the analytic form is right —
#' if the two diverge, something about the null is not what it claims to be.
#'
#' The Monte Carlo shares one focal-label draw across all windows per replicate
#' (that is what "a random subgroup" means, and it is what induces the
#' between-window correlation a per-window analytic form cannot see). The
#' per-window marginal is unaffected, which is why the two agree.
#'
#' @param M          Integer/logical incidence matrix, windows x cluster
#'                   members: 1 = that member is called introgressed at that
#'                   window. Row names are window ids; column order defines the
#'                   member indexing used by `focal_cols`.
#' @param focal_cols Integer column indices of the focal members.
#' @param nperm      Monte-Carlo replicates.
#' @param seed       RNG seed, so the run reproduces.
#' @return tibble(WINDOW, n_called_cluster, focal_support, background_support,
#'                p_perm, p_exact), one row per window, in the row order of `M`.
focal_enrichment <- function(M, focal_cols, nperm = 1000L, seed = 20260828L) {
  stopifnot(is.matrix(M), length(focal_cols) > 0)
  n_cluster <- ncol(M)
  n_focal   <- length(focal_cols)
  if (n_focal >= n_cluster) {
    stop("[focal_enrichment] the focal group is the whole cluster (", n_focal,
         " of ", n_cluster, ") — there is nothing to contrast against")
  }
  if (any(!focal_cols %in% seq_len(n_cluster))) {
    stop("[focal_enrichment] focal_cols out of range for a ", n_cluster,
         "-column matrix")
  }
  storage.mode(M) <- "integer"

  # unname(): rowSums() carries the window ids through as vector names, which
  # would ride into the tibble's numeric columns and out into the TSV header.
  n_called  <- unname(rowSums(M))
  focal_sup <- unname(rowSums(M[, focal_cols, drop = FALSE]))

  # One matrix product runs every replicate: column b of Fmat is one random
  # size-n_focal subgroup, so (M %*% Fmat)[w, b] is that subgroup's support at
  # window w.
  set.seed(seed)
  Fmat <- matrix(0L, nrow = n_cluster, ncol = nperm)
  for (b in seq_len(nperm)) Fmat[sample.int(n_cluster, n_focal), b] <- 1L
  null_sup <- M %*% Fmat

  tibble::tibble(
    WINDOW             = rownames(M) %||% as.character(seq_len(nrow(M))),
    n_called_cluster   = as.integer(n_called),
    focal_support      = as.integer(focal_sup),
    background_support = as.integer(n_called - focal_sup),
    # (b + 1) / (B + 1): the standard estimator, never returns an impossible 0.
    p_perm  = unname(1 + rowSums(null_sup >= focal_sup)) / (nperm + 1),
    p_exact = stats::phyper(focal_sup - 1L, n_called, n_cluster - n_called,
                            n_focal, lower.tail = FALSE))
}
