#!/usr/bin/env bash
# migrate_to_sampleset_layout.sh — one-time upgrade of an EXISTING outputs/
# tree to the {sampleset} layout (docs/clonality.md)
# --------------------------------------------------------------------------
# WHY THIS EXISTS. The clonality increment parameterises the four frequency
# stages by a `{sampleset}` wildcard (`full` | `unique`), so their outputs now
# live under `outputs/<stage>/<sampleset>/…`. A cohort that has never been run
# needs nothing from this script — Snakemake just builds into the new layout.
# A cohort that HAS been run already has hours of results (ADMIXTURE alone is
# ~2 h on the Indo cohort) sitting at the old paths, and Snakemake would
# happily recompute every one of them. This moves those results to where the
# new rules expect them so nothing is recomputed.
#
# WHAT MOVES, AND WHAT DOES NOT. Only artifacts that genuinely depend on which
# samples are in the analysis move into `full/`. The upstream structure prep is
# sample-set-INDEPENDENT (it runs before any sample selection) and stays
# exactly where it is — that is the shared seam both sample sets branch from:
#
#   stays at outputs/structure/   snps.*.vcf.gz, chrom_update.txt,
#                                 Pk.{bed,bim,fam,log,nosex}, Pk.smiss,
#                                 Pk.afreq, Pk.dups
#   moves to outputs/structure/full/
#                                 cleaned.*, admixture/, best_k.txt,
#                                 admix_clusters*.tsv, Pk.eigen*, Pk.dist*,
#                                 pca_variance.tsv
#
# Logs move too, because Snakemake treats a rule's `log:` as one of its
# outputs — `admixture_run` declares its log, so `full` and `unique` would
# collide on the same filename if the logs were not namespaced as well.
#
# SAFE TO RE-RUN. Every move is guarded by "does the source still exist", so a
# second invocation is a no-op. Nothing is deleted: a file that is already at
# its destination is left alone, and the script never overwrites.
#
# REVERSIBLE. `--undo` moves everything back.
#
# Usage:
#   bash scripts/sh/migrate_to_sampleset_layout.sh [--outputs DIR] [--logs DIR]
#                                                  [--reports DIR] [--dry-run] [--undo]
# --------------------------------------------------------------------------
set -euo pipefail

OUTPUTS="outputs"
LOGS="logs"
REPORTS="reports"
DRY_RUN=0
UNDO=0

while [ $# -gt 0 ]; do
  case "$1" in
    --outputs) OUTPUTS="$2"; shift 2 ;;
    --logs)    LOGS="$2";    shift 2 ;;
    --reports) REPORTS="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1;    shift ;;
    --undo)    UNDO=1;       shift ;;
    -h|--help) sed -n '2,44p' "$0"; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

SET="full"
moved=0
skipped=0

say() { printf '[migrate] %s\n' "$*"; }

# Move $1 -> $2, creating the parent. No-op if the source is gone or the
# destination already exists (idempotent, never overwrites).
move() {
  local src="$1" dst="$2"
  if [ ! -e "$src" ]; then skipped=$((skipped + 1)); return 0; fi
  if [ -e "$dst" ]; then
    say "SKIP  $src -> $dst  (destination already exists)"
    skipped=$((skipped + 1)); return 0
  fi
  if [ "$DRY_RUN" -eq 1 ]; then
    say "would move  $src -> $dst"
  else
    mkdir -p "$(dirname "$dst")"
    mv "$src" "$dst"
    say "moved $src -> $dst"
  fi
  moved=$((moved + 1))
}

# --------------------------------------------------------------------------
# The manifest: (relative path under its root, root kind)
# --------------------------------------------------------------------------
structure_items=(
  cleaned.bed cleaned.bim cleaned.fam cleaned.log cleaned.nosex
  cleaned.ld.bed cleaned.ld.bim cleaned.ld.fam cleaned.ld.log cleaned.ld.nosex
  cleaned.prune.in cleaned.prune.out
  admixture best_k.txt admix_clusters.tsv admix_clusters_gis.tsv
  Pk.eigenvec Pk.eigenval Pk.dist Pk.dist.id pca_variance.tsv
)
# Stage-3 rules whose logs are per-sample-set.
structure_logs=(
  admixture_cv_plot.log admixture_cv_table.log assign_clusters.log
  distance_matrix.log final_filters.log gis_join.log ld_prune.log
  pca.log pca_variance.log select_best_k.log
  plot_admixture_bars.log plot_admixture_bars_by_country.log
  plot_admixture_bars_by_geography.log plot_pca.log plot_pca_by_geography.log
  plot_nj_tree.log plot_province_map.log
)
moi_items=( fws_by_country.tsv fws_by_geography.tsv )
# Stage-3 and Stage-5 figures are per-sample-set too — `full` and `unique`
# would otherwise write the same PNG. Stage-2 density panels and Stage-4b
# figures are NOT here: they are per-sample distributions / full-set-only.
figure_items=(
  admixture_bars.png admixture_bars.svg
  admixture_bars_by_country.png admixture_bars_by_country.svg
  admixture_bars_by_geography.png admixture_bars_by_geography.svg
  admixture_cv.png admixture_cv.svg
  pca.png pca.svg pca_by_geography.png pca_by_geography.svg
  njt.png njt.svg admix_map.png admix_map.svg
  introgression_shoulder.png introgression_shoulder.svg
)

do_migrate() {
  say "migrating to the {sampleset} layout (sampleset = '$SET')"

  # -- Stage 3 structure ---------------------------------------------------
  for item in "${structure_items[@]}"; do
    move "$OUTPUTS/structure/$item" "$OUTPUTS/structure/$SET/$item"
  done
  for log in "${structure_logs[@]}"; do
    move "$LOGS/structure/$log" "$LOGS/structure/$SET/$log"
  done
  # ADMIXTURE's per-K logs are an output of admixture_run.
  if [ -d "$LOGS/structure" ]; then
    for f in "$LOGS"/structure/admixture_K*.log; do
      [ -e "$f" ] || continue
      move "$f" "$LOGS/structure/$SET/$(basename "$f")"
    done
  fi

  # -- Stage 2 group summaries (per-sample Fws/MOI is NOT sample-set
  #    dependent and deliberately stays put) --------------------------------
  for item in "${moi_items[@]}"; do
    move "$OUTPUTS/moi/$item" "$OUTPUTS/moi/$SET/$item"
  done
  move "$LOGS/moi/fws_summary.log" "$LOGS/moi/$SET/fws_summary.log"

  # -- Stage 3 + Stage 5 figures -------------------------------------------
  for item in "${figure_items[@]}"; do
    move "$REPORTS/figures/$item" "$REPORTS/figures/$SET/$item"
  done
  # The focal-introgression figure is named after the configured focal group.
  if [ -d "$REPORTS/figures" ]; then
    for f in "$REPORTS"/figures/introgression_focal_*.png "$REPORTS"/figures/introgression_focal_*.svg; do
      [ -e "$f" ] || continue
      move "$f" "$REPORTS/figures/$SET/$(basename "$f")"
    done
  fi

  # -- Stage 5 introgression: everything under the stage dir except the
  #    sampleset dirs themselves --------------------------------------------
  if [ -d "$OUTPUTS/introgression" ]; then
    for f in "$OUTPUTS"/introgression/*; do
      [ -e "$f" ] || continue
      base="$(basename "$f")"
      case "$base" in full|unique) continue ;; esac
      move "$f" "$OUTPUTS/introgression/$SET/$base"
    done
  fi
  if [ -d "$LOGS/introgression" ]; then
    for f in "$LOGS"/introgression/*; do
      [ -e "$f" ] || continue
      base="$(basename "$f")"
      case "$base" in full|unique) continue ;; esac
      move "$f" "$LOGS/introgression/$SET/$base"
    done
  fi

  # -- Stage 6 selection: one directory per model ---------------------------
  if [ -d "$OUTPUTS/selection" ]; then
    for f in "$OUTPUTS"/selection/*; do
      [ -e "$f" ] || continue
      base="$(basename "$f")"
      case "$base" in full|unique) continue ;; esac
      move "$f" "$OUTPUTS/selection/$SET/$base"
    done
  fi
  if [ -d "$LOGS/selection" ]; then
    for f in "$LOGS"/selection/*; do
      [ -e "$f" ] || continue
      base="$(basename "$f")"
      case "$base" in full|unique) continue ;; esac
      move "$f" "$LOGS/selection/$SET/$base"
    done
  fi
}

do_undo() {
  say "UNDO: moving the '$SET' layout back to the flat one"
  for item in "${structure_items[@]}"; do
    move "$OUTPUTS/structure/$SET/$item" "$OUTPUTS/structure/$item"
  done
  for log in "${structure_logs[@]}"; do
    move "$LOGS/structure/$SET/$log" "$LOGS/structure/$log"
  done
  if [ -d "$LOGS/structure/$SET" ]; then
    for f in "$LOGS"/structure/"$SET"/admixture_K*.log; do
      [ -e "$f" ] || continue
      move "$f" "$LOGS/structure/$(basename "$f")"
    done
  fi
  for item in "${moi_items[@]}"; do
    move "$OUTPUTS/moi/$SET/$item" "$OUTPUTS/moi/$item"
  done
  for item in "${figure_items[@]}"; do
    move "$REPORTS/figures/$SET/$item" "$REPORTS/figures/$item"
  done
  if [ -d "$REPORTS/figures/$SET" ]; then
    for f in "$REPORTS/figures/$SET"/introgression_focal_*; do
      [ -e "$f" ] || continue
      move "$f" "$REPORTS/figures/$(basename "$f")"
    done
  fi
  move "$LOGS/moi/$SET/fws_summary.log" "$LOGS/moi/fws_summary.log"
  for stage in introgression selection; do
    if [ -d "$OUTPUTS/$stage/$SET" ]; then
      for f in "$OUTPUTS/$stage/$SET"/*; do
        [ -e "$f" ] || continue
        move "$f" "$OUTPUTS/$stage/$(basename "$f")"
      done
    fi
    if [ -d "$LOGS/$stage/$SET" ]; then
      for f in "$LOGS/$stage/$SET"/*; do
        [ -e "$f" ] || continue
        move "$f" "$LOGS/$stage/$(basename "$f")"
      done
    fi
  done
}

if [ "$UNDO" -eq 1 ]; then do_undo; else do_migrate; fi

say "done: $moved moved, $skipped skipped/absent"
if [ "$DRY_RUN" -eq 1 ]; then
  say "(dry run — nothing was changed)"
elif [ "$UNDO" -eq 0 ]; then
  # Snakemake decides staleness by mtime. The moved results keep their
  # original timestamps, which are the SAME MINUTE as some of the upstream
  # seam files they depend on (they were produced in one run). Bump them so
  # every migrated output is strictly newer than the seam it reads, otherwise
  # the first `snakemake` call recomputes the lot.
  say "bumping mtimes so migrated outputs are newer than the shared seam"
  find "$OUTPUTS/structure/$SET" "$OUTPUTS/moi/$SET" \
       "$OUTPUTS/introgression/$SET" "$OUTPUTS/selection/$SET" \
       "$REPORTS/figures/$SET" \
       -exec touch {} + 2>/dev/null || true
  say "next: snakemake --rerun-triggers mtime -n   (expect NO Stage-3 re-runs)"
fi
