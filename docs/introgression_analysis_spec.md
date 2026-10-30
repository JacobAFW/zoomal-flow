# ZOOMAL-Flow — Stage 5: Introgression detection

**Analysis specification & design.** Written 2026-08-17. This is the committed,
user-facing spec for the introgression stage *and* the design record for the V1→agnostic
revamp. It doubles as the long-form documentation referenced from the script header and the
GitHub README.

> **Validation status — read first.** Unlike Stages 0–4, the V1 introgression stage was
> **never validated against an HPC/reference output** (the original `find_introgression.R`
> header states reference outputs were unavailable and "validation diffs are pending"). There
> is therefore *no* bit-for-bit ground truth to reproduce. Correctness here rests on the
> **synthetic positive-control tests** (below), not on matching V1. Treat real-cohort outputs
> as a method result to be reviewed, not a validated number.

> **⚠ Method status — updated 2026-08-20 (post Stage-5b benchmark). SUPERSEDES both the
> "preserve the density method" framing below AND the 2026-08-17 status block, whose two
> motivating findings have since been shown to be measurement artifacts.**
>
> A third, **deterministic** `detection_rule: "distance"` now exists (§2.3): it thresholds the
> raw per-window match-fraction distances and fits no kernel surface. It was added because two
> findings appeared to show the density surface was cohort-unstable. Benchmarking all three
> rules against the published Malaysian result (`make benchmark-malay`) **did not support
> either finding, and did not support switching the default**:
>
> 1. **The "19% window reproduction" was a WINDOW-id artifact.** The 41/217 figure compared
>    `WINDOW` *ids*, which are row numbers over the sorted (CHROM, bin) pairs and therefore
>    shift whenever the variant set shifts — of 138 coordinate-matched windows shared between
>    the published truth and the bug-fixed rerun, exactly **1** carries the same id. Compared
>    the correct way (CHROM + window start), the density method reproduces **138/217 = 64%**
>    of its own published windows, not 19%. This is precisely the failure mode the agnostic
>    port already fixed with a coordinate-derived `window_id()`.
> 2. **The "inverse cluster-size bias" is mostly the per-cluster FILTER, not the rule.**
>    `per_cluster_min_pct` is a *percentage* of cluster n, so its threshold scales with cluster
>    size — 5% of a 410-sample cluster is 21 samples, but only 2 of a 35-sample one.
>    Swapping it for an equal absolute floor (`per_cluster_min_pct: 0`,
>    `per_cluster_min_samples: 6`) and changing nothing else collapses the effect on both
>    cohorts:
>
>    | cohort | per-cluster filter | median windows/sample (large → small cluster) | max/min |
>    |---|---|---|---|
>    | Malay, `absolute` | 5% of cluster n | Mf 10 · Mn 23 · Pen 36 | 3.6 |
>    | Malay, `absolute` | absolute floor n ≥ 6 | Mf 19 · Mn 18 · Pen 19.5 | **1.08** |
>
>    Under an equal floor the density `absolute` rule is essentially **flat** by cluster size.
>
> **Consequence.** The shipped default stays **`absolute` @ 5e-4**. `distance` is available,
> tested and wired, but on the Malay benchmark it wins none of the four measures — it is the
> *least* stable under perturbation (Jaccard 0.65 / 0.50 at a 10% / 20% sample drop, vs
> 0.66 / 0.61 for `absolute`) and shows the *largest* residual size dependence under the equal
> floor (max/min 6.0). Removing the density surface does not buy reproducibility, because the
> **cluster consensus alleles are themselves cohort-derived** — that dependence survives every
> rule. The open question this leaves is a filter question, not a rule question.
>
> **Decision (Jacob, 2026-08-20): the shipped default per-cluster filter becomes the absolute
> floor** — `per_cluster_min_samples` with `per_cluster_min_pct: 0` — because the headline
> claims are cross-cluster comparisons where the size confound must go.
> The percentage stays available for within-cluster framing. **The floor `N` is derived, not
> guessed:** set it by a null-permutation / FDR argument — permute the introgression calls
> within cluster, recompute per-window sample-support under the null, take the smallest `N`
> that essentially no null window reaches, and report the null distribution plus a stability
> sweep around it. An FDR-safe floor is cohort-specific, so `config/config.yaml` ships the
> mechanism and a generic minimum rather than a number, and points at the derivation script.
> This also closes the diagnosis's flagged "no FDR anywhere" gap.
>
> **Cluster-level and focal-level are two different tests, and the floor belongs to the first.**
> An FDR-safe floor scales with cluster size, so a single global floor is set by the largest
> cluster and will exclude the smallest. A focal subgroup therefore cannot clear it: a subgroup
> of *n* samples can never reach support above *n*, so any floor above that makes a
> focal-unique window structurally unreachable. This is a mis-posed test, not a threshold to
> tune. The cluster floor stays as derived, and the focal group gets **its own test at its own
> scale** — a focal-vs-rest contrast inside the parent cluster, with a subgroup-scaled
> permutation null and BH control (§6). V1's "windows unique to <focal>" output is retained as
> **descriptive only**.
>
> **Three V1 defects found + fixed during the port** (keep): (1) the GFF↔reference contig map
> covered only 3 of 14 chromosomes (exact-length matching broke on the A1.H.1 Icor
> corrections) → now exact-name/exact-length/nearest-length-within-1%, logged; masked windows
> 87→384. (2) `point.in.polygon` was run over concatenated nested rings (undefined) → now
> per-ring, deepest containing level. (3) `level > 5e-4` dropped the outermost ring when breaks
> landed on 5e-4 → implemented as `>=` per spec.

---

## 1. The question

For each genomic window, does a sample carry the genetic signature of a cluster *other than
the one it's assigned to* — evidence of introgression (gene flow) between clusters? The
headline product is windows introgressed **uniquely in a focal group** relative
to the rest of its cluster.

## 2. The method (unchanged in spirit, generalised to pairs)

The V1 method — which we are deliberately keeping, as it has performed well and had good
feedback — works like this:

1. **Cluster consensus.** At each SNP, compute each cluster's dominant (most common) allele.
2. **Per-window distance.** For each sample × window (default 10 kb), compute its genetic
   distance to each cluster's consensus (percent of SNPs where the sample ≠ that consensus).
3. **Density-cloud membership.** Plot sample-windows in distance space, draw a 2D kernel
   density cloud per cluster, and flag a window as **introgressed** when a sample's window
   lands inside a *different* cluster's cloud than its own (point-in-polygon on the density
   contours).
4. **Filter false signal:** hypervariable gene families (SICAvar/KIR), windows shared across
   many clusters, and low-count windows.
5. **Headline:** windows introgressed uniquely in the focal group vs the rest of its cluster.

### 2.1 Why pairwise, and what changes

The V1 implementation is intrinsically **three-cluster**: the density plot has two distance
axes (Mf, Mn) with three clusters welded in by name. It cannot scale as ADMIXTURE K climbs.

The revamp keeps the exact density-contour method but applies it **per cluster pair**
(Kx, Ky): a clean 2-cluster 2D plot — distance-to-Kx × distance-to-Ky — flagging a Kx sample
whose window falls in the Ky cloud (and vice versa). This is *more* principled than V1's
two-axes-three-clouds setup, removes all hardcoded cluster names (they come from
`admix_clusters.tsv`), and scales to any K. The old three-way "ambiguous between Mf and Mn"
case is recovered at the aggregation step as "flagged in both the Kx–Mf and Kx–Mn pairs."

### 2.2 Fidelity — follow the original for the two computational cores

The revamp generalises *structure*, not the maths. CC must preserve the legacy computation
**verbatim** for these two pieces; only the cluster naming, pairing, and thresholds change:

- **Distance computation** — 10 kb window binning and the percent-mismatch-to-cluster-
  consensus distance. Follow `scripts/legacy/find_introgression_updated_framework.R` and its
  V1 port `scripts/R/find_introgression.R` (`window_bin()` + the per-(sample, window)
  distance block, ~lines 182–223).
- **2D kernel density contour step** — building the `geom_density_2d` plot,
  `ggplot_build()`-ing it, and extracting the contour polygons + points for the
  point-in-polygon test. Follow `find_introgression.R` ~lines 226–252. Do not re-derive the
  density estimation; reuse the ggplot-build extraction pattern.

### 2.3 Detection call — isolated, pluggable, tunable

The yes/no introgression decision is a **single isolated function**: the sample-window's
position plus the own-cluster and other-cluster contours in → introgression boolean out.
Three selectable rules — two read a fitted density surface, one does not:

- **`absolute` (default — reproduces V1).** Introgressed when the point is inside the *other*
  cluster's contour at density level ≥ `contour_level_other` **and** *not* inside its *own*
  cluster's contour at level ≥ `contour_level_own`. Both default to `5e-4`, which is exactly
  V1 (`find_introgression.R:241, 257–259`). Raising `contour_level_other` (toward the other's
  core) and tuning `contour_level_own` moves toward "deep in the other, shallow in its own."
- **`relative` (lead alternative).** Introgressed when the point is *more core to the other
  cluster than to its own* — i.e. it reaches a higher-density layer in the other than in its
  own. Self-calibrating (no absolute cutoff that shifts with cluster spread); documented as
  the preferred alternative and wired into the test sweep.
- **`distance` (deterministic; added 2026-08-20, Stage 5b).** Introgressed when the window
  matches the *other* cluster's consensus clearly better than its own, judged **directly on
  the two per-window match-fraction distances** — the same `d_own` / `d_other` numbers the
  density rules use as plot coordinates. No kernel surface is fitted and
  `contour_max_level()` is never called, so the answer for a sample-window depends only on
  that sample-window, not on who else is in the cohort. Two margin modes:
  - **fixed** (`distance_adaptive: false`, default) — `d_own - d_other >= distance_margin`,
    in percentage points (config default 15).
  - **adaptive** (`distance_adaptive: true`) — the margin comes from the two clusters' own
    within-cluster spread instead of a config number: introgressed when `d_own` is at or
    above the `distance_adaptive_quantile` (default 0.9) quantile of the own cluster's
    self-distance distribution, `d_other` is at or below the same quantile of the *other*
    cluster's self-distance distribution, and `d_own > d_other`. Data-calibrated but still
    deterministic — quantiles of a fixed input, not a fitted surface.
- The "point in both clusters' peripheries / overlap zone" definition was considered and
  **rejected** — that's ambiguity, not introgression, and it's where hypervariable noise lands.

Rule + thresholds are config (`detection_rule`, `contour_level_other`, `contour_level_own`,
`distance_margin`, `distance_adaptive`, `distance_adaptive_quantile`), so the synthetic test
can sweep them and let the data pick the winner. Keeping the call isolated means swapping or
adding a rule never touches the distance or contour code: `distance` was added without
editing `introgression_core.R` at all, and under it the pipeline skips the contour build
entirely rather than fitting a surface nothing reads.

## 3. Pipeline structure

```
per-pair detection   (Snakemake {pair} wildcard, one pair per invocation, fanned across cores)
        │
        ▼
aggregation          (stitch all pairs → cross-pair summaries)
        │
        ▼
focal-group headline (windows unique to focal_group vs rest of its cluster)
```

- **Per-pair detection** is a single-purpose script run once per (Kx, Ky). Snakemake expands
  the wildcard over all pairs and parallelises — no R-level parallelism.
- **Pairs** default to **all C(K,2)** combinations, derived from the distinct cluster labels.
- **Aggregation** combines per-pair window calls, applies the dataset-level and gene-family
  filters, and writes the summary tables.
- **Headline** generalises V1's hardcoded focal comparison to a configurable `focal_group`.

## 4. Guardrail (all-pairs scales quadratically)

C(K,2) grows fast: K=3→3 pairs, K=6→15, K=8→28, K=10→45, each a full density pass. A
**runtime warning** fires above a configurable pair-count threshold and states **both**:

- **Compute:** "N pairwise comparisons will run" — the cost implication.
- **Interpretation:** many pairwise contrasts carry a multiple-comparison burden, and not
  every cluster pair is biologically meaningful — the user should consider a focal subset.

The same caveat is documented in three places: this spec, the script header, and a dedicated
**Introgression** section in the GitHub README.

## 5. Inputs & config

Inputs (as V1): combined genotype table (Stage 4), `admix_clusters.tsv`, sample metadata
(role columns), reference `.fai` / `contig_map.tsv`, optional GFF for the gene-family filter.

New/changed `introgression` config block:

```yaml
introgression:
  window_size_bp: 10000
  min_snps_per_window: 5
  min_cluster_n: 5
  min_samples_per_window: 2    # strict >: drop windows carried by <= 2 samples (V1's n > 2)
  per_cluster_min_pct: 0.05
  per_cluster_min_samples: 0   # optional ABSOLUTE floor on top of the % (0 = percentage only)
  pairs: "all"                 # "all" (C(K,2)) | explicit list of [Kx, Ky] pairs
  pair_warn_threshold: 15      # runtime warning above this many pairs
  detection_rule: "absolute"   # "absolute" (V1 default) | "relative" | "distance" (§2.3)
  contour_level_other: 5.0e-4  # density level for "in the other cluster" (absolute rule)
  contour_level_own:   5.0e-4  # density level for "in its own cluster" (absolute rule)
  distance_margin: 15          # percentage points (distance rule, fixed mode)
  distance_adaptive: false     # derive the margin from within-cluster spread instead
  distance_adaptive_quantile: 0.9
  focal_group: null            # geography/group value for the unique-windows headline; null = skip headline
  gene_family_filters:         # attribute keywords masked from the GFF (hypervariable families)
    - "SICA"
    - "KIR"
  gff: null                    # optional; gene-family filter skipped if unset
```

Graceful degradation throughout: no GFF → skip the gene-family filter (with a logged note);
`focal_group: null` → skip the headline; `pairs` list lets a user avoid the quadratic blowup.

## 6. Outputs

Per-pair: introgressed-window calls for (Kx, Ky). Aggregated: cross-pair filtered windows,
per-cluster and per-sample summaries, per-province/geography summary, plus the two artifact
masks as window lists (`gene_family_masked_windows.tsv`, `hypervariable_masked_windows.tsv`)
so the focal test can inherit them.

When `focal_group` is set, two more, and the distinction between them is the point:

- **the result** — `focal_<focal>_enriched_windows.tsv`, windows whose introgression is
  enriched in the focal subgroup against a subgroup-scaled permutation null, BH-adjusted at
  `focal_fdr`; plus `focal_<focal>_window_tests.tsv`, the same statistics for *every* tested
  window so a non-significant window's `p_adj` can still be quoted;
- **descriptive only, superseded** — `unique_windows_in_<focal>_with_freq_and_coords.tsv` +
  per-chromosome counts. A set difference with no null and no multiple-testing control; not to
  be cited as a result.

Plus the shoulder + focal-unique figures (role-driven, as in Stage 3b).

## 7. Generalisations from V1 (summary)

| V1 hardcode | Agnostic |
|---|---|
| Mn/Mf/Peninsular cluster names throughout | from `admix_clusters.tsv`; pairwise |
| 2 fixed distance axes / 3-way logic | per-pair 2D density |
| hardcoded focal-group headline | config `focal_group` |
| SICAvar/KIR literal filter | config `gene_family_filters` keyword list |
| `ordered_PKNH_` gsub, PlasmoDB name map | `contig_map.tsv` |
| lab/manual sample exclusions | pre-removed in input-prep (clean input) |

## 8. Testing strategy

No validated reference exists, so ground truth is **manufactured**:

1. **Synthetic positive controls (anchor).** Extend `tests/tiny_cohort/`: take a few
   cluster-A samples and overwrite one window's genotypes with cluster-B's consensus. The
   A–B pair **must** flag exactly that window; untouched samples must **not**. Include a
   clean negative control (no injected event → no/near-zero false positives). These same
   fixtures drive the **detection-rule / threshold sweep**: run `absolute` (across
   `contour_level_*` values) vs `relative` and compare detection-rate against false-positives
   to pick the default — the empirical answer to "which threshold is best."
2. **Component unit tests:** dominant-allele calling, per-window distance, window binning,
   point-in-polygon membership, all three detection rules, the adaptive rule's calibration,
   and the both-directions call assembly — each in isolation
   (`tests/R/test_introgression_units.R`).
3. **External benchmark against a published result (added Stage 5b).** The Malaysian
   *P. knowlesi* dataset in `data/benchmark_malay/` ships its own inputs *and* its published
   introgression result. `scripts/R/introgression_benchmark_malay.R` (target:
   `make benchmark-malay`) runs the whole Stage-5 chain on it once per detection rule and
   scores four things: sample-level concordance (primary), window overlap with the published
   windows (secondary — compare by **CHROM + window start**, never by `WINDOW` id), the
   **cluster-size test** (the discriminator), and perturbation stability under a random 10%
   and 20% sample drop. The fixture is real embargoed data and its outputs are gitignored;
   the harness reuses the pipeline's own core, detector and aggregate scripts rather than
   reimplementing them.
4. **Real-data sanity check:** run on a real cohort and confirm the focal-group story
   survives the reformulation. Documented as a sanity comparison, **not** a bit-match.

## 9. Implementation record

The dated build log — what each increment changed, the benchmark scoreboards, the floor
derivation and the per-cohort validation runs — is kept with the cohort that produced it,
not in this repo. The durable design decisions it reached are stated in §1–8 above.

## 10. Out of scope / open items

- The N-dimensional density generalisation (approach "C") is explicitly rejected — pairwise
  is the chosen path.
- All three detection rules (`absolute`, `relative`, `distance`) are in scope and wired to
  the sweep; the overlap-zone rule is rejected. Adding further rules later is a drop-in
  (§2.3) — `distance` was added in Stage 5b without touching `introgression_core.R`.
- A **Mahalanobis / covariance-aware** fourth rule is deferred: only worth building if the
  scoreboard shows the distance rule's fixed axes are the limiting factor.
- The deeper *detection-limit* sweep (varying injected-segment size/divergence to find the
  method's sensitivity floor) is deferred; the fixed synthetic anchor + rule/threshold sweep
  are sufficient for now.
- The density-contour *statistic* itself is preserved — we keep V1's estimation + extraction
  (§2.2); only the decision rule on top of it is configurable.
- A **rate-stratified focal null** (drawing the permuted subgroup within call-count strata
  rather than uniformly over cluster members) is deferred. It is the fix for the call-rate
  imbalance the focal test documents, and it is only worth building for a cohort whose diagnostic comes
  out on the anti-conservative side — the Indo run's 0.77 is conservative, so stratifying would
  only make the existing calls stronger.
