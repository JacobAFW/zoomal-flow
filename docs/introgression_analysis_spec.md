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
>    port already fixed with a coordinate-derived `window_id()` (§9.1).
> 2. **The "inverse cluster-size bias" is mostly the per-cluster FILTER, not the rule.**
>    `per_cluster_min_pct` is a *percentage* of cluster n, so its threshold scales with cluster
>    size — on the Indo cohort 5% is 21 samples for Mf (n=410) but 2 for Peninsular (n=35).
>    Swapping it for an equal absolute floor (`per_cluster_min_pct: 0`,
>    `per_cluster_min_samples: 6`) and changing nothing else collapses the effect on both
>    cohorts:
>
>    | cohort | per-cluster filter | median windows/sample (large → small cluster) | max/min |
>    |---|---|---|---|
>    | Malay, `absolute` | 5% of cluster n | Mf 10 · Mn 23 · Pen 36 | 3.6 |
>    | Malay, `absolute` | absolute floor n ≥ 6 | Mf 19 · Mn 18 · Pen 19.5 | **1.08** |
>    | Indo, `absolute`  | 5% of cluster n | Mf 3 · Mn 13 · Pen 24 | 8.0 |
>    | Indo, `absolute`  | absolute floor n ≥ 6 | Mf 6 · Mn 10 · Pen 11 | **1.83** |
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
> claims are cross-cluster / Aceh-vs-Peninsular comparisons where the size confound must go.
> The percentage stays available for within-cluster framing. **The floor `N` is derived, not
> guessed:** set it by a null-permutation / FDR argument — permute the introgression calls
> within cluster, recompute per-window sample-support under the null, take the smallest `N`
> that essentially no null window reaches, and report the null distribution plus a stability
> sweep around it. This also closes the diagnosis's flagged "no FDR anywhere" gap.
>
> **Derivation run 2026-08-24 (§9.5): `N = 32` on the Indo cohort, `N = 44` on the Malay
> benchmark — and the result is that the FDR argument and the focal-group headline cannot both
> hold on this cohort.** The FDR-safe floor scales with cluster size (Indo: Mf n=410 → 32,
> Mn n=108 → 16, Peninsular n=35 → 10), so a single global floor is set by the largest cluster
> and excludes the smallest. At `N = 32` the Peninsular cluster retains **zero** windows and
> the Aceh headline is **empty**; because Aceh is 10 samples inside Peninsular, *any* floor
> above 10 makes an Aceh-unique window structurally unreachable, and the single surviving
> window (chr4 ≈0.95 Mb, 5 Aceh samples) dies at `N = 6`. The shipped configs carry the derived
> values with that consequence spelled out; **the resolution is a science call, not a config
> tweak** — see §9.5. Full Stage-5b record: §9.4.
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
headline product is windows introgressed **uniquely in a focal group** (e.g. Aceh) relative
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
- **Headline** generalises V1's Aceh-vs-Peninsular step to a configurable `focal_group`.

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
  window so a non-significant window's `p_adj` can still be quoted (§9.6);
- **descriptive only, superseded** — `unique_windows_in_<focal>_with_freq_and_coords.tsv` +
  per-chromosome counts. A set difference with no null and no multiple-testing control; not to
  be cited as a result.

Plus the shoulder + focal-unique figures (role-driven, as in Stage 3b).

## 7. Generalisations from V1 (summary)

| V1 hardcode | Agnostic |
|---|---|
| Mn/Mf/Peninsular cluster names throughout | from `admix_clusters.tsv`; pairwise |
| 2 fixed distance axes / 3-way logic | per-pair 2D density |
| `State == "Aceh"` focal headline | config `focal_group` |
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
4. **Real-data sanity check:** run on the Indo cohort; confirm the Aceh-vs-Peninsular story
   survives the reformulation. Documented as a sanity comparison, **not** a bit-match.

## 9. Implementation record (built 2026-08-17)

Where the code lives:

| Piece | File |
|---|---|
| Preserved cores: window binning, consensus, distance, density contours | `scripts/R/introgression_core.R` |
| The yes/no call, both rules, point-in-polygon | `scripts/R/introgression_detect.R` |
| Per-pair driver (`{pair}` wildcard) | `scripts/R/introgression_pair.R` |
| Cross-dataset filters + summaries | `scripts/R/introgression_aggregate.R` |
| Focal-group headline | `scripts/R/introgression_headline.R` |
| Shoulder + focal-unique figures | `scripts/R/plot_introgression_{shoulder,focal}.R` |
| Rule/threshold sweep (not a pipeline rule) | `scripts/R/introgression_rule_sweep.R` |
| Rules | `workflow/rules/05_introgression.smk` |
| Positive control + units | `tests/test_introgression.py`, `tests/R/test_introgression_units.R` |

### 9.1 Departures from V1 — structure, not statistics

1. **One density build per pair, not one per sample.** V1 rebuilt the whole
   ggplot inside a per-sample function. The contour layer's data is the full
   table every time and `coord_cartesian()` is a zoom (it does not touch scale
   limits), so the kde2d grid limits — the union of both layers' data ranges —
   are identical whether the point layer holds one sample or all of them. The
   density estimate is therefore unchanged; the cost drops from O(n_samples)
   ggplot builds to one.
2. **Point-in-polygon per contour ring.** V1 passed every contour vertex above
   the level cut to a single `sp::point.in.polygon()` call, concatenating
   nested rings and separate pieces into one vertex sequence — undefined for
   nested contours. We test each ring separately and take the deepest
   containing level. "Inside at level ≥ L" is then exact, and it is what the
   `relative` rule needs.
3. **`≥` not `>` on the level cut.** V1 wrote `filter(level > 5e-4)`. The
   default contour breaks are `pretty(range(density), 10)`, whose lowest break
   is frequently 5e-4 exactly — so V1's strict `>` silently dropped the
   outermost ring. §2.3 specifies `≥`, which is implemented. At the same
   nominal number the agnostic rule is therefore marginally more permissive
   than V1; raise `contour_level_other` one break to recover V1's region.
4. **Window ids are `w<chrom>_<bin>`, not a row index.** V1's `row_number()`
   over sorted (CHROM, bin) pairs only stays consistent if every run sees the
   same variant set; a coordinate-derived id joins across pairs without a
   shared index table. The `w` prefix is load-bearing: a bare `1:5000` is
   parsed as a clock time by readr's type guesser, which silently rewrites the
   id on round-trip.
5. **GFF↔reference contig mapping now covers the whole karyotype.** V1 mapped
   PlasmoDB `LT727…` accessions onto `ordered_PKNH_NN_v2` by exact sequence
   length. On this cohort's PlasmoDB-68 GFF vs the A1.H.1 **Icor** reference
   only 3 of 14 chromosomes have identical lengths (the Icor corrections shift
   the rest by 99–7138 bp), so V1's gene-family filter masked SICAvar/KIR on
   3 chromosomes and let them through on the other 11. The agnostic version
   tries exact name, then exact length, then **mutual nearest length** within
   1%, logs every inexact pairing, and writes `gff_contig_map.tsv` for audit.
   All 14 map, in karyotype order; masked windows went 87 → 384.
6. **Per-cluster percentage denominator is the cluster's true size** from
   `admix_clusters.tsv`. V1 used "samples with any call at all", which is
   always ≤ cluster size and so slightly more permissive.
7. **Window coordinates are the window's own bounds** (BIN ± window_size/2)
   rather than V1's observed min/max SNP position inside it, so reported
   coordinates don't move with marker density.
8. **`min_samples_per_window` is config**, defaulting to V1's hardcoded `n > 2`.

### 9.2 What the first real run showed (Indo cohort, K=3, absolute @ 5e-4)

Reported as a method result for review, **not** a validated number.

- 3 pairs; 21,944 raw calls over 457 windows → 4,625 calls / 195 windows after
  the four filters. The hypervariable filter is the big one (198 windows).
- Headline: **3 windows unique to Aceh** vs the other 25 Peninsular samples
  (chr4 ≈0.95 Mb in 5 of 10 Aceh samples; chr10 ≈1.43 Mb and chr14 ≈3.14 Mb in
  3 each). V1's own Indo run produced an **empty** headline table, so there is
  nothing to diff against; the qualitative story (chr12/13/14 carrying most of
  the genome-wide signal) matches the legacy script's margin notes.
- **Cluster-size bias, flagged for review — LARGELY EXPLAINED, see §9.4.3.** The Stage-5b
  benchmark showed most of this is the *percentage* per-cluster filter scaling its threshold
  with cluster n, not the density surface; under an equal absolute floor the same Indo run
  gives Mf 6 · Mn 10 · Peninsular 11 (max/min 1.83) instead of the 8.0× gradient below.
  Original observation, kept for the record: windows per sample run *inversely*
  to cluster size: Peninsular (n=35) median 24, Mn (n=108) median 13, Mf
  (n=404) median 3. A large cluster's density cloud is tall and tight — easy
  for another cluster's samples to fall inside — while a small cluster's is low
  and broad, so its own members often sit below their own contour level and
  satisfy the "not in its own cloud" half of the rule. V1's numbers ran the
  other way (Mf median 58), which is consistent with its three-way logic rather
  than with either being right. This is a real property of the method, not a
  bug, and it needs a science call before the numbers are used comparatively.
- `relative` on the same data gives ~20% more calls (26,305 vs 21,944) and
  *widens* the size asymmetry rather than correcting it.
- The shoulder plot shows no shoulder on this cohort: the per-window sample
  count rises smoothly, so `min_samples_per_window: 2` removes only 11 of 457
  windows. The low-n filter is doing almost nothing here.

### 9.3 Sweep result

On the tiny-cohort fixture the sweep is **flat**: every `absolute` threshold
from 1e-4 to 1e-3 and `relative` all recover 6/6 injected sample-windows with
the same 2 (singleton) false positives; 2e-3 loses everything. The fixture's
injected event sits deep enough in the other cluster's core that it does not
discriminate between rules. Discriminating would need the *detection-limit*
sweep (varying injected-segment size and divergence), which §10 defers.

### 9.4 Stage 5b — deterministic rule + method benchmark (built 2026-08-20)

Two things landed: a third detection rule that fits no density surface, and an external
benchmark that scores every rule against a published result so the default can be chosen from
evidence rather than from inheritance.

### 9.4.1 What was built

| File | Role |
|---|---|
| `scripts/R/introgression_detect.R` | `distance` rule + `distance_calibration()` + `pair_calls()` |
| `scripts/R/introgression_benchmark_malay.R` | the scoreboard harness |
| `Makefile` (`make benchmark-malay`) | thin entry point; refuses to run without the fixture |
| `tests/R/test_introgression_units.R` | 25 new unit tests (72 total, all passing) |
| `tests/test_introgression.py` | `distance` recovery gate on the tiny cohort |

Two refactors, both to avoid duplicating logic into the benchmark. `pair_calls()` moved into
`introgression_detect.R` and now owns the both-directions orientation, the adaptive rule's
per-pair calibration and the output schema, so the pipeline and the benchmark call
introgression through the same code path. And `introgression_pair.R` now **skips the contour
build entirely** when the chosen rule needs no density surface, rather than fitting one
nothing reads.

One config addition beyond the rule's own parameters: **`per_cluster_min_samples`** (default
`0`, i.e. no behaviour change). The filter-2 threshold became
`max(ceiling(cluster_n * per_cluster_min_pct), per_cluster_min_samples)`. It exists because the
published Malay chain used an absolute `n > 5` per cluster, which a percentage cannot express
across clusters of different sizes — and, as §9.4.3 shows, the choice between the two turned
out to matter far more than the choice of detection rule.

### 9.4.2 The benchmark harness

`data/benchmark_malay/` (real, embargoed, gitignored) ships the Malaysian inputs *and* the
published result. The harness builds the four fixture tables the pipeline expects (clusters
from `cleaned.fam` + `cleaned.3.Q`, contig map + pseudo-`.fai` from the coordinate bed, a
minimal metadata table), then per scenario and per cluster pair computes the distance table
**once** — it is rule-independent — and calls every rule against it, before shelling out to
the pipeline's own `introgression_aggregate.R` for the filter chain. Nothing is
reimplemented. Whole run: 4 rules × 3 scenarios × 3 pairs in **4.7 minutes**, 5.2 GB peak.
Outputs land in gitignored `outputs/benchmark_malay/`.

Fixture parameters matched: window 10 kb, per-window `n > 5`, per-cluster `n > 5` (via
`per_cluster_min_pct: 0` + `per_cluster_min_samples: 6`), SICA/KIR mask, hypervariable filter.
Cohort: 558 samples — Mf 410, Mn 118, Peninsular 30. Truth: 521 samples / 217 windows.

### 9.4.3 The scoreboard (2026-08-20)

**1. Sample-level concordance (nominally primary) — degenerate on this fixture.** 521 of 558
samples are introgressed in truth, so "call every sample" already scores precision 0.934,
recall 1.0. `absolute`, `relative` and `distance_adaptive` all call **all 558** and therefore
score *exactly* the baseline. Only `distance` declines to call anyone (460 samples, precision
0.930, recall 0.821, Jaccard 0.774) — and its precision is no better than the baseline's, so
it is not enriching for truth samples either. **This measure cannot rank these rules.**

**2. Window overlap vs the published windows (secondary, caveated).** By CHROM + window start:

| rule | windows called | shared with truth (of 217) | Jaccard vs truth |
|---|---|---|---|
| `regen_density` (published fix, reference) | 361 | **138 (64%)** | 0.314 |
| `distance_adaptive` | 496 | 106 (49%) | 0.175 |
| `absolute` | 485 | 81 (37%) | 0.130 |
| `relative` | 356 | 30 (14%) | 0.055 |
| `distance` | 161 | 16 (7%) | 0.044 |

The `regen_density` row is the same density method run three-cluster with its own filter chain,
so it is not perfectly comparable — but the correction it carries is: **64%, not the 19% the
previous status block reported.** That 19% counted `WINDOW`-id matches, and ids are row numbers
that move with the variant set (1 of 138 coordinate-matched windows shares an id). Within the
harness the ranking is `distance_adaptive` > `absolute` > `relative` > `distance`: the rule
least like the density method reproduces the density method's published windows least, which is
what you would expect and is not by itself evidence of being wrong.

**3. Cluster-size test (the intended discriminator) — dominated by the filter.** Median windows
per sample, zero-filled over every cluster member:

| rule | per-cluster filter | Mf (410) | Mn (118) | Pen (30) | max/min |
|---|---|---|---|---|---|
| `absolute` | floor n ≥ 6 | 19 | 18 | 19.5 | **1.08** |
| `distance_adaptive` | floor n ≥ 6 | 11 | 13 | 8 | 1.62 |
| `relative` | floor n ≥ 6 | 8 | 18 | 25.5 | 3.19 |
| `distance` | floor n ≥ 6 | 1 | 6 | 4 | **6.00** |
| `absolute` | 5% of cluster n | 10 | 23 | 36 | 3.60 |
| `distance_adaptive` | 5% of cluster n | 5 | 16.5 | 15.5 | 3.30 |
| `distance` | 5% of cluster n | 1 | 8 | 6 | 8.00 |
| `relative` | 5% of cluster n | 4 | 22 | 41 | 10.2 |

Changing *only* the per-cluster filter moves `absolute` from a 3.6× inverse gradient to
essentially flat (1.08×). The same swap on the **Indo** cohort moves it from 8.0× (Mf 3 · Mn 13
· Peninsular 24) to 1.83× (Mf 6 · Mn 10 · Peninsular 11). So the cluster-size bias recorded in
§9.2 and flagged for a science call is **largely an artifact of `per_cluster_min_pct` scaling
its threshold with cluster n**, not a property of the density surface. It is *not* removed by
dropping the density surface: `distance` shows the largest residual size dependence of any rule.

Report `size_bias_rho` and `size_bias_ratio` together — with K = 3 the Spearman rho has only
four possible values, so direction alone reads a rank swap between two near-equal medians as a
real effect.

**4. Perturbation stability** — window-set Jaccard (CHROM + start) against the full run, after
dropping a random 10% / 20% of samples (seed 20260820):

| rule | drop 10% | drop 20% |
|---|---|---|
| `distance_adaptive` | 0.709 | 0.577 |
| `relative` | 0.668 | 0.583 |
| `absolute` | 0.661 | **0.609** |
| `distance` | 0.647 | **0.500** |

No rule is stable and the deterministic one is the *least* stable. The reason is structural:
the **cluster consensus alleles are recomputed from whoever is in the cohort**, so the
distances themselves move under a cohort shift. `distance` removes the density surface's
contribution to cohort dependence; it does not remove cohort dependence. Its smaller window
set (161 vs 485) also makes its Jaccard more sensitive to a fixed number of borderline flips.

### 9.4.4 Recommendation (Jacob's call — nothing was changed)

**Keep `absolute` @ 5e-4.** The evidence that motivated replacing it does not survive contact
with the benchmark, and the deterministic rule wins none of the four measures. The real finding
is one level down: the **per-cluster filter** is the parameter driving the cluster-size
behaviour, and swapping `per_cluster_min_pct: 0.05` for `per_cluster_min_samples: 6` flattens it
on both cohorts. That is the change worth considering, and it needs a science call — a
percentage is defensible ("a window must be common within its cluster"), an absolute floor is
defensible ("a window needs the same evidentiary support everywhere"), and they answer
different questions. `distance` stays shipped, tested and unselected: it is the honest control
that told us the surface was not the problem.

### 9.4.5 Fixture note

The tiny cohort still cannot discriminate between rules. `distance` recovers 6/6 injected
sample-windows at every margin from 5 to 30 pp and at every adaptive quantile, with the same
2–4 singleton false positives the density rules produce, all of which the low-n filter removes.
The sweep (`detection_rule_sweep.tsv`) now carries `distance` rows for the record.

### 9.5 Stage 5c — deriving the per-cluster floor by permutation (built 2026-08-24)

Stage 5b settled that the per-cluster filter should be an absolute floor rather than a
percentage. This derives the floor instead of guessing it, and closes the "no FDR anywhere"
gap the Malay diagnosis flagged. Two scripts, neither of which re-runs detection — both read
the cached per-pair calls:

- `scripts/R/introgression_floor_derivation.R` → `outputs/introgression/floor_derivation.tsv`
- `scripts/R/introgression_floor_sweep.R` → `outputs/introgression/floor_stability_sweep.tsv`

#### The null

Filter 2 asks whether a window has enough same-cluster support to be a shared event rather
than a coincidence. The null therefore has to be "what does coincidence look like": per
cluster, each sample keeps its *number* of introgressed windows and its *eligible* window set,
but redraws **which** windows it hits, uniformly without replacement. That destroys genomic
co-location between samples while preserving per-sample call rate, per-sample testability and
cluster size — the three things that would otherwise confound a comparison between clusters.
1,000 replicates; each replicate passes through the dataset-level filter (filter 1) exactly as
the real data does before support is measured.

**The eligible set is computed, not assumed, and it is the whole argument.** A sample can only
be called where the pipeline computed a distance for it, i.e. where it has more than
`min_snps_per_window` non-missing calls. The derivation recomputes that eligibility table from
the genotype table (window binning + a per-(sample, window) SNP count — no consensus, no
distances, no density, ~6 s on the Indo cohort). Assume "any 10 kb window in the genome"
instead and the universe roughly triples, the null thins out, and `N` comes out far too small.
On the Indo cohort only **700** of ~2,300 genomic windows carry more than 5 SNPs. One stated
caveat, and it biases `N` downward: the pipeline also drops sample-windows whose two distances
tie, which needs the consensus, so the eligible set here is marginally too large.

#### Result — Indo cohort (counts only)

Universe 700 windows; 19,372 distinct (sample, window) calls; 553 samples.

| cluster | n | null mean support | null p95 | null p99 | null max |
|---|---|---|---|---|---|
| Mf | 410 | 17.95 | 25 | 28 | 42 |
| Mn | 108 | 6.73 | 11 | 13 | 21 |
| Peninsular | 35 | 3.14 | 6 | 7 | 13 |

Target: expected null windows < 1 per cluster **and** FDR < 0.05. Smallest floor meeting it:

| cluster | n | FDR-safe N | expected null windows | observed windows | FDR |
|---|---|---|---|---|---|
| Mf | 410 | **32** | 0.994 | 120 | 0.0083 |
| Mn | 108 | 16 | 0.724 | 109 | 0.0066 |
| Peninsular | 35 | 10 | 0.368 | 86 | 0.0043 |

**Chosen global `N = 32`**, binding cluster Mf. The Malay benchmark reproduces the pattern
independently: universe 1,363 windows, `N = 44` (Mf 44 / Mn 20 / Peninsular 10).

#### Result — stability sweep (Indo, counts only)

Windows surviving the full filter chain, by floor:

| N | Mf | Mn | Peninsular | all | Aceh headline |
|---|---|---|---|---|---|
| 1 | 9 | 6 | 22 | 37 | 1 |
| 3 | 35 | 26 | 40 | 101 | 1 |
| 5 | 55 | 51 | 47 | 153 | 1 |
| 6 | 70 | 63 | 43 | 176 | **0** |
| 10 | 99 | 70 | 23 | 192 | 0 |
| 16 | 127 | 59 | 7 | 193 | 0 |
| 20 | 118 | 49 | 1 | 168 | 0 |
| 25 | 115 | 39 | 1 | 155 | 0 |
| 30 | 109 | 32 | 0 | 141 | 0 |
| **32** | **106** | **28** | **0** | **134** | **0** |
| 35 | 107 | 23 | 0 | 130 | 0 |

Two things to read here. The total is **not monotone** in `N` — it rises to 210 around N = 12
before falling — because the hypervariable filter (filter 4) runs *after* the floor: raising the
floor strips a cluster's weak support for a window, and a window that was called in two clusters
becomes called in one and survives filter 4. "Raise the floor to be safer" is therefore not a
one-way lever. And the Mf/Mn columns are locally flat around N = 32 (106 → 105 → 106 → 107), so
the derived floor is not perched on a cliff *for the large clusters*. Peninsular has been at
zero since N = 30.

#### The finding: the FDR floor and the headline are incompatible here

Aceh is **10 samples inside the 35-sample Peninsular cluster**. An Aceh-unique window can
therefore never have support above 10, so **no floor above 10 can produce an Aceh headline at
all**, and the single window that does survive under a low floor — chr4 ≈0.95 Mb, carried by 5
of the 10 Aceh samples — dies at `N = 6`. The FDR-safe floor is 32.

This is not a tuning problem. It says the detection step calls too liberally for a
support-based filter to rescue at these cluster sizes: 4–9% of every cluster's sample-windows
are called introgressed, so with 410 Mf samples drawing from 700 windows, chance alone puts
~18 samples on every window. A floor strong enough to beat that is stronger than any
10-sample subgroup can supply.

Options, in the order they should be considered. **Options 1 and 2 were both taken** in
Stage 5d (§9.6): the cluster floor stays exactly as derived here, and the focal group gets its
own test at its own scale rather than being filtered through that floor. Options 3 and 4 remain
untaken.

1. **Accept the floor and retire the focal-group headline** as a window-level claim. The
   cross-cluster result (134 windows, Mf/Mn) stands; "windows unique to Aceh" does not.
2. **Test the focal group directly** rather than filtering it through its cluster's floor —
   an Aceh-vs-rest-of-Peninsular contrast with its own permutation null over 10 samples,
   which is a different and much better-posed test than "survive a cluster-wide floor".
   **BUILT — see §9.6.**
3. **Reduce the call rate first** (a stricter detection threshold), then re-derive. The null
   support scales with the call rate, so halving it roughly halves the required floor.
4. **Per-cluster floors** (Mf 32 / Mn 16 / Peninsular 10), each FDR-safe in its own cluster.
   Explicitly rejected in the Stage 5c brief as reintroducing size-dependence — recorded here
   because it is the only option that keeps a Peninsular result while controlling FDR.

The shipped default is the derived global floor, with the consequence documented in
`config/cohort.example.yaml`. `config/config.yaml` does **not** inherit the number: an
FDR-safe floor is cohort-specific (Indo 32, Malay 44) and a 410-sample cluster's threshold
would return nothing on a small cohort, so the generic config ships the mechanism plus the
minimum that keeps the per-cluster test no weaker than the pooled one
(`min_samples_per_window + 1` = 3) and points at the derivation script.

### 9.6 Stage 5d — the focal-group enrichment test (built 2026-08-28)

§9.5 ended on a finding, not a setting: a focal-group claim cannot survive a cluster-wide FDR
floor, because a 10-sample subgroup inside a 35-sample cluster can never reach a floor
calibrated on a 410-sample cluster. The conclusion Stage 5d draws from that is architectural.
The cluster-level question and the focal-level question are **two different tests at two
different scales, and each needs its own null**:

| | question | null | control | output |
|---|---|---|---|---|
| **cluster level** | which windows show gene flow between clusters | each sample redraws its calls over its own eligible windows (§9.5) | derived per-cluster support floor, FDR target | `introgressed_windows_filtered.tsv` |
| **focal level** | is a window's introgression *enriched* in this subgroup vs the rest of its cluster | size-preserving focal/background **label** permutation inside the cluster | BH across the cluster's windows, `focal_fdr` | `focal_<group>_enriched_windows.tsv` |

**Nothing about the cluster floor changed.** `per_cluster_min_samples` is still the
permutation/FDR-derived number of §9.5, the cross-cluster result is still the 134 windows it
produced, and the focal test does not feed back into it. The two are independent on purpose —
this is not a route to relax the floor.

This is a **general** design, not an Aceh rescue. Any agnostic cohort with a focal subgroup —
a district, a host species, a sampling year — has the same shape of problem: the subgroup is
smaller than its parent cluster, so a cluster-scaled threshold cannot speak about it. The
subgroup gets a null sized to the subgroup. `focal_group` names the value, `focal_role` names
the metadata role it is a value of (default `geography`), and the parent cluster is derived
from the cluster assignment. No name literals anywhere.

#### The test

`scripts/R/introgression_focal_test.R`, statistic in `scripts/R/introgression_focal_core.R`.
Reads the **cached per-pair calls** — detection is not re-run.

1. **Scope.** Resolve the focal group's parent cluster from `admix_clusters.tsv`; if the group
   straddles more than one cluster the run errors, because "the rest of its cluster" is then
   undefined and a null built on the wrong comparison set gives a wrong p-value rather than a
   rough one. (The descriptive headline step resolves the same ambiguity by majority vote; the
   test does not.) Everything below is inside that one cluster.
2. **Filters: the artifact masks, and no support floor.** The gene-family mask (filter 3) and
   the hypervariable/multi-cluster mask (filter 4) are inherited from the aggregate step, which
   now writes both as window lists so the dependency is explicit rather than recomputed. Those
   two excise things that are not introgression at all, so applying them is not a support
   judgement. Neither `min_samples_per_window` nor `per_cluster_min_samples` is applied: both
   are support floors, and this test carries its own multiple-testing control. "Hypervariable"
   is taken at the cluster scale on purpose — that is the scale at which "this window is
   variable in more than one cluster" is a meaningful claim.
3. **Statistic.** Per window, the number of focal samples called introgressed.
4. **Null.** Hold each window's set of introgressed samples in the cluster fixed, and randomly
   reassign which `n_focal` of the cluster's members carry the focal label. 1,000 replicates at
   a stated seed. The question is *given the cluster's own introgression pattern, would a random
   subgroup of this size show this much support at this window by chance?* — so a window
   introgressed cluster-wide shows **no** focal enrichment, and one carried mostly by the focal
   group does. Enrichment **subsumes** V1's "unique" as the degenerate case where background
   support is 0, and unlike "unique" it survives a single background carrier.
5. **Call.** BH across the cluster's windows; enriched at `focal_fdr` (default 0.05).

**Two p-values, one null — this matters for reading the table.** Permuting the label with the
called set held fixed makes each window's marginal null exactly hypergeometric,
`focal_support ~ Hyper(cluster_n, n_called, n_focal)`. The permutation and `phyper()` therefore
describe the *same* distribution — the second evaluates what the first samples. Both are
reported. `p_perm` is the Monte-Carlo estimate, `(b+1)/(B+1)`; its resolution floor is
`1/(B+1)` ~ 1e-3 at B = 1,000, and after BH across a few hundred windows that floor alone can
put significance out of reach no matter how strong the signal. So the **exact** tail is what BH
adjusts and what the call is made on (`p_raw` -> `p_adj`), and `p_perm` is kept as a live check
that the analytic form is the null the method claims to run. Each run prints
`max |p_perm - p_raw|` over the windows where the Monte Carlo has resolution; on the Indo run
that is 0.034, about 2 Monte-Carlo standard errors, i.e. the two agree.

#### Stated limitation — call-rate imbalance, and which way it bites

The label permutation treats every cluster member as equally likely to be called anywhere, so
it does **not** preserve per-sample call rate: the null assumes the focal group's call rate is
the cluster average. When it is not, the direction decides whether the result is threatened or
reinforced, and each run prints the diagnostic and names the direction:

- focal called **more** often than average -> the null under-states the support a random
  subgroup would show, the test **over-calls**. The dangerous direction; calls are provisional
  and a rate-stratified null is the fix.
- focal called **less** often than average -> the null over-states it, the test **under-calls**.
  Surviving windows cleared a higher bar than their p-value implies, and an empty result may be
  a power limit rather than an absence.

#### Result — Indo cohort, `focal_group: Aceh` (counts only)

Aceh is 10 clustered samples inside the 35-sample Peninsular cluster. Universe after the two
artifact masks (384 gene-family windows, 12 hypervariable) and no support floor: **294 windows**
carrying at least one Peninsular call. Call-rate diagnostic: Aceh mean 47.7 windows/sample vs
background 62.0, ratio **0.77** — the **conservative** direction, so these five cleared a bar
slightly higher than their p-values imply.

| window | contig:coords | Aceh | rest of Peninsular | p_raw | p_adj |
|---|---|---|---|---|---|
| w12_2635000 | PKNH_12:2,630,000-2,639,999 | 7/10 | 0/25 | 1.8e-05 | 0.0052 |
| w14_1225000 | PKNH_14:1,220,000-1,229,999 | 10/10 | 6/25 | 4.4e-05 | 0.0064 |
| w11_385000 | PKNH_11:380,000-389,999 | 10/10 | 7/25 | 1.1e-04 | 0.0095 |
| w9_2025000 | PKNH_09:2,020,000-2,029,999 | 7/10 | 1/25 | 1.3e-04 | 0.0095 |
| **w4_955000** | **PKNH_04:950,000-959,999** | **5/10** | **0/25** | **7.8e-04** | **0.046** |

**The chr4 ~0.95 Mb window — the one §9.5 recorded as dying at floor N = 6 — is
focal-enriched, at p_adj = 0.046.** It is the weakest of the five and sits just inside the 0.05
line; it should be read as a marginal call, not a robust one. Under the shipped cluster floor
(N = 32) the descriptive `unique_windows_in_Aceh_*.tsv` table is **empty**, so all five of
these are windows the old rule could not have reported. Two of them (w14, w11) have substantial
background support and were never "unique" under any floor — they are enrichment findings a set
difference cannot express.

Read with §9.5's standing caveat: the *detection* step calls liberally (4-9% of each cluster's
sample-windows), and this stage has never been validated against a reference output. The focal
test controls the multiple-testing burden of asking 294 window questions; it does not fix a
liberal detector, and it inherits whatever the detector called.

#### Synthetic ground truth

`tests/tiny_cohort/` injects a known event in six group-A samples at chr2:45,000, and those six
are also the fixture's `RegionA2` geography — so they double as a focal subgroup with a known
right answer. The focal test flags **that window and nothing else** (6/6 focal, 0/19 background,
p_adj = 1.1e-05), checked independently of the cluster floor. The statistic itself is unit-tested
against hand-computable cases in `tests/R/test_introgression_units.R`, including the two
properties the design rests on: a cluster-wide window scores p = 1 (no enrichment), and a single
background carrier does not destroy a signal the "unique" rule would have discarded.

#### The old headline is superseded

`unique_windows_in_<focal>_with_freq_and_coords.tsv` and its per-chromosome companion are still
written, and the figure still renders, but they are **descriptive only**. A set difference has
no null, no p-value and no multiple-testing control, and it is brittle in exactly the wrong way
— one background sample carrying a window removes it. The script header, the Snakemake rule
docstring, the run's own stdout and the config comments all now say so explicitly, and the
significance claim lives only in the enrichment table. This closes the diagnosis's FDR gap at
the focal level: before Stage 5d there was no multiple-testing control anywhere in the focal
step.

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
  imbalance §9.6 documents, and it is only worth building for a cohort whose diagnostic comes
  out on the anti-conservative side — the Indo run's 0.77 is conservative, so stratifying would
  only make the existing calls stronger.
