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
> rule. The open question this leaves is a filter question, not a rule question: whether
> `per_cluster_min_pct` should be replaced by `per_cluster_min_samples` as the shipped default
> for cohorts with very unequal cluster sizes. That is Jacob's call. Full record: §9.4.
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
per-cluster and per-sample summaries, per-province/geography summary. Headline (if
`focal_group` set): `unique_windows_in_<focal>_with_freq_and_coords.tsv` + per-chromosome
counts. Plus the shoulder + focal-unique figures (role-driven, as in Stage 3b).

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
