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
Two selectable rules:

- **`absolute` (default — reproduces V1).** Introgressed when the point is inside the *other*
  cluster's contour at density level ≥ `contour_level_other` **and** *not* inside its *own*
  cluster's contour at level ≥ `contour_level_own`. Both default to `5e-4`, which is exactly
  V1 (`find_introgression.R:241, 257–259`). Raising `contour_level_other` (toward the other's
  core) and tuning `contour_level_own` moves toward "deep in the other, shallow in its own."
- **`relative` (lead alternative).** Introgressed when the point is *more core to the other
  cluster than to its own* — i.e. it reaches a higher-density layer in the other than in its
  own. Self-calibrating (no absolute cutoff that shifts with cluster spread); documented as
  the preferred alternative and wired into the test sweep.
- The "point in both clusters' peripheries / overlap zone" definition was considered and
  **rejected** — that's ambiguity, not introgression, and it's where hypervariable noise lands.

Rule + thresholds are config (`detection_rule`, `contour_level_other`, `contour_level_own`),
so the synthetic test can sweep them and let the data pick the winner. Keeping the call
isolated means swapping/adding a rule never touches the distance or contour code.

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
  per_cluster_min_pct: 0.05
  pairs: "all"                 # "all" (C(K,2)) | explicit list of [Kx, Ky] pairs
  pair_warn_threshold: 15      # runtime warning above this many pairs
  detection_rule: "absolute"   # "absolute" (V1 default) | "relative" (lead alternative, §2.3)
  contour_level_other: 5.0e-4  # density level for "in the other cluster" (absolute rule)
  contour_level_own:   5.0e-4  # density level for "in its own cluster" (absolute rule)
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
   point-in-polygon membership — each in isolation.
3. **Real-data sanity check:** run on the Indo cohort; confirm the Aceh-vs-Peninsular story
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
- **Cluster-size bias, flagged for review.** Windows per sample run *inversely*
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

## 10. Out of scope / open items

- The N-dimensional density generalisation (approach "C") is explicitly rejected — pairwise
  is the chosen path.
- Both detection rules (`absolute`, `relative`) are in scope and wired to the sweep; the
  overlap-zone rule is rejected. Adding further rules later is a drop-in (§2.3).
- The deeper *detection-limit* sweep (varying injected-segment size/divergence to find the
  method's sensitivity floor) is deferred; the fixed synthetic anchor + rule/threshold sweep
  are sufficient for now.
- The density-contour *statistic* itself is preserved — we keep V1's estimation + extraction
  (§2.2); only the decision rule on top of it is configurable.
