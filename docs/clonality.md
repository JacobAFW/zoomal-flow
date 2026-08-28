# Clonality handling — why every frequency result is produced twice

Written 2026-08-28. The user-facing explanation of the `{sampleset}` axis, and the design
record for it.

## The problem

Clonal samples are near-identical genotypes: the same parasite lineage sampled more than once.
A cohort with a clonal block of three carries **one** genotype's worth of independent evidence
at those three rows, not three. Treating them as three is *pseudo-replication*, and it biases
every result that is fundamentally a **frequency or a count**:

| what it distorts | how |
|---|---|
| allele frequencies | a clonal block votes once per replicate, so its alleles look commoner than they are — and MAF filters then keep or drop variants on an inflated count |
| ADMIXTURE ancestry proportions | a component supported mostly by one clonal block is supported by one genotype; K itself can be propped up by replication |
| per-group polyclonality rates | three clonal polyclonal samples from one province count that province's polyclonality three times |
| introgression window support | the per-cluster support floor asks "how many samples carry this window" — replicates answer three times, and the cluster-size denominator is inflated too |
| iHS case/control | a haplotype-frequency statistic, so clonal replicates in either arm distort precisely the quantity it measures |

*Plasmodium* cohorts are frequently clonal — the Indo cohort has 15 clonal groups covering 33
samples — so a pop-gen pipeline that detects clonality (which this one already did, in Stage 4)
and then ignores it downstream is scientifically incomplete.

What clonality does **not** distort is anything computed per sample in isolation. Per-sample
Fws/MOI is unaffected and is deliberately never recomputed: a sample's within-host complexity
does not depend on which other samples are in the file.

## Why both outputs, always

The de-clonalized result is not automatically the "right" one and the full-set result is not
automatically wrong — they answer slightly different questions, and which you want depends on
the claim. So the pipeline **always produces both, side by side**, plus a comparison table per
stage. The comparison is the feature: a reader who does not want to blind-trust the pipeline can
see, per claim, how much of it survives collapsing clonal replicates.

- A metric that barely moves is a result clonality did not drive.
- A metric that moves a lot was substantially pseudo-replication.

Neither is a verdict on its own. The point is that the difference is never invisible.

## The circularity, and the one-pass architecture

Clonality is **detected from the data** (Stage 4 hmmIBD over per-cluster IBD), and detection
needs every sample — you cannot find clonal groups in a set you have already thinned. That is a
genuine circularity, and it is resolved by ordering rather than iteration:

```
QC ─► MOI ─► structure PREP (shared seam: Pk.bed/bim/fam, Pk.smiss, Pk.dups)
                    │
                    ├─► structure/full/ ──► Stage 4 IBD ──► clonal_clusters.tsv
                    │        (backbone)                            │
                    │                                              ▼
                    │                              clonality/exclude_unique.txt
                    │                                              │
                    └─► structure/unique/ ◄────────────────────────┘
                             (reporting)
```

1. The **full set** runs through Stage 4 exactly as it always did.
2. The **full-set clusters remain the operational backbone.** Stage 4's own clonal detection and
   the Stage-5 introgression cluster definitions key off `structure/full/admix_clusters.tsv` in
   *both* arms. De-clonalization changes *who is inside* a cluster, never *what the clusters are*.
3. After Stage 4, the unique-genotypes list is derived.
4. The frequency stages run a **second, one-pass** branch on it.

**One pass, deliberately.** Clonality is not re-derived on the de-clonalized clusters. Doing so
would be circular, and clonal groups are robust to it — this is the field norm. The `unique`
structure output is therefore a *reporting and comparison* product: it must not, and does not,
feed Stage 4 or the Stage-5 cluster definitions.

## The representative rule

One representative per clonal group, chosen **deterministically**:

1. lowest per-sample missingness (`F_MISS` from PLINK's `.smiss`), then
2. alphabetical sample id, as a tie-break.

This is the same rule `find_duplicates.R` uses to pick between replicate pairs in Stage 3 — same
question, same answer, so the pipeline is not making two different "which copy is best"
judgements. Re-running picks the same sample. Without a `.smiss` the rule degrades to
alphabetical only (still deterministic) and says so loudly in the log.

Every sample in **no** clonal group is kept, untouched.

### Straddling clonal groups

A clonal group whose members sit in more than one ADMIXTURE cluster is possible in principle.
The representative is taken from the group's **majority cluster**; a tie between clusters goes to
the alphabetically-first cluster name. The count of straddling groups is reported in the log and
flagged per group in `declonalization_audit.tsv` — never silent. On this pipeline's own output it
should not arise, because clonality is derived per cluster; the script does not assume its input
came from this pipeline.

## A drop-list, not a keep-list

This is the one implementation detail worth knowing, because it is what makes the `full` arm
trustworthy.

The stages have **different natural sample universes**. On the Indo cohort: Stage 2 scores Fws for
978 QC-surviving samples; Stage 3 works on the 807-sample PLINK set; Stage 5 on the 553 clustered
samples. If de-clonalization were expressed as a keep-list, intersecting it against each stage
would silently shrink the `full` arm to whichever universe the list was drawn from — and `full`
would stop reproducing the pre-clonality pipeline.

So the wired artifact is an **exclusion list**:

| file | contents |
|---|---|
| `outputs/clonality/exclude_full.txt` | empty — so `full` is a guaranteed no-op in every stage |
| `outputs/clonality/exclude_unique.txt` | the non-representative members of each clonal group |

"Remove these specific samples" is universe-independent. The readable keep-lists
(`all_genotypes.txt`, `unique_genotypes.txt`) are written too — they are what you read to see
what the sample sets *are* — but the rules consume the drop-list.

## What you get

```
outputs/clonality/
  unique_genotypes.txt              the de-clonalized sample set (readable)
  all_genotypes.txt                 the full sample set (readable)
  exclude_unique.txt                what the rules consume for `unique`
  exclude_full.txt                  empty; what the rules consume for `full`
  declonalization_audit.tsv         one row per clonal group: members,
                                    representative, who was dropped, whether
                                    the group straddled clusters

outputs/structure/{full,unique}/        ADMIXTURE, PCA, NJ distance, clusters
outputs/introgression/{full,unique}/    pairs, filtered windows, focal test
outputs/moi/{full,unique}/              Fws by country / geography
outputs/selection/{full,unique}/<model>/ iHS scan + candidate regions
reports/figures/{full,unique}/          the Stage-3 and Stage-5 figures

outputs/<stage>/declonalization_comparison.tsv    full vs unique, per metric
reports/figures/declonalization_comparison.{png,svg}
```

The per-sample Fws/MOI table (`outputs/moi/fws_MOI.tsv`), the QC outputs, Stage 4, and the shared
structure prep are **not** namespaced — they are sample-set-independent by construction.

## The two arms are not like-for-like at the variant level

Worth knowing before you compare anything, because it is easy to assume otherwise: the
de-clonalized arm is **not** "the same analysis with fewer samples". Stage 3's filter chain is
re-derived per sample set, and two of its steps are *fractions of the cohort* rather than fixed
counts — `--geno` (per-variant missingness) and `--maf`. Shrinking the sample set shrinks their
denominators, so the surviving VARIANT set genuinely differs.

That is deliberate and it is the whole point: allele frequencies are exactly what clonal
pseudo-replication distorts, so the de-clonalized arm has to re-derive them rather than inherit
the full-set variant list. But the effect is larger than the sample change would suggest. On the
Indo cohort:

| step | full | unique |
|---|---|---|
| samples after the duplicate/clonal removal | 784 | 766 (−18) |
| variants after `--geno 0.20` | 75,204 | **60,616 (−19%)** |
| samples after `--mind 0.10` | 553 | 539 |
| variants after `--maf 0.01` | 24,641 | **19,448 (−21%)** |

Removing 2.3% of samples costs 19% of variants, and almost all of that happens at the FIRST
step, `--geno`, not at MAF. The mechanism is a denominator effect: 95% of the 1.44 M raw variants
already fail `--geno 0.20`, so the survivors are packed hard against the threshold. The dropped
clonal replicates have high *overall* missingness (that is why they lost the representative
pick), but at those well-covered surviving variants they are mostly called — so removing them
raises each variant's missing *fraction* and tips a large number over the line.

Two consequences for reading the comparison:

- A cohort-level count falling in the `unique` arm may reflect the smaller variant set as much as
  the smaller sample set. Check `n_samples_analysed` **and** the Stage-3 logs before attributing
  a change to pseudo-replication alone.
- If that trade is not what you want for your cohort, the lever is `structure.max_variant_missing`
  / the `--geno 0.20` first pass, not the clonality config. Nothing here is tuned to hide it.

## How to read the comparison

Every stage's `declonalization_comparison.tsv` has one row per metric with `full`, `unique`,
`delta` (unique − full) and `pct_change`. The first rows are the same everywhere:
`n_samples_analysed`, `n_samples_removed`, `n_clonal_groups_collapsed` — the size of the
intervention. Read that first: if only a handful of samples were removed and a metric still moved
a lot, that metric was resting on very few independent genotypes.

Then per stage:

| stage | metrics |
|---|---|
| structure | `best_k`, `n_clusters`, per-cluster sizes, PC1/PC2 variance |
| introgression | window counts overall and per cluster, focal-enriched windows, the descriptive "unique windows" count |
| moi | per-group polyclonality percentage and group n |
| selection | candidate regions and scanned markers, per model |

Counts usually **fall** — fewer independent genotypes support fewer windows and fewer candidate
regions. A metric that **rises** is worth attention: it normally means a threshold that scales
with cluster size got easier to clear once the clonal padding was removed. Stage 5's per-cluster
support floor is exactly such a threshold (a fraction/floor over cluster *n*), so an
introgression window count going *up* under de-clonalization is expected behaviour, not a bug.

`reports/figures/declonalization_comparison.png` is the same information as a dumbbell per
metric: two coincident dots mean clonality did not drive that number.

## What it found on the Indo cohort (2026-08-28, counts only)

15 clonal groups covering 33 samples, none straddling clusters. 18 replicates dropped
(807 → 789 in the PLINK universe; 553 → 539 clustered).

**Structure.** K stays 3 and the three clusters survive, so the population structure itself is
not an artifact of clonality. The clonality is very unevenly distributed: Mf loses 1 of 410
(−0.2%), Mn 8 of 108 (−7.4%), Peninsular 5 of 35 (−14.3%). PC1 variance rises 48.1% → 52.9%.

**Group summaries.** Aceh is the biggest mover: n 24 → 19 (−21%). Polyclonality rates barely
shift anywhere (≤ 0.4 percentage points), which is expected — dropping a clonal replicate
removes one sample from both the numerator and the denominator of its own group.

**Introgression, cluster level.** Essentially stable: 134 → 133 windows. Per cluster it
redistributes — Mf 106 → 113 (up, because the per-cluster support floor is size-scaled and
Mf barely shrank), Mn 28 → 20 (down 29%, and Mn lost 7.4% of its samples).

**Introgression, focal level — the headline, and it is a negative one.** The Aceh focal
enrichment result of §9.6 (5 windows) goes to **0 windows**. The cause is that **5 of the 10
Aceh samples were clonal replicates of each other** (three Peninsular clonal groups are Aceh /
Sabang samples), so the focal group halves to n = 5.

The two reasons are worth separating, because they are not the same claim:

| window | full (Aceh) | unique (Aceh) | full p_adj | unique p_adj | reading |
|---|---|---|---|---|---|
| w12_2635000 | 7/10 vs 0/25 | 4/5 vs 0/25 | 0.0052 | 0.0551 | proportion held — lost to POWER |
| w14_1225000 | 10/10 vs 6/25 | 5/5 vs 6/25 | 0.0064 | 0.326 | proportion held — lost to POWER |
| w11_385000 | 10/10 vs 7/25 | 5/5 vs 7/25 | 0.0095 | 0.336 | proportion held — lost to POWER |
| w9_2025000 | 7/10 vs 1/25 | 3/5 vs 1/25 | 0.0095 | 0.351 | proportion FELL (70% → 60%) — partly clonal |
| w4_955000 | 5/10 vs 0/25 | 2/5 vs 0/25 | 0.0456 | 0.434 | proportion FELL (50% → 40%) — partly clonal |

So three of the five are *not* shown to be pseudo-replication — the fraction of Aceh carrying
them is preserved or higher, and they fail only because a hypergeometric test on 5 focal samples
cannot clear BH across 302 windows no matter how clean the signal. Two of them additionally lose
proportion, so clonal replicates were part of their support.

The defensible conclusion is the conservative one: **on independent genotypes the Indo cohort has
only 5 Aceh samples, which is too few to support a window-level focal claim at all.** The chr4
≈0.95 Mb window in particular was already marginal (p_adj = 0.046) and is now 0.434 with its
proportion falling — it should not be carried forward.

**Selection.** The published negative result is unchanged: 0 candidate regions in both arms. The
scanned-marker count rises (12 → 286) because rehh discards markers not genotyped in 100% of
haplotypes, and a smaller haplotype set clears that bar more often — a denominator effect, not a
change in signal.

## Turning it off

```yaml
clonality:
  declonalize: false
```

Builds only the `full` arm and skips every comparison. Nothing else changes — `full` is
bit-identical either way, because its exclusion list is empty.

If Stage 4 produced no clonal signal at all, `unique_genotypes.txt` equals `all_genotypes.txt`,
`exclude_unique.txt` is empty, and the `unique` arm simply reproduces `full`. That is a
legitimate result (a cohort with no clonality), not an error, and the comparison table will show
every delta as zero.

## Existing outputs: the layout migration

Cohorts already run before this increment have their results at the old flat paths. Run

```bash
bash scripts/sh/migrate_to_sampleset_layout.sh --dry-run   # inspect
bash scripts/sh/migrate_to_sampleset_layout.sh             # move
```

to move them into `full/` so they are not recomputed. It is idempotent, never overwrites, and
`--undo` reverses it. A fresh cohort needs none of this.

One re-run is unavoidable even after migrating: `final_filters` gained the exclusion-list input,
so Stage 3 re-derives once. That is deliberate — it is the regression check that the `full` arm
still reproduces its previous numbers with an empty drop-list.

**Result of that check on the Indo cohort (2026-08-28).** The re-derived full-set arm is
byte-identical to the pre-clonality outputs: `admix_clusters.tsv`, `best_k.txt` (K = 3),
`pca_variance.tsv`, `Pk.eigenvec` and `clonal_clusters.tsv` all compare equal, and ADMIXTURE
reported the same input dimensions and seed (`Size of G: 553x22785`, `Random seed: 43`). The one
difference is ADMIXTURE's own **cross-validation error at the fifth decimal** on the two K values
that happened to be recomputed (K = 3: 0.29616 → 0.29614; K = 7: 0.32534 → 0.32537). That is
floating-point non-determinism in ADMIXTURE's CV reduction, not an effect of this change — the
`.Q` point estimates those CV values score are themselves byte-identical, and the best-K pick is
unaffected (K = 3 wins by 0.006 over the runner-up).
