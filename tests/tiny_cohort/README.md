# Tiny synthetic cohort — smoke-test fixture + Stage 5 positive control

Self-contained test cohort for exercising the whole WGS path (Stages 0–5) in
seconds-to-minutes, without needing the real 18 GB VCF. It is also the **ground
truth** for introgression detection, which has no validated reference output.

## What's in it

- **Reference.** 3 contigs: `chr1` (200 kb), `chr2` (100 kb), `chrMT` (5 kb).
  `chrMT` is excluded via `reference.exclude_contigs: ["MT"]` in the config,
  so downstream sees only 2 nuclear contigs — enough to test contig-map
  wiring end to end.
- **Cohort.** 54 samples: 25 in group A + 25 in group B + 4 controls
  (`ctrl_1..4`). Controls are dropped via `controls.exclude_patterns: ["^ctrl"]`.
- **Ancestry.** 250 "structure SNPs" on chr1 + 50 on chr2, strongly divergent
  between groups (AF 0.9 in A vs 0.1 in B), plus 250 neutral SNPs on chr1 and
  50 on chr2. ADMIXTURE at K=2 recovers the A/B split reliably.
- **Geography.** Group A splits into `RegionA1` (19) and `RegionA2` (6); group
  B is all `RegionB1`. Cluster labelling is `auto` (majority geography), so the
  clusters come out as `RegionA1` / `RegionB1`.
- **Injected introgression.** The six `RegionA2` samples carry a group-B
  donor's genotypes across one 10 kb window on chr2 (40 000–49 999, 30 SNPs) —
  a transferred haplotype. `data/introgression_truth.tsv` records exactly which
  samples, which donor, and which window. `RegionA2` is also the config's
  `introgression.focal_group`, so the same event drives the headline test.
- **Filter tags.** Every record is `FILTER=PASS`; no failing sites.
- **Metadata.** `sample_id`, `country`, `geography`, `host`, `date` (dates
  are stub ISO strings staggered across 2024–2025).
- **VCF.** ~630 records, ~20 KB bgzipped, indexed with `bcftools index -c`.

Everything is fake — safe to commit.

### Why the fixture looks the way it does

Two design points are forced by Stage 5's detection method, and both were
learned by watching the positive control fail:

1. **Injection copies a donor sample, not the group-B consensus.** A
   consensus-overwritten sample sits 0% away from B's consensus while real B
   samples sit 10–20% away, so a consensus-injected point lands off the *edge*
   of B's density cloud instead of inside it — undetectable by any rule.
   Copying a donor puts the injected sample exactly where that donor already
   is: in the cloud's core. It is also what introgression physically is.
2. **Between-group divergence is uniform across windows.** The density clouds
   are estimated over all windows at once, so a window whose divergence is far
   from the genome-wide norm forms its own sparse lobe that no fixed density
   threshold reaches. Earlier revisions drew structure-SNP frequencies from
   wide Beta distributions; fixed 0.9/0.1 frequencies keep every window
   comparable and the two clouds tight and well separated.

## Regenerate

```bash
# Activate the env first so bcftools is on PATH.
source ../../../envs/activate.sh
python tests/tiny_cohort/generate.py
```

Deterministic seed (`RNG_SEED = 20260803`) — the output is byte-identical
across regenerations on the same platform.

## Run the smoke test

```bash
cd agnostic/
snakemake --configfile tests/tiny_cohort/config.yaml --cores 2
```

That runs Stages 0–5 against the tiny cohort (58 jobs, ~25 s on 4 cores).
Outputs land under `tests/tiny_cohort/{outputs,logs,reports}/` (all
gitignored). Stage 6 (iHS selection) is wired but the config's
`selection.models: []` is empty, so no scans are triggered.

Expected shape:

| Stage | Output                                                     | ~size |
|-------|------------------------------------------------------------|-------|
| 0     | `outputs/setup/{vcf_samples,nuclear_contigs,contig_map}.*` | tiny  |
| 1     | `outputs/qc/snps.qc.vcf.gz`                                | ~16 KB |
| 2     | `outputs/moi/fws_MOI.tsv`                                  | ~1 KB  |
| 3     | `outputs/structure/admix_clusters.tsv`                     | ~1 KB  |
| 3b    | `reports/figures/pca.png`, `njt.png`, …                    | 10s KB |
| 4     | `outputs/ibd/Region{A1,B1}/*.hmm_fract.txt`                | ~few KB |
| 4b    | `reports/figures/ibd_connectivity_*.png`                   | ~100 KB |
| 5     | `outputs/introgression/introgressed_windows_filtered.tsv`  | ~1 KB  |

If you break something in a rule and CI ever gets set up, this is the fixture
to run.

## Verify the Stage 5 positive control

```bash
pytest tests/test_introgression.py -v
```

The pipeline must recover exactly the injected event: all six `RegionA2`
samples called at `chr2:40000-49999`, nobody else called there, and — after
the cross-dataset filters — nothing else called at all. The current fixture
produces 8 raw calls (6 true + 2 singleton false positives, both removed by
the low-n filter) and a headline of exactly that one window.

The same fixture drives the detection-rule sweep; see the Introgression
section of `agnostic/README.md`.

## Design intent

- Every rule in the pipeline is exercised exactly once at a scale that
  fails fast on structural bugs (missing rule wiring, wildcard errors,
  input-function bugs) without waiting for the full 18 GB dataset.
- The cohort is small enough to keep the fixture under 1 MB total but
  large enough to actually recover structure at K=2, run a proper CV
  scan up to K=3, and produce non-trivial IBD networks.
- No cohort-specific hardcodes in the config — the `roles` block uses
  `sample_id`/`country`/`geography` verbatim as column names, which
  exercises the validator's role-rename path.
