# Tiny synthetic cohort — smoke-test fixture

Self-contained test cohort for exercising the whole WGS path (Stages 0–4) in
seconds-to-minutes, without needing the real 18 GB VCF.

## What's in it

- **Reference.** 3 contigs: `chr1` (200 kb), `chr2` (100 kb), `chrMT` (5 kb).
  `chrMT` is excluded via `reference.exclude_contigs: ["MT"]` in the config,
  so downstream sees only 2 nuclear contigs — enough to test contig-map
  wiring end to end.
- **Cohort.** 40 samples: 18 in group A + 18 in group B + 4 controls
  (`ctrl_1..4`). Controls are dropped via `controls.exclude_patterns: ["^ctrl"]`.
- **Ancestry.** 250 "structure SNPs" on chr1 with allele frequencies drawn
  from opposing Beta distributions per group, plus 250 neutral SNPs on chr1
  and 100 more on chr2. ADMIXTURE at K=2 recovers the A/B split reliably.
- **Filter tags.** Every record is `FILTER=PASS`; no failing sites.
- **Metadata.** `sample_id`, `country`, `geography`, `host`, `date` (dates
  are stub ISO strings staggered across 2024–2025).
- **VCF.** ~16 KB bgzipped, indexed with `bcftools index -c`.

Everything is fake — safe to commit.

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

That runs Stages 0–4 against the tiny cohort. Outputs land under
`tests/tiny_cohort/{outputs,logs,reports}/` (all gitignored). Stage 6
(iHS selection) is wired but the config's `selection.models: []` is
empty, so no scans are triggered.

Expected shape:

| Stage | Output                                                     | ~size |
|-------|------------------------------------------------------------|-------|
| 0     | `outputs/setup/{vcf_samples,nuclear_contigs,contig_map}.*` | tiny  |
| 1     | `outputs/qc/snps.qc.vcf.gz`                                | ~16 KB |
| 2     | `outputs/moi/fws_MOI.tsv`                                  | ~1 KB  |
| 3     | `outputs/structure/admix_clusters.tsv`                     | ~1 KB  |
| 3b    | `reports/figures/pca.png`, `njt.png`, …                    | 10s KB |
| 4     | `outputs/ibd/{groupA,groupB}/*.hmm_fract.txt`              | ~few KB |
| 4b    | `reports/figures/ibd_connectivity_*.png`                   | ~100 KB |

The whole thing runs in **a few minutes** on 2 cores. If you break
something in a rule and CI ever gets set up, this is the fixture to run.

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
