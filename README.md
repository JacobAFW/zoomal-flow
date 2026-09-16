# ZOOMAL-Flow

A data-agnostic population-genetics pipeline for zoonotic malaria. The same
six analysis stages (QC → MOI → structure → IBD → introgression → selection)
run for any cohort + reference, with every cohort-specific assumption lifted
to config or derived from inputs. Originally refactored from the V1 Indonesia
*P. knowlesi* pipeline.

This increment ships **Stage 0 (setup)** and **Stage 1 (QC)** only. The
WGS path is implemented; the microhap path is a stubbed seam.

See `DESIGN.md` for the full architecture.

---

## What's included — and what's deliberately not

This repository contains **code only** — workflow rules, scripts, schemas, and
docs. By design it does **not** track:

- raw or processed **data** (VCFs, reference FASTA/index, PLINK sets, any
  sequence/genotype/tabular records)
- the **sample metadata table** (`data/metadata/samples.tsv`) — it carries
  individual-level fields under publication embargo
- **run artefacts** — `outputs/`, `logs/`, `reports/`, `.snakemake/`
- credentials, tokens, or environment files
- the **solved environment** itself (`.pixi/`, ~2.2 GB) — only the `pixi.toml`
  manifest and the `pixi.lock` that pins it are tracked

Data lives outside version control. The pipeline expects it under `agnostic/data/`
at the paths the active config points to (see Quick start). `.gitignore` is
default-deny (all data/genomics/secret/sample-sheet patterns) so a stray input
can't be committed by accident.

**One deliberate exception:** the public example under `examples/` — a VCF, a
sample table, and a mask, all of which are already-public data (see "Pull and
play" below). Those files are un-ignored individually, by name, so default-deny
still applies to everything around them.

---

## Quick start

```bash
cd agnostic

# 1. Build the environment. ZOOMAL-Flow owns its own pinned env — pixi.toml +
#    pixi.lock in this directory. `install` gets the locked conda stack;
#    `setup` adds the six tools no conda channel publishes (moimix, rehh,
#    rnaturalearthhires, and on macOS SeqArray/ADMIXTURE/PLINK2).
pixi install
pixi run setup
pixi run check-env      # prints a version for every tool; fails if any is missing

# 2. Put your cohort inputs under agnostic/data/ (this dir is gitignored):
#      data/vcf/<cohort>.vcf.gz        bgzipped + indexed
#      data/reference/<ref>.fasta      with a .fai beside it
#      data/metadata/samples.tsv       one tidy row per sample (see "Role-based
#                                      metadata" below)
#    For the worked Indo example these already resolve: data/vcf and
#    data/reference are symlinks into the V1 data tree, and
#    data/metadata/samples.tsv is a frozen copy of the V1 sample table.
#    (Large inputs are symlinked, not copied; sample metadata stays gitignored.)

# 3. Configure. Either edit config/config.yaml in place, or start from the
#    worked Indo example:
cp config/cohort.example.yaml config/config.yaml

# 4. Dry-run to validate config + DAG, then run.
pixi run snakemake --cores 1 -n     # validates config; prints the DAG
pixi run snakemake --cores 8        # actually runs
```

To run with a different config file:

```bash
pixi run snakemake --cores 8 --configfile path/to/my_cohort.yaml
```

`pixi run <cmd>` runs a single command inside the environment. To work in it
for a while instead, `pixi shell` once and then call `snakemake` directly.

No pixi at your site? `env/environment.yml` is a conda/mamba fallback — you get
a fresh solve rather than the recorded one. Both paths are documented in
**`env/README.md`**, along with the version table, the one deviation from V1
(hmmIBD), and the verification run behind it.

The V1 Indonesia pipeline keeps its own separate `vvg-box` env in
`../envs/` — nothing here touches it, and nothing here needs it.

---

## Pull and play — try it on public data first

The repository ships a complete, runnable cohort, so you can see the whole
pipeline work before you point it at your own data:

```bash
cd agnostic
pixi install && pixi run setup

# fetch the P. knowlesi A1-H.1 reference from PlasmoDB into
# examples/data/reference/  (navigation path + why there is no direct link:
# examples/README.md), then check you got the right assembly:
pixi run bash examples/check_reference.sh

pixi run snakemake --configfile examples/config.yaml --cores 4
```

**71 publicly-archived ENA isolates of *P. knowlesi* from Malaysia, 18,279
SNPs, ~4 MB.** ADMIXTURE recovers K = 3 by cross-validation, and all 71 samples
land in the same cluster they occupy in the full 981-sample cohort the example
was drawn from. Every stage runs — QC, Fws/MOI, structure, IBD, clonality,
introgression — and both the full and de-clonalized arms are built.

This is the one place the repository deliberately commits a VCF and a sample
table. It is safe because every sample is a published ENA accession and every
metadata field is public; `examples/build_example.sh` re-derives the files and
refuses to write them unless the sample list shares zero IDs with the embargoed
cohort. The `.gitignore` un-ignores those files **by name**, so default-deny
still covers anything else that lands in `examples/data/`.

The reference genome is not shipped — it is ~24 MB that PlasmoDB already
versions properly. `examples/README.md` has the accession list, the provenance,
what is deliberately switched off and why, and how the example's contig names
were reconciled with PlasmoDB's.

Treat it as a demonstration that the pipeline runs and your environment is
correct — not as a result about *P. knowlesi*.

---

## What's in this increment

| Stage | File | Status |
|---|---|---|
| 0  setup            | `workflow/rules/00_setup.smk`       | implemented |
| 1  QC (common)      | `workflow/rules/01_qc.smk`          | implemented |
| 1  QC (WGS)         | `workflow/rules/qc_wgs.smk`         | implemented |
| 1  QC (microhap)    | `workflow/rules/qc_microhap.smk`    | stub (fails fast) |
| 2  MOI / Fws        | `workflow/rules/02_moi.smk`         | implemented (WGS path) |
| 3  structure        | `workflow/rules/03_structure.smk`   | implemented |
| 3b structure figures| `workflow/rules/03_figures.smk`     | implemented |
| 4  IBD              | `workflow/rules/04_ibd.smk`         | implemented |
| 4b IBD figures      | `workflow/rules/04_figures.smk`     | implemented |
| 5  introgression    | `workflow/rules/05_introgression.smk` | implemented (see below) |
| 6  selection (iHS)  | `workflow/rules/06_selection.smk`   | implemented |
| 7  report           | `workflow/rules/99_report.smk`      | implemented |

Stage 0 produces: `outputs/setup/vcf_samples.txt`, `nuclear_contigs.txt`,
`contig_map.tsv`, and a canonical role-renamed `outputs/metadata/samples.tsv`.
Stage 1 (WGS) produces: `outputs/qc/snps.qc.vcf.gz` (the shared seam output),
`outputs/qc/variant_count.txt`, and four QC density panels under
`reports/figures/`.

`docs/WALKTHROUGH.md` (generated from the rule docstrings) is the step-by-step
tour of every rule, its tunables, and something to try at each step.

---

## Clonality — every frequency result is produced twice

Clonal samples are near-identical genotypes. Counting each as independent evidence is
pseudo-replication, and it biases everything that is fundamentally a frequency or a count:
allele frequencies, ADMIXTURE proportions, per-group polyclonality rates, introgression window
support, iHS case/control. *Plasmodium* cohorts are frequently clonal, so the pipeline treats
this as first-class rather than optional.

Four stages — structure, introgression, group summaries, selection — are each built **twice**,
under a `{sampleset}` axis:

```
outputs/<stage>/full/     every sample (the operational backbone)
outputs/<stage>/unique/   one representative per clonal group (the comparison arm)
outputs/<stage>/declonalization_comparison.tsv    the two, side by side
```

Both arms always. The comparison is the feature: you can see, per claim, how much of it survives
collapsing clonal replicates, rather than taking the pipeline on faith. Clonality is detected on
the **full** set (it has to be — you cannot find clonal groups in a set you already thinned), and
the full-set clusters stay authoritative for Stage 4 and the Stage-5 cluster definitions. One
pass, never iterated.

Rationale, the architecture, the representative rule, and how to read the comparison:
**`docs/clonality.md`**. Turn it off with `clonality.declonalize: false`.

Upgrading a cohort that was run before this existed:
`bash scripts/sh/migrate_to_sampleset_layout.sh` (idempotent, `--undo` reverses it).

---

## Introgression (Stage 5) — read before you run it

Stage 5 asks, per genomic window: does a sample carry the genetic signature of
a cluster *other* than the one it is assigned to? It works **per cluster pair**
— for (Kx, Ky) it measures every sample-window's distance to both clusters'
consensus alleles, draws a 2D kernel-density cloud per cluster in that distance
space, and flags a window when a sample lands in the other cluster's cloud
rather than its own.

Stage 5 answers **two** questions, at two different scales, each with its own
null. Do not read one as the other:

- **Cluster level** — which windows show gene flow between clusters. Controlled
  by the per-cluster support floor, derived by permutation with an explicit FDR
  target (`scripts/R/introgression_floor_derivation.R`, spec §9.5). The floor is
  cohort-specific: derive it, don't inherit a number.
- **Focal level** — is a window's introgression *enriched* in a subgroup
  (`focal_group`, e.g. one district) relative to the rest of its own cluster.
  This has its own size-preserving label-permutation null, scaled to the
  subgroup, with its own BH/FDR (`scripts/R/introgression_focal_test.R`, spec
  §9.6). Result: `focal_<group>_enriched_windows.tsv`.

The two are independent by design. A subgroup can never clear a floor
calibrated on a cluster many times its size, so filtering the focal question
through the cluster floor is a mis-posed test, not a threshold to tune.

⚠ The older `unique_windows_in_<focal>_*.tsv` table is **descriptive only** and
is superseded. "Unique to the focal group" is a raw set difference: no null, no
p-value, no multiple-testing control, and one background sample carrying the
window flips it out of the set. Cite the enrichment table instead.

Full method, design decisions and validation status:
**`docs/introgression_analysis_spec.md`**.

Two things to know before running it:

- **This stage has never been validated against a reference output.** V1's
  introgression script was written after the HPC run and its own header records
  that reference outputs were unavailable. Correctness here rests on the
  synthetic positive control in `tests/tiny_cohort/` (a known injected
  introgression event that the pipeline must recover exactly), not on
  reproducing V1. Treat real-cohort numbers as a method result to review, not a
  validated figure.
- **`pairs: "all"` is quadratic, in two different ways.** All-pairs means
  C(K,2) comparisons: K=3 → 3, K=6 → 15, K=8 → 28, K=10 → 45. Each is a full
  density pass, so *compute* grows fast. Just as importantly, each is another
  contrast: run 45 of them and windows will clear any fixed threshold by
  chance, and not every cluster pair is a biologically meaningful comparison in
  the first place. The workflow prints a warning above
  `introgression.pair_warn_threshold`; the fix is to set an explicit
  `introgression.pairs: [[Kx, Ky], …]` list for the contrasts you actually
  intend to interpret.

The detection call is pluggable (`introgression.detection_rule`): `absolute`
(V1's fixed density cutoffs, the shipped default) or `relative` (deeper in the
other cluster's cloud than in its own — no absolute cutoff to drift when one
cluster is more diffuse than another). To see how the two behave on data with
known ground truth:

```bash
Rscript scripts/R/introgression_rule_sweep.R \
  --genotype-table tests/tiny_cohort/outputs/ibd/combined/hmmIBD_input.tsv \
  --clusters       tests/tiny_cohort/outputs/structure/full/admix_clusters.tsv \
  --contig-map     tests/tiny_cohort/outputs/setup/contig_map.tsv \
  --truth          tests/tiny_cohort/data/introgression_truth.tsv \
  --pair           RegionA1__RegionB1 \
  --window-size 10000 --min-snps 5 \
  --out tests/tiny_cohort/outputs/introgression/detection_rule_sweep.tsv
```

It reports detection rate against false positives across a threshold grid. It
changes no defaults — picking a rule is a judgement call, and the sweep is the
evidence for making it.

---

## Small cohorts — what works, and what does not yet

The pipeline is built to run on any cohort, but two of its tools have a hard
**50-sample floor** because below that they cannot estimate what they need to.
One of those is handled; the other is not yet.

**Handled — LD-pruning (Stage 3).** PLINK2 refuses `--indep-pairwise` below 50
samples. Below `structure.min_samples_for_ld_prune` (default 50, PLINK2's own
floor) the pipeline now **skips the prune**, passes the unpruned variant set to
ADMIXTURE, and says so loudly: a banner in
`logs/structure/{sampleset}/ld_prune.log`, a machine-readable
`outputs/structure/{sampleset}/ld_prune_status.txt`, and a warning callout in
the report. Read those before quoting the structure results — linked variants
are correlated evidence, so they can inflate apparent structure, and K
selection in particular should be treated as indicative.

**Not handled — PCA (Stage 3).** PLINK2 `--pca` needs allele frequencies and
refuses to impute them from fewer than 50 samples:

```
Error: This run requires decent allele frequencies, but they aren't being
loaded with --read-freq, and less than 50 samples are available to impute them
from.
```

A 30-sample cohort therefore gets through QC, MOI, LD-prune-skip, ADMIXTURE
(K correctly recovered), IBD, clonality and introgression — 81 of 92 steps —
and then stops at `pca`, taking the PCA figures and everything downstream of
them with it. The remedy PLINK2 itself suggests is `--freq` followed by
`--read-freq`, which is a small change, but it means computing frequencies
from the same small sample PLINK2 is warning you about. That is a judgement
call about what the PCA then means, so it is left open rather than decided
here.

**In short:** under ~50 samples, expect structure (ADMIXTURE), IBD, clonality
and introgression to work with a documented caveat, and PCA to stop the run.


## Known upstream limitation: hmmIBD and long output paths

hmmIBD builds its output filenames in a fixed-size buffer without bounds
checking. Give it a long `-o` prefix and it dies with **SIGTRAP (exit 133,
"Trace/BPT trap: 5") and an empty log** — no error message, nothing to go on.

Measured on `hmmibd 2.1.3` (bioconda, osx-arm64): a 46-character prefix works,
56 fails. The pipeline stays under that by using workspace-relative output
paths, which is why the shipped configs are fine. You can trip it by pointing
`paths.outputs` at a deep absolute directory, or by running the pipeline from
one.

**If Stage 4 dies with exit 133 and an empty `logs/ibd/run_hmmibd_*.log`, this
is why.** Shorten the output path — run from a shallower working directory, or
set `paths.outputs` to something shorter — and it will work unchanged.


## The report

One HTML document that assembles whatever the pipeline actually produced:

```bash
pixi run snakemake --configfile config/config.yaml report --cores 8
```

It lands at `{paths.reports}/report.html`, self-contained — every figure is
embedded, so the single file can be emailed or archived on its own.

It is **role- and config-driven, with no cohort narrative in it**. Titles come
from `cohort.name`; which sections exist is decided by which roles you mapped,
which stages your config switched on, and which output files are actually on
disk. Cluster names, cluster-pair names, focal-group names and model names are
read at render time, never written in.

A stage that did not run is reported as a short note saying *why* and what
would enable it — the same graceful degradation the rules follow. Unset
`metadata.roles.date` and the clonal timeline becomes "no collection dates to
plot against", not a gap or a crash.

Wherever the `{sampleset}` axis applies, the **full and de-clonalized arms are
shown side by side**, and the full-vs-de-clonalized comparison tables get their
own section. That comparison is the point: a claim that looks the same in both
arms survives de-clonalization, and one that changes did not.

`report` is a named target rather than part of `all`, so an existing cohort's
run does not acquire a hard dependency on quarto, and re-rendering after an
edit does not make the pipeline look out of date. Depending on it pulls the
whole pipeline anyway.


## Configuration

Two files under `config/`:

- `config.yaml` — the active config (template; copy or edit in place).
- `cohort.example.yaml` — the V1 Indo cohort, expressed in the agnostic
  schema. Use as a worked template.

`config/schema/config.schema.yaml` is the authoritative schema. The Snakefile
validates the active config against it at parse time and refuses to build a
DAG if validation fails (missing required keys, bad `input_type`, …).

### Role-based metadata

Bring one tidy `samples.tsv`. Map your columns to analytical roles in
`metadata.roles`:

| Role | Required? | If absent |
|---|---|---|
| `sample_id`  | yes | hard error |
| `group`      | no  | downstream group-faceted plots skipped |
| `geography`  | no  | map + province facets skipped |
| `country`    | no  | country facets collapse to "all" |
| `host`       | no  | ignored |
| `date`       | no  | temporal analyses skipped |
| `case_control` | no | selection models reference roles or raw cols |

Optional roles set to `null` (or pointing at a column that doesn't exist) log
a note at validation time and disable their dependent analyses — they never
crash the pipeline.

---

## Tests

```bash
# Everything (config validation, contig derivation, walkthrough parser,
# Stage 5 units + the introgression positive control).
pixi run test

# Individual suites also run standalone, without pytest:
pixi run python tests/test_contigs_from_fai.py     # contigs_from_fai (acceptance §7)
pixi run python tests/test_config_validation.py    # config-schema negative tests
pixi run Rscript tests/R/test_introgression_units.R  # Stage 5 component units
```

(Drop the `pixi run` prefix if you are already inside `pixi shell`.)

`tests/test_introgression.py` is the Stage 5 **positive control**: the tiny
cohort injects a known introgression event and the test asserts the pipeline
recovers exactly it — right window, right samples, nobody else. It builds the
tiny cohort's Stage 5 targets with snakemake if they are missing (seconds when
they are already there), and skips if the toolchain isn't on PATH.

---

## Conventions

- **Don't modify anything outside `agnostic/`.** V1 is frozen.
- **Every rule has a structured docstring** (one-line summary +
  `WHAT/WHY/TUNABLES/OUTPUT/TRY` block). The walkthrough generator in a
  later increment parses these — don't strip them when editing rules.
- **Graceful degradation.** A missing *optional* metadata column disables
  its analysis with a logged note, never a crash.
- **Seam output is shared.** Both WGS and microhap QC paths terminate at
  `outputs/qc/snps.qc.vcf.gz`. Stage 2+ reads from there and is fork-free.

---

## License

MIT — see [`LICENSE`](LICENSE). © 2026 Menzies School of Health Research.
