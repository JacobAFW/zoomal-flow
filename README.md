# ZOOMAL-Flow

A data-agnostic population-genetics pipeline for zoonotic malaria. Six analysis stages (QC, MOI/Fws, structure, IBD, introgression, selection) run on any cohort and reference; every cohort-specific assumption lives in config or is derived from the inputs. It was refactored from the V1 Indonesia *P. knowlesi* pipeline into a standalone, portable tool.

The whole pipeline is implemented and runs end to end. The WGS input path is complete; the microhaplotype path is a stubbed seam that fails fast until it is built.

`DESIGN.md` covers the architecture. `docs/` holds the per-topic method specs referenced below.

---

## Start here: run the public example

The repository ships a small, fully public cohort so you can clone it and watch the whole pipeline run before pointing it at your own data.

```bash
# 1. Clone. The repo root is the workspace; there is nothing further to cd into.
git clone https://github.com/JacobAFW/zoomal-flow.git
cd zoomal-flow

# 2. Get pixi if you do not have it. One binary, no root, no system Python; it installs
#    to ~/.pixi/bin (add that to your PATH when it tells you to). Docs at https://pixi.sh
curl -fsSL https://pixi.sh/install.sh | bash

# 3. Build the environment. `install` gets the locked conda stack; `setup` adds the few
#    tools no conda channel publishes (moimix, rehh, and a handful more). `setup` exits
#    non-zero if anything failed to install, so a clean exit means the env is usable.
pixi install
pixi run setup
pixi run check-env      # prints a version for every tool; fails if any is missing

# 4. Get the reference. Download the P. knowlesi A1-H.1 assembly from PlasmoDB into
#    examples/data/reference/. PlasmoDB needs a free login and publishes no stable direct
#    link, so examples/README.md gives the click-path. Then confirm the right assembly:
pixi run bash examples/check_reference.sh

# 5. Run.
pixi run snakemake --configfile examples/config.yaml --cores 4
```

That runs every stage on 71 publicly-archived ENA isolates of *P. knowlesi* from Malaysia (18,279 SNPs, about 4 MB). ADMIXTURE recovers K = 3 by cross-validation, and the 71 samples fall into the same clusters they occupy in the 981-sample cohort they were drawn from. QC, Fws/MOI, structure, IBD, clonality and introgression all run, in both the full and de-clonalized arms.

Read the result as proof that the pipeline runs and your environment is sound, not as a finding about *P. knowlesi*. `examples/README.md` has the accession list, the provenance (Westaway et al., PLOS Negl Trop Dis 2025, doi:10.1371/journal.pntd.0012885), what the example switches off and why, and how its contig names were reconciled with PlasmoDB's.

The example is safe to commit because every sample is a published ENA accession and every metadata field is public. `examples/build_example.sh` re-derives the files and refuses to write them if the sample list shares any ID with the embargoed cohort. `.gitignore` un-ignores those files by name, so default-deny still covers everything around them.

---

## Run it on your own cohort

```bash
# Put your inputs under data/ (gitignored):
#   data/vcf/<cohort>.vcf.gz      bgzipped + indexed
#   data/reference/<ref>.fasta    with a .fai beside it
#   data/metadata/samples.tsv     one tidy row per sample (see Role-based metadata)

cp config/cohort.example.yaml config/config.yaml    # then edit paths + the role map
pixi run snakemake --cores 1 -n                      # validate config + print the DAG
pixi run snakemake --cores 8                         # run
```

To point at a different config: `pixi run snakemake --cores 8 --configfile path/to/my_cohort.yaml`.

`pixi run <cmd>` runs one command inside the environment. `pixi shell` drops you into it so you can call `snakemake` directly. If pixi is not an option at your site, `env/environment.yml` is a conda/mamba fallback (you get a fresh solve rather than the recorded lock). `env/README.md` documents both paths, the tool-version table, and the one deviation from V1 (hmmIBD).

The V1 Indonesia pipeline keeps its own `vvg-box` env under `../envs/`. Nothing here uses or touches it.

---

## What the repository contains, and what it leaves out

Code only: workflow rules, scripts, schemas, docs. It does not track:

- data of any kind (VCFs, reference FASTA/index, PLINK sets, tabular records)
- the sample metadata table (`data/metadata/samples.tsv`), which carries individual-level fields under publication embargo
- run artefacts (`outputs/`, `logs/`, `reports/`, `.snakemake/`)
- credentials, tokens, environment files
- the solved environment (`.pixi/`, about 2.2 GB); only `pixi.toml` and its `pixi.lock` are tracked

`.gitignore` is default-deny across data, genomics, secret and sample-sheet patterns, so a stray input cannot be committed by accident. The public example under `examples/` is the one exception, un-ignored file by file.

---

## Stages

| Stage | File | Status |
|---|---|---|
| 0  setup             | `workflow/rules/00_setup.smk`         | implemented |
| 1  QC (common)       | `workflow/rules/01_qc.smk`            | implemented |
| 1  QC (WGS)          | `workflow/rules/qc_wgs.smk`           | implemented |
| 1  QC (microhap)     | `workflow/rules/qc_microhap.smk`      | stub (fails fast) |
| 2  MOI / Fws         | `workflow/rules/02_moi.smk`           | implemented (WGS path) |
| 3  structure         | `workflow/rules/03_structure.smk`     | implemented |
| 3b structure figures | `workflow/rules/03_figures.smk`       | implemented |
| 4  IBD               | `workflow/rules/04_ibd.smk`           | implemented |
| 4b IBD figures       | `workflow/rules/04_figures.smk`       | implemented |
| 5  introgression     | `workflow/rules/05_introgression.smk` | implemented (read the caveats below) |
| 6  selection (iHS)   | `workflow/rules/06_selection.smk`     | implemented |
| 7  report            | `workflow/rules/99_report.smk`        | implemented |

`docs/WALKTHROUGH.md`, generated from the rule docstrings, is the step-by-step tour of every rule, its tunables, and something to try at each step.

---

## Clonality: every frequency result is produced twice

Clonal samples are near-identical genotypes. Counting each as independent evidence is pseudo-replication, and it biases anything that is fundamentally a frequency or a count: allele frequencies, ADMIXTURE proportions, per-group polyclonality rates, introgression window support, iHS case/control. *Plasmodium* cohorts are frequently clonal, so the pipeline handles this as a first-class step rather than an option.

Four stages (structure, introgression, group summaries, selection) each build twice, under a `{sampleset}` axis:

```
outputs/<stage>/full/     every sample (the operational backbone)
outputs/<stage>/unique/   one representative per clonal group (the comparison arm)
outputs/<stage>/declonalization_comparison.tsv    the two, side by side
```

Both arms build on every run, so you can see per claim how much of it survives collapsing clonal replicates rather than taking the pipeline on faith. Clonality is detected on the full set (you cannot find clonal groups in a set you already thinned), and the full-set clusters stay authoritative for Stage 4 and the Stage-5 cluster definitions. It runs one pass, never iterated.

`docs/clonality.md` has the rationale, the architecture, the representative rule, and how to read the comparison. Turn it off with `clonality.declonalize: false`. To upgrade a cohort run before this existed: `bash scripts/sh/migrate_to_sampleset_layout.sh` (idempotent; `--undo` reverses it).

---

## Introgression (Stage 5): read before you run it

Stage 5 asks, per genomic window, whether a sample carries the genetic signature of a cluster other than the one it is assigned to. It works per cluster pair: for (Kx, Ky) it measures every sample-window's distance to both clusters' consensus alleles, draws a 2D kernel-density cloud per cluster in that distance space, and flags a window when a sample lands in the other cluster's cloud rather than its own.

It answers two questions, at two scales, each with its own null. Do not read one as the other:

- **Cluster level.** Which windows show gene flow between clusters. Controlled by the per-cluster support floor, derived by permutation with an explicit FDR target (`scripts/R/introgression_floor_derivation.R`, spec §9.5). The floor is cohort-specific: derive it, do not inherit a number.
- **Focal level.** Whether a window's introgression is enriched in a subgroup (`focal_group`, e.g. one district) relative to the rest of its own cluster. This has its own size-preserving label-permutation null, scaled to the subgroup, with its own BH/FDR (`scripts/R/introgression_focal_test.R`, spec §9.6). Output: `focal_<group>_enriched_windows.tsv`.

These are independent by design. A subgroup can never clear a floor calibrated on a cluster many times its size, so filtering the focal question through the cluster floor is a mis-posed test rather than a threshold to tune.

The older `unique_windows_in_<focal>_*.tsv` table is descriptive only and is superseded. "Unique to the focal group" is a raw set difference with no null, no p-value, and no multiple-testing control, and one background sample carrying the window flips it out of the set. Cite the enrichment table instead. Full method, design decisions and validation status: `docs/introgression_analysis_spec.md`.

Two things to know before running it:

- **This stage has never been validated against a reference output.** V1's introgression script was written after the HPC run, and its own header records that reference outputs were unavailable. Correctness here rests on the synthetic positive control in `tests/tiny_cohort/`, a known injected introgression event the pipeline must recover exactly, not on reproducing V1. Treat real-cohort numbers as a method result to review, not a validated figure.
- **`pairs: "all"` is quadratic, in two ways.** All-pairs means C(K,2) comparisons: K=3 gives 3, K=6 gives 15, K=8 gives 28, K=10 gives 45. Each is a full density pass, so compute grows fast. Each is also another contrast, so run 45 of them and some windows clear any fixed threshold by chance, and not every cluster pair is a biologically meaningful comparison anyway. The workflow warns above `introgression.pair_warn_threshold`. Set an explicit `introgression.pairs: [[Kx, Ky], ...]` list for the contrasts you actually intend to interpret.

The detection call is pluggable (`introgression.detection_rule`): `absolute` (V1's fixed density cutoffs, the shipped default) or `relative` (deeper in the other cluster's cloud than in its own, so no absolute cutoff drifts when one cluster is more diffuse). To compare the two on data with known ground truth:

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

It reports detection rate against false positives across a threshold grid, and changes no defaults. Picking a rule is a judgement call, and the sweep is the evidence for it.

---

## Small cohorts: what works and what does not yet

The pipeline runs on any cohort, but two of its tools have a hard 50-sample floor because below that they cannot estimate what they need. One case is handled; the other is not yet.

**Handled: LD-pruning (Stage 3).** PLINK2 refuses `--indep-pairwise` below 50 samples. Below `structure.min_samples_for_ld_prune` (default 50, PLINK2's own floor) the pipeline skips the prune, passes the unpruned variant set to ADMIXTURE, and says so loudly: a banner in `logs/structure/{sampleset}/ld_prune.log`, a machine-readable `outputs/structure/{sampleset}/ld_prune_status.txt`, and a warning callout in the report. Read those before quoting the structure results. Linked variants are correlated evidence and can inflate apparent structure, so treat K selection in particular as indicative.

**Not handled: PCA (Stage 3).** PLINK2 `--pca` needs allele frequencies and refuses to impute them from fewer than 50 samples:

```
Error: This run requires decent allele frequencies, but they aren't being
loaded with --read-freq, and less than 50 samples are available to impute them
from.
```

A 30-sample cohort therefore gets through QC, MOI, the LD-prune skip, ADMIXTURE (K recovered correctly), IBD, clonality and introgression (81 of 92 steps) and then stops at `pca`, taking the PCA figures and everything downstream of them with it. PLINK2's own suggested remedy is `--freq` then `--read-freq`, a small change that computes frequencies from the same small sample PLINK2 is warning about. That is a judgement call about what the PCA then means, so it is left open rather than decided in code.

Under about 50 samples, then: structure (ADMIXTURE), IBD, clonality and introgression run with a documented caveat, and PCA stops the run.

---

## Known upstream limitation: hmmIBD and long output paths

hmmIBD builds its output filenames in a fixed-size buffer with no bounds check. Give it a long `-o` prefix and it dies with SIGTRAP (exit 133, "Trace/BPT trap: 5") and an empty log, with no error message to go on.

Measured on `hmmibd 2.1.3` (bioconda, osx-arm64): a 46-character prefix works, 56 fails. The pipeline stays under that with workspace-relative output paths, which is why the shipped configs are fine. You can trip it by pointing `paths.outputs` at a deep absolute directory, or by running from one.

If Stage 4 dies with exit 133 and an empty `logs/ibd/run_hmmibd_*.log`, this is why. Shorten the output path (run from a shallower directory, or set `paths.outputs` to something shorter) and it runs unchanged.

---

## The report

One self-contained HTML document that assembles whatever the pipeline actually produced:

```bash
pixi run snakemake report --configfile config/config.yaml --cores 8
```

The order matters. `report` comes before `--configfile`, because `--configfile` takes one or more files and a target written after it is read as a second config, so the run dies with `FileNotFoundError: 'report'`.

It lands at `{paths.reports}/report.html`, with every figure embedded, so the single file can be emailed or archived on its own. It is role- and config-driven with no cohort narrative written in: titles come from `cohort.name`, and which sections exist is decided by which roles you mapped, which stages your config switched on, and which output files are on disk. Cluster names, cluster-pair names, focal-group names and model names are read at render time.

A stage that did not run appears as a short note saying why and what would enable it, the same graceful degradation the rules follow. Unset `metadata.roles.date` and the clonal timeline becomes "no collection dates to plot against" rather than a gap or a crash. Wherever the `{sampleset}` axis applies, the full and de-clonalized arms sit side by side with their own comparison section, so you can see which claims hold across both arms and which change.

`report` is a named target rather than part of `all`, so an existing run does not acquire a hard dependency on quarto and re-rendering after an edit does not make the pipeline look out of date. Depending on it pulls the whole pipeline anyway.

---

## Configuration

Two files under `config/`:

- `config.yaml` is the active config (a template; copy it or edit in place).
- `cohort.example.yaml` is the V1 Indo cohort expressed in the agnostic schema, a worked template.

`config/schema/config.schema.yaml` is the authoritative schema. The Snakefile validates the active config against it at parse time and refuses to build a DAG on a bad config (missing required keys, an unknown `input_type`, and so on).

### Role-based metadata

Bring one tidy `samples.tsv` and map your columns to analytical roles in `metadata.roles`:

| Role | Required? | If absent |
|---|---|---|
| `sample_id`  | yes | hard error |
| `group`      | no  | group-faceted plots skipped |
| `geography`  | no  | map + region facets skipped |
| `country`    | no  | country facets collapse to "all" |
| `host`       | no  | ignored |
| `date`       | no  | temporal analyses skipped |
| `case_control` | no | selection models fall back to role/raw columns |

An optional role set to `null`, or pointing at a column that does not exist, logs a note at validation time and disables its dependent analyses. It never crashes the pipeline.

---

## Tests

```bash
# Everything: config validation, contig derivation, walkthrough parser,
# Stage 5 units, and the introgression positive control.
pixi run test

# Individual suites also run standalone, without pytest:
pixi run python tests/test_contigs_from_fai.py       # contig derivation
pixi run python tests/test_config_validation.py      # config-schema negative tests
pixi run Rscript tests/R/test_introgression_units.R  # Stage 5 component units
```

Drop the `pixi run` prefix if you are already inside `pixi shell`.

`tests/test_introgression.py` is the Stage 5 positive control: the tiny cohort injects a known introgression event and the test asserts the pipeline recovers exactly it (right window, right samples, nobody else). It builds the tiny cohort's Stage 5 targets with snakemake if they are missing (seconds when they already exist), and skips if the toolchain is not on PATH.

---

## Conventions for contributors

- Do not modify anything outside this repository. The V1 Indonesia pipeline is frozen.
- Every rule carries a structured docstring (a one-line summary plus a `WHAT/WHY/TUNABLES/OUTPUT/TRY` block). `render_walkthrough.py` parses these into `docs/WALKTHROUGH.md`, so keep them intact when editing rules.
- Graceful degradation: a missing optional metadata column disables its analysis with a logged note, never a crash.
- Shared seam output: both the WGS and microhap QC paths terminate at `outputs/qc/snps.qc.vcf.gz`, so Stage 2 onward reads from one place and stays fork-free.

---

## License

MIT, see [`LICENSE`](LICENSE). © 2026 Menzies School of Health Research.
