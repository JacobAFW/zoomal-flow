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

Data lives outside version control. The pipeline expects it under `agnostic/data/`
at the paths the active config points to (see Quick start). `.gitignore` is
default-deny (all data/genomics/secret/sample-sheet patterns) so a stray input
can't be committed by accident.

---

## Quick start

```bash
cd agnostic

# 1. Activate the V1 env (carried over unchanged; not duplicated here)
source ../envs/activate.sh

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
snakemake --cores 1 -n         # validates config; prints the DAG
snakemake --cores 8            # actually runs
```

To run with a different config file:

```bash
snakemake --cores 8 --configfile path/to/my_cohort.yaml
```

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
| Quarto report       | —                                    | later increment |

Stage 0 produces: `outputs/setup/vcf_samples.txt`, `nuclear_contigs.txt`,
`contig_map.tsv`, and a canonical role-renamed `outputs/metadata/samples.tsv`.
Stage 1 (WGS) produces: `outputs/qc/snps.qc.vcf.gz` (the shared seam output),
`outputs/qc/variant_count.txt`, and four QC density panels under
`reports/figures/`.

`docs/WALKTHROUGH.md` (generated from the rule docstrings) is the step-by-step
tour of every rule, its tunables, and something to try at each step.

---

## Introgression (Stage 5) — read before you run it

Stage 5 asks, per genomic window: does a sample carry the genetic signature of
a cluster *other* than the one it is assigned to? It works **per cluster pair**
— for (Kx, Ky) it measures every sample-window's distance to both clusters'
consensus alleles, draws a 2D kernel-density cloud per cluster in that distance
space, and flags a window when a sample lands in the other cluster's cloud
rather than its own. The headline output is the set of windows introgressed
uniquely in a configurable `focal_group`.

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
  --clusters       tests/tiny_cohort/outputs/structure/admix_clusters.tsv \
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
pytest tests/

# Individual suites also run standalone, without pytest:
python tests/test_contigs_from_fai.py     # contigs_from_fai (acceptance §7)
python tests/test_config_validation.py    # config-schema negative tests
Rscript tests/R/test_introgression_units.R  # Stage 5 component units
```

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
