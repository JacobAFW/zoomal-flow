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
| 0  setup       | `workflow/rules/00_setup.smk`     | implemented |
| 1  QC (common) | `workflow/rules/01_qc.smk`        | implemented |
| 1  QC (WGS)    | `workflow/rules/qc_wgs.smk`       | implemented |
| 1  QC (microhap) | `workflow/rules/qc_microhap.smk` | stub (fails fast) |
| 2–6, walkthrough generator, tiny-cohort VCF | — | later increments |

Stage 0 produces: `outputs/setup/vcf_samples.txt`, `nuclear_contigs.txt`,
`contig_map.tsv`, and a canonical role-renamed `outputs/metadata/samples.tsv`.
Stage 1 (WGS) produces: `outputs/qc/snps.qc.vcf.gz` (the shared seam output),
`outputs/qc/variant_count.txt`, and four QC density panels under
`reports/figures/`.

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
# Unit + integration test for contigs_from_fai (acceptance criterion §7).
python tests/test_contigs_from_fai.py

# Negative tests for config schema validation.
python tests/test_config_validation.py
```

Both run without pytest if it isn't installed; they shell out to the helpers
they exercise.

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
