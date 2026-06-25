# Claude Code brief — Stage 0–1 of the data-agnostic pop-gen pipeline

Paste this into Claude Code, run from inside `Pop-gen_pipeline/` (the V1 repo root).
This is the first build increment. Scope: scaffold the agnostic tree + implement
**Stage 0 (setup)** and **Stage 1 (QC)** for the WGS path, with the microhap seam stubbed.

---

## Before you touch anything — read these

1. `agnostic/DESIGN.md` — the approved architecture. Authoritative. Especially §2 (layout),
   §3 (config + metadata schema), §4 (WGS/microhap seams), §5 (generalisation backlog),
   §6 (docstring→walkthrough convention).
2. `workflow/rules/00_setup.smk` and `workflow/rules/01_qc.smk` — the V1 rules you are
   generalising. Replicate their bcftools/plink logic faithfully; only lift the
   cohort-specific assumptions named below.
3. `workflow/config.yaml` — V1 config to carry forward (thresholds, compute, paths blocks).

## Hard constraints

- **Do not modify anything outside `agnostic/`.** V1 is frozen (it's a shipped bundle).
  Read V1 freely; write only under `agnostic/`.
- All new code is config-driven and argv/CLI-driven — no hardcoded paths, contig names,
  thresholds, or cohort labels anywhere in rules or scripts.
- **Graceful degradation is a rule:** a missing *optional* metadata column disables its
  analysis with a logged note, never a crash. Only `sample_id` is hard-required.
- Every rule gets a structured docstring now (convention below) even though the
  walkthrough generator is a later increment — so we don't backfill.
- Use `git` properly: branch `agnostic-stage-00-01`, one commit per logical step, leave it
  unmerged for review.

---

## 1. Scaffold to create (under `agnostic/`)

```
agnostic/
  config/
    config.yaml                 # see §2
    cohort.example.yaml         # the Indo cohort filled in as a worked template
    schema/
      config.schema.yaml        # validate config.yaml at workflow start
      metadata.schema.yaml      # validate samples.tsv (required/optional roles)
  workflow/
    Snakefile                   # loads config, validates it, includes rules, defines `all`
    rules/
      00_setup.smk
      01_qc.smk                 # common QC rules + the input_type include seam
      qc_wgs.smk                # WGS-specific variant filter (implemented)
      qc_microhap.smk           # microhap stub (fails fast, see §5)
  scripts/
    py/
      validate_config.py        # jsonschema/yaml validation, called from Snakefile
      contigs_from_fai.py       # derive contig list + contig→int map from <fasta>.fai
    R/
      validate_metadata.R       # load samples.tsv, check schema + VCF concordance
      plot_qc_distributions.R   # port from V1 scripts/R/plot_qc_distributions.R
  tests/
    tiny_cohort/                # leave empty for now; populated in a later increment
  data/                         # gitignored; config points here
```

Carry `envs/` behaviour from V1 unchanged (don't copy the env, just assume
`source ../envs/activate.sh` for now — note it in `agnostic/README.md`).

## 2. Config (`config/config.yaml`)

Implement exactly the structure in DESIGN.md §3a. The blocks Stage 0–1 actually consume:

```yaml
cohort:
  name: "..."             # output paths / report title
  input_type: "wgs"       # "wgs" | "microhap"  — drives the include seam (§4)
reference:
  fasta: "data/reference/....fasta"   # .fai must sit beside it
  exclude_contigs: ["MIT", "API"]     # regex, matched against .fai contig names
  mask_regions: "data/reference/regions_to_mask.list"  # optional; QC mask skipped if null
metadata:
  table: "data/metadata/samples.tsv"
  roles:
    sample_id: "Sample"   # REQUIRED
    group: "Cluster"
    geography: "State"
    country: "Country"
    host: "Host"
    date: "EnrolDate"
    case_control: null
controls:
  exclude_patterns: ["ctrl", "cpos", "cneg"]
qc:
  min_mac: 2
  filter_pass: true
compute:
  threads_heavy: 8
  threads_light: 2
paths:
  outputs: "outputs"
  logs: "logs"
  reports: "reports"
```

`cohort.example.yaml` = this, filled with the real Indo values (PKA1H1 fasta, the V1 VCF
path, `samples.tsv` from V1's Stage-0 output for reference). The Snakefile validates
`config.yaml` against `schema/config.schema.yaml` on load and **errors with a clear message**
if a required key is missing or `input_type` isn't one of the two allowed values.

## 3. Metadata schema + Stage 0 setup

V1's six-file `build_metadata.R` is **dropped** — do not port it. The agnostic input is one
tidy `samples.tsv`. Stage 0 becomes:

- `rule extract_vcf_samples` — `bcftools query -l {vcf}` → `outputs/setup/vcf_samples.txt`
  (unchanged from V1).
- `rule validate_metadata` — run `scripts/R/validate_metadata.R`:
  - load `metadata.table`; confirm the `sample_id` role column exists (hard error if not);
  - check every other declared role maps to a real column, else log "role X column 'Y'
    absent → dependent analyses will be skipped" and continue;
  - concordance vs `vcf_samples.txt`: report samples in VCF but not the table, and vice
    versa, to the log; write the validated, role-renamed table to
    `outputs/metadata/samples.tsv` (canonical role names as columns).

`metadata.schema.yaml` encodes required (`sample_id`) vs optional (all others) per
DESIGN.md §3b, with the column→role mapping read from config.

## 4. Stage 1 QC — common rules (port from V1 `01_qc.smk`)

Replicate faithfully, with the one generalisation noted:

- `rule filter_samples_list` — drop controls via `grep -viE` on `controls.exclude_patterns`
  (V1 logic, verbatim).
- `rule subset_nuclear` — **generalised:** contigs are no longer hand-listed. Call
  `scripts/py/contigs_from_fai.py {fasta}.fai --exclude <patterns>` to produce the contig
  list at workflow-parse time (a checkpoint or a plain Python helper imported in the
  Snakefile, your call — simplest that works). Then `bcftools view -S keep --force-samples
  -r <contigs> -Oz` exactly as V1. `contigs_from_fai.py` also emits the contig→integer map
  (col1 order) to `outputs/setup/contig_map.tsv` for later stages.
- `rule prepare_mask_bed` — V1's awk colon-dash→BED, **only if** `reference.mask_regions`
  is set; otherwise skip masking entirely (graceful degradation).
- `rule mask_regions` — `bcftools view -T ^{bed}` (V1 verbatim); bypassed if no mask.

## 5. Stage 1 QC — the WGS/microhap seam (DESIGN §4)

In `01_qc.smk`, after the mask step, branch on `config["cohort"]["input_type"]`:

```python
if config["cohort"]["input_type"] == "wgs":
    include: "qc_wgs.smk"
elif config["cohort"]["input_type"] == "microhap":
    include: "qc_microhap.smk"
```

- `qc_wgs.smk` — implement the V1 final filter:
  - `rule biallelic_pass_mac` — `bcftools view -m2 -M2 -v snps -f PASS | bcftools view -e
    'MAC<{min_mac}'` → `outputs/qc/snps.biallelic.vcf.gz` (V1 verbatim).
  - `rule count_biallelic_variants` — single-line variant count (V1 verbatim).
  - `rule compute_maf_missingness` — `plink --freq --missing --double-id --allow-extra-chr`
    (V1 verbatim).
  - `rule plot_qc_distributions` — port `scripts/R/plot_qc_distributions.R` unchanged.
- `qc_microhap.smk` — **stub only.** One rule producing the same downstream output name
  (`outputs/qc/snps.biallelic.vcf.gz` — or a neutral `snps.qc.vcf.gz`; pick one name both
  paths share so Stage 2+ is fork-free) whose shell body is:
  `echo "MicrohapNotImplemented: this rule must emit a QC'd VCF + loci.tsv describing
  multiallelic allele structure (see DESIGN.md §4). Microhap path not yet built." >&2; exit 1`
  Add a comment block stating the contract (locus-level call-rate QC, no -m2 -M2).

The shared output filename is the seam contract: whichever path runs, Stage 2 reads the
same file. Decide the shared name once and document it at the top of `01_qc.smk`.

## 6. Docstring convention (apply to every rule now)

```python
rule subset_nuclear:
    """
    Restrict the VCF to nuclear contigs + non-control samples.

    WHAT: bcftools view -S keep -r <contigs> ; contigs derived from <fasta>.fai
          minus reference.exclude_contigs.
    WHY:  MIT/API are haploid and break diploid-SNP assumptions downstream.
    TUNABLES: reference.exclude_contigs, controls.exclude_patterns
    OUTPUT: outputs/qc/snps.nuclear.vcf.gz
    TRY:  widen exclude_contigs to drop a noisy contig and re-run from here to
          see the effect on the Stage-1 variant count.
    """
```

The `WHAT/WHY/TUNABLES/OUTPUT/TRY` block is what the later walkthrough generator parses.
`TRY` is the "things to play with" hook (DESIGN §6) — include a sensible one per rule.

## 7. Acceptance criteria

- `snakemake -n` (dry run) from `agnostic/` builds a clean DAG for `input_type: wgs` and
  reports the rules in order, with no reference to hand-listed contigs anywhere.
- Switching `input_type: microhap` and dry-running shows the stub rule wired in; actually
  running it fails fast with the `MicrohapNotImplemented` message.
- Config validation rejects a config missing `cohort.name` / bad `input_type` with a clear
  error (add a quick negative test).
- `contigs_from_fai.py` unit-tested against V1's `.fai`: it must reproduce the 14
  `ordered_PKNH_NN_v2` names and drop MIT/API.
- A short `agnostic/README.md` quick-start (env activate, edit config, `snakemake --cores N`).
- Branch `agnostic-stage-00-01` committed, unmerged, ready for review.

## Out of scope (later increments)

Stages 2–6, the walkthrough generator script, the synthetic tiny-cohort VCF, generic
cluster labelling. Don't start them. Stop at a green Stage 0–1 WGS dry run + real run on
the Indo VCF if the env is active.
