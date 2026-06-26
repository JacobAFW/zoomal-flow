# Claude Code brief — Stage 2 (MOI / Fws) of the data-agnostic pipeline

Paste into Claude Code from `Pop-gen_pipeline/agnostic/` with the V1 env active
(`source ../envs/activate.sh`). Scope: port V1 Stage 2 (MOI/Fws) into the agnostic tree,
behind the input_type seam, with the country/geography summaries + plots lifted to the
metadata role-map. WGS path implemented and validated against V1; microhap path stubbed.

Pre-work housekeeping (from the Stage-1 review): correct the stale Session-11 line in V1
`MEMORY.md` (~line 481) that still describes `agnostic/data` as a blanket symlink to
`../data` and `cohort.example.yaml` pointing at `../outputs/...` — the curated `data/`
layout replaced both.

---

## Read first

1. `agnostic/DESIGN.md` §3 (role-map + graceful degradation), §4 (the MOI/Fws seam — this
   is **seam 2**), §5 item 3 (geography generalisation).
2. V1 `workflow/rules/02_moi.smk` — the rules you're porting.
3. V1 `scripts/R/{vcf_to_gds.R, run_moimix.R, fws_summary.R, plot_fws_density.R}` — the
   logic. `vcf_to_gds.R` and `run_moimix.R` are already cohort-agnostic (port near-verbatim);
   `fws_summary.R` and `plot_fws_density.R` hardcode `Country`/`State` and an Indonesia-only
   province panel — those get the role treatment (below).
4. How Stage 1 was built (`01_qc.smk` + the `qc_wgs.smk`/`qc_microhap.smk` seam) — mirror
   that exact pattern for Stage 2.

## Hard constraints (same as Stage 1)

- Write only under `agnostic/`. V1 frozen.
- Config-driven, argv-driven; no hardcoded cohort labels/columns/thresholds.
- **Graceful degradation:** a missing optional role disables its summary/plot with a logged
  note, never a crash.
- Every rule keeps the `WHAT/WHY/TUNABLES/OUTPUT/TRY` docstring.
- Reads the canonical role-renamed metadata `outputs/metadata/samples.tsv` (columns are the
  canonical role names: `sample_id`, `country`, `geography`, `group`, `host`, `date`, …) —
  not the raw cohort table. Join Fws to it on `sample_id`.

## 1. Config additions

Add a `moi` block to `config/config.yaml`, `config/cohort.example.yaml`, and the schema
`config/schema/config.schema.yaml` (carry V1's values):

```yaml
moi:
  min_maf: 0.05                 # stricter MAF for moimix (biallelic BAF needs it)
  fws_polyclonal_cutoff: 0.95   # Fws below = polyclonal (Manske 2012 convention)
  fws_exclusion_cutoff: 0.95    # Fws below = excluded from Stage 3+ (single-genotype methods)
  seed: 2023                    # moimix set.seed — match V1 for reproducibility
```

Schema: `moi` object, all four required, `min_maf` number 0–1, cutoffs number 0–1, `seed`
integer. Keep `additionalProperties: true` on the block (later sub-params may land).

## 2. The MOI seam (DESIGN §4, seam 2)

Structure `02_moi.smk` exactly like `01_qc.smk`: common rules, then
`include: "moi_wgs.smk"` or `"moi_microhap.smk"` on `cohort.input_type`.

**Seam contract — `outputs/moi/fws_MOI.tsv`:** two columns, `sample` (matches VCF/`sample_id`)
and `Proportion` (the Fws / within-host-complexity score; *lower = more polyclonal*). Both
input-type paths MUST emit this exact shape so the common rules and Stage 3 are fork-free.
Document the contract at the top of `02_moi.smk`.

**Common rules** (run regardless of input_type; read the seam output):
- `fws_high_moi_list` — `awk 'NR>1 && $2 < {fws_exclusion_cutoff}'` → `outputs/moi/exclude_high_moi.txt`.
  This is the file Stage 3 consumes. Port V1 verbatim.
- `fws_summary` — role-driven (see §4).
- `plot_fws_density` — role-driven (see §4).

**`moi_wgs.smk`** (implemented — port V1):
- `moi_filter_vcf` — `bcftools +fill-tags -- -t AF | bcftools view -e 'MAF<{min_maf}'` on the
  Stage-1 seam `outputs/qc/snps.qc.vcf.gz` → `snps.moi.vcf.gz`. (Biallelic-MAF — WGS-specific,
  hence behind this seam.)
- `vcf_to_gds` — `Rscript vcf_to_gds.R` → `snps.moi.gds` (port verbatim; it's generic).
- `run_moimix` — `Rscript run_moimix.R <gds> <fws_MOI.tsv> <BAF_dataframe.tsv> {seed}` →
  the seam output `fws_MOI.tsv` (+ BAF). Port verbatim; it's generic.

**`moi_microhap.smk`** (stub): one rule producing `outputs/moi/fws_MOI.tsv`, shell body
`echo "MicrohapNotImplemented: must emit fws_MOI.tsv (sample, complexity score) from a
multiallelic-aware MOI estimator — moimix's biallelic-BAF model does not apply. The Fws
cutoff semantics differ for microhap; see DESIGN §4 seam 2." >&2; exit 1`. Comment the
contract.

## 3. (covered above)

## 4. Role-driven summaries + plots (the agnostic lift)

Rewrite the two cohort-specific scripts to read canonical role columns and degrade
gracefully. Resolve which roles exist from the metadata header (the validator already noted
absent roles) — or pass the active role names in as argv from the rule.

- **`fws_summary.R`** → group by the `country` role and the `geography` role instead of
  literal `Country`/`State`. Emit `outputs/moi/fws_by_country.tsv` only if the `country`
  role is present; `outputs/moi/fws_by_geography.tsv` only if `geography` is present. If a
  role is absent, write nothing for it and log `[fws_summary] no <role> column → skipping`.
  Keep the same summary columns (n, median, min, max, n_polyclonal, pct_polyclonal).
- **`plot_fws_density.R`** → panel A: Fws density coloured by the `country` role; panel B:
  coloured by the `geography` role. **Drop the hardcoded `filter(Country == "Indonesia")`
  province panel** — panel B is just "density by geography" across the cohort. Each panel is
  produced only if its role exists. Output names: `reports/figures/fws_density_by_country.{png,svg}`
  and `reports/figures/fws_density_by_geography.{png,svg}`.
- Theme: match whatever you did for `plot_qc_distributions.R` in Stage 1 (you stripped the
  external `theme_pub.R` source there). Be consistent — either port a self-contained theme
  into `agnostic/scripts/R/` or inline a minimal one. No `source("scripts/R/theme_pub.R")`
  with a V1-relative path.

## 5. Snakefile wiring

- `include: "rules/02_moi.smk"`.
- Add Stage 2 targets to `FINAL_TARGETS`: always `fws_MOI.tsv`, `exclude_high_moi.txt`;
  add `fws_by_country.tsv` + `fws_density_by_country.*` only when the `country` role is set,
  and the `geography` equivalents only when `geography` is set (guard with the same
  role-presence check, computed at parse time from `ROLES`).

## 6. Acceptance + validation

Run end-to-end on the Indo example and **validate against V1** `outputs/moi/`:

```bash
snakemake --configfile config/cohort.example.yaml --cores 4 -p
```

- `outputs/moi/fws_MOI.tsv` produced; with `seed: 2023` and the same moimix it should match
  V1's `fws_MOI.tsv` (per-sample Fws). Diff them — flag any sample whose Fws differs beyond
  floating-point noise.
- `outputs/moi/exclude_high_moi.txt` sample count matches V1's.
- `fws_by_country.tsv` `pct_polyclonal` per country matches V1.
- Both density panels render.
- Confirm Stage 3's future input is in place: `exclude_high_moi.txt` exists and is the
  filter Stage 3 will apply.
- Negative-degradation check: temporarily set `metadata.roles.geography: null`, dry-run, and
  confirm the geography summary/plot drop out of the DAG with a logged note and no error.
  (Revert after.)

Report back: the Fws diff result vs V1, the high-MOI exclusion count, and per-country
pct_polyclonal.

## 7. Commit (repo-scaffold discipline)

Stage explicitly — never `git add -A`. Confirm `git status` shows only code/docs (no
`outputs/`, `reports/`, `*.gds`, `*.vcf*`, `samples.tsv`). Commit on `agnostic/main`, no
remote, message summarising the Stage 2 port + the V1 validation numbers.

## Out of scope

Stage 3+ proper, the walkthrough generator, the tiny-cohort fixture. Stop at a green,
V1-validated Stage 2 + commit. The high-MOI exclusion list it produces is the handoff to
Stage 3 (generic cluster labelling), which is the next increment.
