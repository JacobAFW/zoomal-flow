# ZOOMAL-Flow (data-agnostic pop-gen pipeline) — design doc

Draft for review, 25 Jun 2026. Companion to `agnostic_pipeline_kickoff_250626.md`.
**Nothing is built yet.** All four §7 questions resolved (see §7); architecture is settled
pending your go-ahead on the build sequence in §8.
Grounded in a full read of the V1 `Pop-gen_pipeline` (config + all eight
`workflow/rules/*.smk` + the cohort-specific R scripts), which stays frozen.

---

## 0. Decisions locked this session

From the kickoff back-and-forth:

1. **Microhap = design the seams now, implement WGS only.** Branch points are built
   into the architecture as explicit, named extension points. The WGS path is the only
   one wired and tested; the microhap path is stubbed with a clear contract at each seam.
   No microhap data is validated against this session.
2. **One source of truth = generate the manual MD from the rules.** Each rule carries a
   structured docstring; an extraction script renders the manual step-through MD. The
   rules are authoritative; the MD is always regenerated, never hand-edited. V1's
   hand-maintained `WALKTHROUGH.md` drift problem cannot recur.
3. **This session produces this design doc only.** Build starts after you approve it.

Proceeding-unless-you-object defaults (settled from reading the code):

- Nested dir at `Pop-gen_pipeline/agnostic/`, self-contained for later lift into its own repo.
- V1 `Pop-gen_pipeline` is not refactored in place; the Edy bundle stays frozen.
- Nuclear contig list derived from the reference `.fai`, not hand-listed in config.
- Metadata schema starts from V1's 8-column `samples.tsv` with a required-vs-optional
  split and a config block declaring which columns drive the auto-run.

---

## 1. What "data-agnostic" means here

V1 answers one cohort's question (990 Malaysian + Indonesian *P. knowlesi*, PKA1H1
reference, slides 6–10). The agnostic version runs the **same six analysis stages** for
*any* cohort and reference, with every cohort-specific assumption lifted to config or
derived from inputs. The science is unchanged; the wiring is generalised.

Three classes of cohort-specificity to remove (all enumerated in §5):

- **Hard-coded biology of the cohort** — contig names, reference labels, geography,
  control-naming patterns, clonal-sample mappings.
- **Hard-coded shape of the metadata** — six bespoke source files joined by a
  cohort-specific `build_metadata.R`.
- **Hard-coded input assumptions** — biallelic-SNP WGS baked into QC, MOI/Fws, and
  structure-input prep (the microhap fork).

---

## 2. Directory layout

```
Pop-gen_pipeline/agnostic/          # self-contained; lift-to-repo ready
  DESIGN.md                         # this doc
  README.md                         # quick-start (generated section + hand intro)
  config/
    config.yaml                     # analysis params + paths + auto-run knobs
    cohort.example.yaml             # worked example (the Indo cohort, as a template)
    schema/
      config.schema.yaml            # validates config.yaml at workflow start
      metadata.schema.yaml          # validates the sample table (required/optional cols)
  workflow/
    Snakefile
    rules/                          # per-stage .smk, same stage numbering as V1
      00_setup.smk … 06_selection.smk, 99_report.smk
    docstrings/                     # (mechanism, §6) — not separate files; see note
  scripts/
    R/                              # argv/CLI-driven, cohort-agnostic
    py/                             # docstring→MD extractor, schema validators
    sh/
  docs/
    WALKTHROUGH.md                  # GENERATED from rule docstrings — do not edit
  data/                            # cohort inputs (gitignored); pointed at by config
  tests/
    tiny_cohort/                    # small synthetic WGS VCF for CI / smoke tests
  envs/                             # carried from V1 (pixi vvg-box), unchanged
```

The `agnostic/` tree never imports from `../` at runtime — paths resolve inside it or
to config-declared data locations, so it lifts cleanly.

---

## 3. Two-layer configuration

V1 already config-drives paths + thresholds. The agnostic version keeps that and adds
two things: an **input-type declaration** (drives the WGS/microhap fork) and an
**auto-run column map** (tells the pipeline which metadata columns mean what).

### 3a. `config.yaml` — additions over V1

```yaml
cohort:
  name: "indonesia_pk"          # used in output paths / report title
  input_type: "wgs"             # "wgs" | "microhap"  ← the fork switch (§4)

reference:
  fasta: "data/reference/strain_A1_H.1.Icor.fasta"
  # nuclear_contigs no longer hand-listed: derived from <fasta>.fai, minus `exclude_contigs`
  exclude_contigs: ["MIT", "API"]   # matched against .fai names (regex)
  gff: "data/reference/PknowlesiA1H1.gff"   # optional; introgression gene lookups

metadata:
  table: "data/metadata/samples.tsv"   # ONE pre-built tidy table — the only supported input (§3b)
  # Column-role map: which columns drive which analysis. This is the heart of agnostic.
  roles:
    sample_id:    "Sample"     # must match VCF sample names
    group:        "Cluster"    # structure / IBD grouping column
    geography:    "State"      # province/region facet for plots + maps
    country:      "Country"    # top-level facet
    host:         "Host"
    date:         "EnrolDate"  # temporal analyses; null disables them
    case_control: null         # optional; selection models reference roles or raw cols

structure:
  cluster_labelling: "auto"    # "auto" (majority-geography) | "numbered" | "reference"
  reference_labels: null       # used only when labelling == "reference" (V1's Mn/Mf/Pen)

controls:
  exclude_patterns: ["ctrl", "cpos", "cneg"]   # already config in V1; kept
```

Existing V1 blocks (`qc`, `moi`, `structure` thresholds, `ibd`, `introgression`,
`selection`, `compute`, `disk`) carry over largely unchanged.

### 3b. Metadata schema — required vs optional

The agnostic version reads **one tidy sample table** (`samples.tsv`) — and nothing else.
V1's runtime six-file assembly does **not** carry over: users are expected to arrive with
clean data, and a convoluted "here's how we stitched messy sources together" builder adds
no value to a cohort that isn't ours. The Indo `build_metadata.R` stays in frozen V1 only.
The agnostic contract is simply: supply a table that matches the schema below.

Schema (validated at workflow start against `metadata.schema.yaml`):

| Column role | Required? | Drives | If absent |
|---|---|---|---|
| `sample_id` | **required** | everything; must match VCF | hard error |
| `group` | optional* | structure colouring, IBD cluster loop, selection | falls back to ADMIXTURE-derived clusters |
| `country` | optional | top-level facet (Fws, ADMIXTURE) | facet collapses to "all" |
| `geography` | optional | province facets, maps, connectivity colour | map + geography facets skipped |
| `host` | optional | host filtering | ignored |
| `date` | optional | temporal/clonal plots | temporal stage skipped, not errored |
| `case_control` | optional | selection models | model must define filters inline |

\* `group` is "required for the grouped views to render, optional for the pipeline to
run" — absence degrades gracefully rather than crashing. **Graceful degradation is a
design rule: a missing optional column disables its analysis with a logged note, never a
crash.** This is what lets a cohort with only `sample_id` + `country` still get QC, MOI,
ADMIXTURE, and PCA.

---

## 4. WGS vs microhap branch-point map

`cohort.input_type` selects the path. "Both are VCF" is necessary, not sufficient —
multiallelic microhap loci violate the biallelic-SNP assumption in three specific
places. Each is an explicit seam: WGS implemented, microhap stubbed with a contract.

| Stage | V1 rule(s) | WGS path (implemented) | Microhap seam (stubbed) |
|---|---|---|---|
| **QC** | `01_qc.smk` `biallelic_pass_mac` (`-m2 -M2 -v snps`) | biallelic SNP filter as today | keep multiallelic loci; locus-level call-rate + per-locus allele-count QC instead of `-m2 -M2`. Contract: emit a QC'd VCF + a `loci.tsv` describing allele structure. |
| **MOI / Fws** | `02_moi.smk` `run_moimix` (moimix BAF from biallelic depth ratios) | moimix Fws as today | moimix's biallelic-BAF model doesn't apply. Seam calls an MOI estimator that takes multiallelic genotype counts (e.g. heterozygosity- or `THE REAL McCOIL`-style). Contract: emit `fws_MOI.tsv`-shaped table (`sample`, complexity score) so Stage 3 exclusion is identical downstream. |
| **Structure input prep** | `03_structure.smk` `normalise_vcf` → `vcf_to_plink` → ADMIXTURE/PLINK | split-multiallelic + biallelic PLINK as today | ADMIXTURE/PLINK need biallelic. Seam: either (a) decompose microhaps to biallelic SNPs, or (b) swap to a distance/PCA method that accepts multiallelic loci. Contract: emit the cleaned bfile **or** a distance matrix that Stage 3's PCA/NJT/ADMIXTURE-or-substitute consumes. |

Stages 4–6 (IBD, introgression, selection) consume Stage-3 outputs and the genotype
table; once the three seams above produce schema-compatible outputs, downstream stages
need no fork. The seams are encoded as Snakemake rule selection keyed on `input_type`
(e.g. `01_qc.smk` includes `qc_wgs.smk` *or* `qc_microhap.smk`), so the microhap path is
a future file drop, not a pipeline rewrite.

**Stub contract pattern:** each microhap rule exists as a rule that fails fast with a
clear `MicrohapNotImplemented: <what it must produce>` message naming its output schema.
This makes the seam visible and testable without implementing the biology.

---

## 5. Generalisation backlog (V1 cohort-assumptions → agnostic)

Each item names the exact V1 location and the lift. Ordered by leverage.

1. **Contig names** — V1 `config.yaml` hand-lists 14 `ordered_PKNH_NN_v2`; `03_structure`
   `make_chrom_update` and `genotype_table.R` string-strip `ordered_PKNH_`/`_v2`.
   → Derive contigs from `<fasta>.fai` minus `reference.exclude_contigs`; replace the
   hard-coded `gsub` with a contig→integer map built from the `.fai` order.
2. **Cluster labelling** — `assign_clusters.R` majority-votes ADMIXTURE components onto
   the literal labels `Mn`/`Mf`/`Peninsular` from `Pk_clusters_metadata.csv`, hard-codes
   `K=3`, and the IBD loop hard-codes `IBD_CLUSTERS = ["Mf","Mn","Peninsular"]`.
   → `structure.cluster_labelling`: `numbered` (Cluster_1…K), `auto` (label by majority
   `geography` role), or `reference` (V1 behaviour via `reference_labels`). Best-K chosen
   from CV (or `admixture_k` override) rather than fixed 3. IBD loop reads cluster names
   from the assignment TSV, not a literal list.
3. **Geography** — province map + `State` categories are Malaysian/Indonesian
   (`plot_province_map.R`; `STATE_PAL` / `PROVINCE_PAL_MAP` in `palettes.R`).
   → Map + facets driven by the `geography` role; palette generated for whatever levels
   exist; map rule skipped if no `geography` column / no GIS coords.
4. **Stage-4 specifics** — `plot_clonal_temporal.R` hard-codes an 8-sample clonal mapping
   + Sabang dates read from legacy sources; `plot_ibd_connectivity.R` colours by `State`.
   → Clonal/temporal driven by the `date` role + computed clonal pairs (no hand mapping);
   connectivity colours by the `group`/`geography` role. Whole stage skipped if `date` absent.
5. **Control naming** — already config (`qc.exclude_patterns`) in V1; carry over verbatim.
6. **Selection models** — V1 `selection.models` filters reference literal `Cluster`/`State`
   values. → Keep the inline-filter mechanism but document that filters reference role
   columns; ship the Aceh example as `cohort.example.yaml`, not a default in core config.
7. **GIS coords** — V1 borrows `data/validation/admix_clusters_gis.tsv`.
   → Optional `metadata.gis` file (sample_id, lat, long); map skipped if absent.

Already config-driven in V1 and carried over as-is: VCF path, all thresholds, the single
metadata-table + VCF-matching join.

---

## 6. Source of truth: rules → manual MD

Mechanism (decision 2):

- Every rule's docstring uses a small structured convention — a one-line summary, a
  `WHAT`/`WHY`/`TUNABLES`/`OUTPUT` block — that reads fine as a Python docstring *and*
  parses cleanly.
- `scripts/py/render_walkthrough.py` walks `workflow/rules/*.smk` in stage order, pulls
  each docstring, resolves the config values it references, and writes
  `docs/WALKTHROUGH.md` with a header stamping the commit + config it was generated from.
- `docs/WALKTHROUGH.md` is **committed as raw Markdown** — not knitted/rendered to HTML.
  The raw MD *is* the walkthrough: a plain step-by-step a user reads to play around at any
  step. A CI check (or a `snakemake walkthrough` target diff) fails if the committed copy
  is stale relative to the rules, so the commit can't drift.
- The "manual step-through so I can change things along the way" need is met by the
  walkthrough listing, per stage: the exact `snakemake <target>` to run that stage alone,
  the config knobs in play, **and a "Things to try" block** of suggested experiments. E.g.
  for Stage 3: *"Re-run ADMIXTURE at a different K than the CV-auto pick and check the
  `admixture_bars_by_province` figure — does the extra component map onto a real geographic
  region, or is it noise? If it tracks geography cleanly, that K may be telling you
  something the CV minimum smoothed over."* Each stage gets a couple of these, written
  into the docstring so they regenerate with the rest.

This replaces a hand-maintained MD with a rendered one; there is exactly one source.

---

## 7. Resolved decisions (settled 25 Jun 2026)

1. **Best-K selection → auto-pick from CV error.** ADMIXTURE best-K is the CV-error
   minimum, computed automatically; the CV plot is an automatic pipeline output. Users who
   want to play pick a different K by hand off that plot (the manual walkthrough shows how
   and suggests it as an experiment, §6). No fixed K=3, no required override knob — though
   the manual path can pass a K explicitly.
2. **`WALKTHROUGH.md` → committed, raw MD, not knitted.** The plain Markdown is the
   walkthrough; it's committed so it's viewable directly. A staleness check keeps the
   commit honest. It carries "things to try" suggestions at each step (§6).
3. **Metadata builder → dropped from agnostic.** No six-file assembly in the agnostic tree.
   Users bring clean data: one tidy `samples.tsv` matching the schema (§3b). The messy-data
   stitching script lives only in frozen V1.
4. **Tiny test cohort → yes.** A small synthetic WGS VCF (a few contigs, a few dozen
   samples) under `tests/tiny_cohort/` gives the WGS path a fast smoke test independent of
   the real ~18 GB VCF.

---

## 8. Proposed build sequence (after approval)

1. Scaffold `agnostic/` + `config/schema/*` + the two-layer config with validation.
2. Port Stage 0–1 (setup + QC) with contig auto-detect from `.fai` and the WGS/microhap
   `include:` seam (microhap stubs failing fast).
3. Generic cluster labelling (Stage 3 `assign_clusters` rewrite) — highest-leverage lift.
4. Geography/`group`-role-driven plots (Stages 3–4) with graceful degradation.
5. `render_walkthrough.py` + docstring convention across all ported rules.
6. Tiny-cohort smoke test green end-to-end on the WGS path.

Stages are ported in V1's existing order; each is mergeable on its own.
