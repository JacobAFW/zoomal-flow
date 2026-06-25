# Claude Code brief — Stage-1 fixup, end-to-end run, and git init

Paste into Claude Code running from `Pop-gen_pipeline/agnostic/` with the V1 env active
(`source ../envs/activate.sh`). Three jobs: (1) fix the missing index dependency, (2) run
Stage 0–1 end-to-end on the real Indo VCF to prove the tools fire, (3) initialise git
cleanly inside `agnostic/`. Plus a small cleanup of artefacts left by the Cowork session.

Context: Stage 0–1 was reviewed in Cowork. The DAG is correct and all unit tests pass, but
no real bcftools/plink/R step has actually run yet (the prior attempt was interrupted at
`subset_nuclear`). Self-containment + repo/bundle isolation were already handled in Cowork
(curated `agnostic/data/`, `agnostic/` excluded from the root `.gitignore` and the Edy
share bundle) — don't redo those.

---

## 0. Cleanup first (artefacts from the Cowork session, which couldn't delete files)

```bash
cd Pop-gen_pipeline/agnostic
rm -f  data_old_symlink         # old blanket data->../data symlink, renamed aside
rm -f  data/_perm_test          # stray probe file parked under gitignored data/
rm -rf .git                     # half-initialised repo with a stuck index.lock — start clean
```

Confirm `data/` is the curated layout: real `data/metadata/samples.tsv` (987 lines) +
`data/vcf` and `data/reference` dir-symlinks into `../../data/`. Leave it as is.

## 1. Fix the missing index dependency (Cowork review item 1)

`rule subset_nuclear` runs `bcftools view -r <contigs>`, which requires the **source VCF to
be indexed**, but the rule doesn't declare the `.csi` as an input. It works today only
because `merged_popgen.vcf.gz.csi` happens to sit beside the VCF; a cohort supplied without
a sitting index fails cryptically, and Snakemake neither tracks nor rebuilds it.

Do both:

- Add an `index_vcf` rule that produces `<vcf>.csi` from the configured `cohort.vcf`
  (`bcftools index -c`), guarded so it no-ops / isn't forced if the index already exists
  and is newer than the VCF. (Snakemake will simply use the existing one as long as it's
  declared.)
- Declare the index as an input to `subset_nuclear`:
  ```python
  input:
      vcf = COHORT["vcf"],
      idx = COHORT["vcf"] + ".csi",
      keep = rules.filter_samples_list.output.keep,
  ```
  Point `idx` at the `index_vcf` output so the dependency is tracked. Keep the
  `WHAT/WHY/TUNABLES/OUTPUT/TRY` docstring convention; note the index requirement in `WHY`.

Mirror the same pattern anywhere else a `-r`/region read assumes an index (check
`mask_regions` — it reads the local `snps.nuclear.vcf.gz`, which `subset_nuclear` already
indexes as an output, so that one is fine).

## 2. Run Stage 0–1 end-to-end on the real VCF (Cowork review item 2)

The point is to prove the actual toolchain produces outputs, not just that the DAG resolves.

```bash
snakemake --configfile config/cohort.example.yaml --cores 4 -p
```

Acceptance:
- `outputs/setup/{vcf_samples.txt, contig_map.tsv, nuclear_contigs.txt}` written;
  `nuclear_contigs.txt` lists the 14 `ordered_PKNH_NN_v2` and no MIT/API.
- `outputs/metadata/samples.tsv` written by `validate_metadata`; its log reports VCF
  concordance and notes any optional-role columns it couldn't find (none expected for the
  Indo example — all roles map).
- `outputs/qc/snps.qc.vcf.gz` (+ `.csi`) produced via the WGS path; `outputs/qc/variant_count.txt`
  holds a plausible integer.
- The four `reports/figures/*_density_*.png` render.
- Capture the final variant count + sample count in the commit message / a short run note.

If anything fails, fix forward — this is the increment's real validation. Sanity-check the
post-QC variant + sample counts against V1's known Stage-1 numbers (see V1
`outputs/qc/variant_count.txt` / logs) and flag any large divergence; small drift from
tool-version differences is expected, an order-of-magnitude gap is not.

## 3. Initialise git cleanly inside `agnostic/` (Cowork item c)

This is a standalone repo for the agnostic pipeline — **separate** from any Indo repo. The
root `.gitignore` already excludes `agnostic/`, so a future root repo won't see it.

Follow the **repo-scaffold** skill's discipline (Jacob's standard for new repos):
default-deny, **stage the INCLUDE set explicitly — never `git add -A`/`git add .`**, and
leave `git init`/commit/remote/push for Jacob to run, not auto-executed. The `.gitignore`
was already regenerated in Cowork from the skill's default (aggressive: all data/genomics/
secret/sample-sheet patterns), so it backstops the explicit staging.

INCLUDE set (code + docs + safe config only):

```bash
cd Pop-gen_pipeline/agnostic
git init
git add .gitignore DESIGN.md README.md \
        briefs/ config/ scripts/ tests/ workflow/
git status            # MUST show only the above. NO data/, outputs/, logs/,
                      # reports/, .snakemake/, __pycache__, samples.tsv, *.fai, *.vcf*
```

If `git status` shows anything from the EXCLUDE set, stop and fix `.gitignore` — do not
commit. Once clean:

```bash
git commit -m "Stage 0-1: data-agnostic pipeline scaffold + setup + QC

WGS path implemented end-to-end; microhap QC stubbed (fails fast on the
seam contract). Contigs auto-derived from the reference .fai; config +
metadata-role schema validated at parse time. Self-contained data/ layout.
Verified: <N> samples, <M> post-QC biallelic SNPs on the Indo example."
```

Fill in `<N>`/`<M>` from the §2 run. Do **not** add a remote or push — Jacob does that.

Both repo-scaffold flags were resolved in Cowork — no action needed:
- `config/cohort.example.yaml` → **include** (confirmed: all 986 sample IDs from
  samples.tsv were grepped against every committable file; none appear — it holds only
  cohort name, input filenames, role→column mappings, thresholds).
- `README.md` → already updated to the curated `data/` layout + carries the skill's
  "what's included / deliberately excluded" section.

The INCLUDE `git add` line above is the final set. Stage exactly those paths.

## Out of scope

Stages 2–6, the walkthrough generator, generic cluster labelling, the synthetic tiny-cohort
VCF. Stop after a green end-to-end Stage 0–1 run and the initial commit. Report back the
run numbers + anything that diverged from the design so it can be reviewed before Stage 2.
