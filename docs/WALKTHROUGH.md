# ZOOMAL-Flow — Walkthrough

*Generated from `workflow/rules/*.smk` docstrings by `scripts/py/render_walkthrough.py`.*  
*Do NOT hand-edit — regenerate via `snakemake walkthrough` or `python scripts/py/render_walkthrough.py --write`.*

- Config: `config/cohort.example.yaml`
- Commit: `f1b3b7d-dirty`

Each rule below carries a WHAT/WHY block, its resolved TUNABLES (current values from the config above), its OUTPUT path(s), and a TRY suggestion — a concrete experiment you can run by editing the config and re-invoking that stage's target.

---

## Stage 0 — Setup (VCF sample list, contig map, metadata validation)

Source: `workflow/rules/00_setup.smk`

### `index_vcf`

Build the .csi index for the cohort VCF if it isn't already there.

**WHAT.** bcftools index -c <vcf>

**WHY.** bcftools view -r <region> (used by subset_nuclear) requires the
input VCF to be indexed. Without this rule, a cohort supplied
without a sitting .csi fails cryptically at subset_nuclear. With
it, the index is a tracked dependency that Snakemake will rebuild
if the VCF mtime moves past the index's.

**TUNABLES.**

- `cohort.vcf` = `'data/vcf/merged_popgen.clean.vcf.gz'`

**OUTPUT.** `<cohort.vcf>.csi`

**TRY.** `touch -d "1 hour ago" <cohort.vcf>` then re-run — the rule
re-builds the index because its declared input is now newer
than its declared output.

---

### `extract_vcf_samples`

Extract one VCF sample ID per line; downstream rules read this rather
than re-opening the (multi-GB) VCF.

**WHAT.** bcftools query -l <vcf>

**WHY.** Every QC/metadata step needs the sample roster; doing it once,
here, avoids redundant VCF reads.

**TUNABLES.**

- `cohort.vcf` = `'data/vcf/merged_popgen.clean.vcf.gz'`

**OUTPUT.** `{outputs}/setup/vcf_samples.txt`

**TRY.** swap cohort.vcf to a subsetted VCF (bcftools view -s s1,s2,...)
to dry-run the pipeline on a sample subset without re-deriving anything.

---

### `build_contig_map`

Persist the nuclear-contig list and contig→integer map derived from the
reference .fai (minus reference.exclude_contigs).

**WHAT.** scripts/py/contigs_from_fai.py <fasta>.fai --exclude <patterns>

**WHY.** Replaces V1's hand-listed nuclear_contigs array. The integer map
is the input later stages need for PLINK chromosome codes.

**TUNABLES.**

- `reference.fasta` = `'data/reference/strain_A1_H.1.Icor.fasta'`
- `reference.exclude_contigs` = `['MIT', 'API']`

**OUTPUT.** `{outputs}/setup/nuclear_contigs.txt, {outputs}/setup/contig_map.tsv`

**TRY.** add a contig regex to reference.exclude_contigs (e.g. "_14_") and
re-run from here to see Stage 1 drop that contig from the VCF.

---

### `validate_metadata`

Validate samples.tsv against the configured role map + VCF concordance.
Emit a canonical, role-renamed copy at {outputs}/metadata/samples.tsv
so downstream rules read by role (sample_id, group, geography, ...) and
never reach back into the cohort-specific column names.

**WHAT.** scripts/R/validate_metadata.R <samples.tsv> <vcf_samples.txt>
<roles.json> <out_tsv> <out_log>

**WHY.** Hard-errors only on missing sample_id; missing optional roles log
a note and let dependent analyses skip gracefully (DESIGN §3b).

**TUNABLES.**

- `metadata.table` = `'data/metadata/samples.tsv'`
- `metadata.roles` = `{'sample_id': 'Sample', 'group': 'Cluster', 'geography': 'State', 'country': 'Country', 'host': 'Host', 'date': 'EnrolDate', 'case_control': None}`

**OUTPUT.** `{outputs}/metadata/samples.tsv`

**TRY.** null out metadata.roles.geography in your config — the rule logs
"geography role skipped" and downstream geography facets vanish
without a crash.

---

## Stage 1 — QC (common)

Source: `workflow/rules/01_qc.smk`

### `filter_samples_list`

Build a keep-list by dropping anything in the VCF roster matching the
controls.exclude_patterns regex (ctrl|cpos|cneg by default).

**WHAT.** grep -viE '<patterns>' vcf_samples.txt

**WHY.** Controls in a pop-gen VCF break allele-frequency stats; dropping
them up front means every downstream count is over the analysis
set.

**TUNABLES.**

- `controls.exclude_patterns` = `['ctrl', 'cpos', 'cneg']`

**OUTPUT.** `{outputs}/qc/keep_samples.txt`

**TRY.** add a cohort-specific pattern (e.g. "blank") to the list and
re-run from here to see the sample count drop accordingly.

---

### `subset_nuclear`

Restrict the VCF to nuclear contigs + non-control samples.

**WHAT.** bcftools view -S keep -r <contigs> ; contigs derived from
<fasta>.fai minus reference.exclude_contigs at parse time.

**WHY.** Mitochondrial / apicoplast contigs are haploid and break
diploid-SNP assumptions in MOI / structure / IBD downstream.
bcftools `-r` requires a .csi/tbi beside the VCF; that index is
declared as an explicit input below (built by rule index_vcf)
so Snakemake tracks it and rebuilds it when stale, rather than
relying on one being silently sitting on disk.

**TUNABLES.**

- `reference.exclude_contigs` = `['MIT', 'API']`
- `controls.exclude_patterns` = `['ctrl', 'cpos', 'cneg']`

**OUTPUT.** `{outputs}/qc/snps.nuclear.vcf.gz`

**TRY.** widen reference.exclude_contigs to drop a noisy contig and
re-run from here to see the effect on the Stage-1 variant count.

---

## Stage 1 — QC (WGS path)

Source: `workflow/rules/qc_wgs.smk`

### `biallelic_pass_mac`

Keep biallelic SNPs with FILTER=PASS and MAC ≥ qc.min_mac. Writes the
Stage-1 seam output that Stage 2+ reads (snps.qc.vcf.gz).

**WHAT.** bcftools view -m2 -M2 -v snps -f PASS | bcftools view -e 'MAC<N'

**WHY.** Downstream stages (MOI/Fws, ADMIXTURE, PCA) assume biallelic
diploid SNPs. MAC ≥ 2 removes singletons (mostly sequencing error).

**TUNABLES.**

- `qc.min_mac` = `2`
- `qc.filter_pass` = `True`

**OUTPUT.** `{outputs}/qc/snps.qc.vcf.gz   (the seam — see top of 01_qc.smk)`

**TRY.** bump qc.min_mac to 5 to look at only common variants; or set
qc.filter_pass: false to keep non-PASS calls and re-run from
here to see how much sequencing noise the PASS filter removes.

---

### `count_biallelic_variants`

Single-line variant count for the Stage-1 seam output.

**WHAT.** bcftools view -H | wc -l

**WHY.** The report renders the post-Stage-1 count without re-running
bcftools at knit time.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/qc/variant_count.txt`

**TRY.** after changing qc.min_mac, diff this number against an earlier
run to see the impact of the MAC cutoff in one integer.

---

### `compute_maf_missingness`

Compute MAF + per-sample/per-variant missingness with PLINK 1.9.

**WHAT.** plink --vcf <qc.vcf.gz> --double-id --allow-extra-chr --freq --missing

**WHY.** Cheap diagnostic numbers + density panels that surface problem
contigs or low-coverage samples before later stages run.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/qc/plink/Pk.{{frq,imiss,lmiss}}`

**TRY.** after a config change, eyeball maf_density_overall.png — a kink
on the low end usually means the MAC cutoff should be higher.

---

### `plot_qc_distributions`

Render the four QC density panels (MAF overall + by chrom; missingness
per sample + per variant) used by the Quarto chapter.

**WHAT.** scripts/R/plot_qc_distributions.R <frq> <imiss> <lmiss> <4 png paths>

**WHY.** Lets a reader eyeball the QC distributions before trusting any
downstream number.

**TUNABLES.** *(none)*

**OUTPUT.** `{reports}/figures/maf_density_overall.png + 3 others`

**TRY.** open maf_density_by_chrom.png — if one contig sits noticeably
higher than the rest, that's a candidate for reference.exclude_contigs.

---

## Stage 1 — QC (microhap seam stub)

Source: `workflow/rules/qc_microhap.smk`

### `microhap_qc_stub`

STUB. Produces the same seam output as qc_wgs.smk's terminal rule,
so Stage 2+ wiring is identical regardless of input_type. Running it
fails fast with the contract message.

**WHAT.** echo MicrohapNotImplemented >&2 ; exit 1

**WHY.** Forces explicit work to implement the microhap path; refuses to
let a half-built path produce silently-wrong output.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/qc/snps.qc.vcf.gz  (never actually written)`

**TRY.** swap cohort.input_type to "wgs" to run the real WGS path; this
rule disappears from the DAG.

---

## Stage 2 — MOI / Fws (common)

Source: `workflow/rules/02_moi.smk`

### `fws_high_moi_list`

Samples with Fws < moi.fws_exclusion_cutoff are polyclonal and violate
the single-genotype assumption of ADMIXTURE / PCA / hmmIBD in Stage 3+.
This list is the Stage 2 → Stage 3 handoff.

**WHAT.** awk 'NR>1 && $2 < <cutoff>' fws_MOI.tsv → exclude list

**WHY.** Single-genotype methods downstream silently produce wrong answers
on polyclonal samples. Better to filter them up front than rely on
each downstream rule to remember the threshold.

**TUNABLES.**

- `moi.fws_exclusion_cutoff` = `0.95`

**OUTPUT.** `{outputs}/moi/exclude_high_moi.txt`

**TRY.** raise fws_exclusion_cutoff to 0.85 (V1's HPC value) to keep
mildly-polyclonal samples in the Stage-3 input and see how it
shifts ADMIXTURE cluster counts.

---

## Stage 2 — MOI / Fws (WGS path)

Source: `workflow/rules/moi_wgs.smk`

### `moi_filter_vcf`

Stricter MAF filter on the Stage-1 seam VCF before moimix.

**WHAT.** bcftools +fill-tags -- -t AF | bcftools view -e 'MAF<min_maf'

**WHY.** moimix estimates BAF from read-depth ratios; low-MAF sites have
unstable B-allele frequencies and add noise. WGS-specific (the
microhap path doesn't use biallelic MAF), hence behind this seam.

**TUNABLES.**

- `moi.min_maf` = `0.05`

**OUTPUT.** `{outputs}/moi/snps.moi.vcf.gz`

**TRY.** raise moi.min_maf from 0.05 to 0.10 to keep only well-sampled
variants; usually trims a small % of the Fws tail.

---

### `vcf_to_gds`

Convert the moimix-filtered VCF to SeqArray GDS.

**WHAT.** SeqArray::seqVCF2GDS(storage.option = "LZ4_RA")

**WHY.** GDS = random-access by variant + sample. moimix's BAF/Fws calls
do per-sample passes over the variant matrix; tabix VCF is
orders of magnitude slower for that pattern.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/moi/snps.moi.gds`

**TRY.** run `Rscript scripts/R/vcf_to_gds.R <vcf> <out>` standalone on
a subset VCF to see GDS conversion timing per million variants.

---

### `run_moimix`

Compute per-sample BAF matrix + Fws via moimix. Writes the seam output.

**WHAT.** bafMatrix(gds) + getFws(gds), with set.seed(moi.seed) for
reproducibility.

**WHY.** Fws is the within-host complexity score — lower = more polyclonal.
The seam contract requires this output to be (sample, Proportion)
so Stage 3+ doesn't branch on input_type.

**TUNABLES.**

- `moi.seed` = `2023`

**OUTPUT.** `{outputs}/moi/fws_MOI.tsv         (the seam)
{outputs}/moi/BAF_dataframe.tsv   (long-format BAF, for the
Quarto narrative / QC eyes)`

**TRY.** change moi.seed and compare fws_MOI.tsv pre/post — should be
bit-identical, since moimix uses a deterministic estimator;
divergence here would suggest a non-deterministic code path
(none expected with the current moimix version).

---

## Stage 2 — MOI / Fws (microhap seam stub)

Source: `workflow/rules/moi_microhap.smk`

### `microhap_moi_stub`

STUB. Emits the seam path so downstream rules wire up identically;
actually running it fails fast with the contract message.

**WHAT.** echo MicrohapNotImplemented >&2 ; exit 1

**WHY.** Refuses to let a half-built microhap path silently produce a
biallelic-shaped Fws file from a non-biallelic dataset.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/moi/fws_MOI.tsv  (never actually written)`

**TRY.** switch cohort.input_type to "wgs" — this rule disappears, the
real moimix chain wires in instead.

---

## Stage 3 — Structure (analysis core)

Source: `workflow/rules/03_structure.smk`

### `admixture_run`

Run ADMIXTURE for a single K with --cv cross-validation. Operates on
the LD-pruned bfile from the prep seam. K values run in parallel via
Snakemake's wildcard expansion.

**WHAT.** admixture --cv N cleaned.bed K  (staged in a per-stage dir so .Q
and .P land under outputs/structure/admixture/)

**WHY.** ADMIXTURE is the slide-7 ancestry-bar method. CV error vs K is
the standard model-selection diagnostic.

**TUNABLES.**

- `structure.admixture_cv_folds` = `5`

**OUTPUT.** `{outputs}/structure/admixture/cleaned.{K}.{Q,P}
+ {logs}/structure/admixture_K{K}.log`

**TRY.** bump admixture_k_max to 12 to scan further if the CV-vs-K curve
hasn't visibly bottomed out by K=10 (rare; usually flatlines).

---

### `admixture_cv_table`

Collect CV error from every K's log into one tidy TSV. The pattern in
ADMIXTURE's log is `CV error (K=N): <err>` — awk parses it out.

**WHAT.** awk extract over the admixture_K*.log files

**WHY.** The CV curve is the model-selection input + the slide-7 panel.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/structure/admixture/cv_error.tsv  (K, CV_error)`

**TRY.** grep `CV error` admixture_K*.log to spot-check parsing.

---

### `select_best_k`

Pick the best K — minimum CV error from cv_error.tsv unless overridden
by structure.admixture_k. Writes a one-line file consumed by every
downstream rule that needs to know which K to read.

**WHAT.** awk min over CV_error column, OR echo of the override.

**WHY.** Per DESIGN §7 decision 1: best-K is the CV minimum, with a manual
override knob. Centralising the choice in one file keeps every
downstream rule consistent.

**TUNABLES.**

- `structure.admixture_k`: *(not set in this config)*

**OUTPUT.** `{outputs}/structure/best_k.txt`

**TRY.** set structure.admixture_k to an integer one step away from the
auto pick — re-run from assign_clusters down and inspect whether
the extra component tracks a real geographic / host split.

---

### `admixture_cv_plot`

CV-error-vs-K diagnostic plot (the slide-7 small panel).

**WHAT.** scripts/R/plot_admixture_cv.R

**WHY.** Lets a user read off the best K visually and notice if the curve
is flat (multiple K's nearly equivalent — adds a wrinkle to the
interpretation that the single min-K choice doesn't show).

**TUNABLES.** *(none)*

**OUTPUT.** `{reports}/figures/admixture_cv.png (+ .svg)`

**TRY.** eyeball the plot — if K=2 and K=3 are within 0.005 of each
other, try running assign_clusters at K=2 manually and compare.

---

### `assign_clusters`

Read the best-K .Q, assign each sample to its dominant ancestry
component, label components per structure.cluster_labelling.

Declared as a `checkpoint` (not a plain rule) because the labels it
emits drive Stage 4's per-cluster IBD wildcard expansion — the set of
cluster names isn't known until this file exists. Snakemake re-
evaluates the DAG once the checkpoint's output is on disk.

**WHAT.** scripts/R/assign_clusters.R — label-mode-driven (numbered |
auto | reference). Reads `best_k.txt` so it never hardcodes a K.

**WHY.** The core lift of the agnostic refactor (DESIGN §5 item 2). V1
hardcoded K=3 and literal `Mn/Mf/Peninsular` labels — neither
generalises. The agnostic version:
numbered   — components → Cluster_1..K (zero-metadata default).
auto       — components labelled by majority geography role.
reference  — components labelled by majority group role
(this is V1's behaviour, applied via the role map).

**TUNABLES.**

- `structure.cluster_labelling` = `'reference'`
- `structure.admixture_k`: *(not set in this config)*

**OUTPUT.** `{outputs}/structure/admix_clusters.tsv
(cols: sample_id, Component, Proportion, Cluster)`

**TRY.** flip cluster_labelling between numbered/auto/reference on the
same .Q — the component memberships are identical, only the
human-readable Cluster column changes. Confirms the lift is
label-only, never re-runs ADMIXTURE.

---

### `pca`

Compute PCA from the cleaned (non-LD-pruned) bfile.

**WHAT.** plink2 --pca on cleaned.{bed,bim,fam}

**WHY.** PC1/PC2 separation is the slide-9 panel. Run on the cleaned set
rather than the LD-pruned one to match V1.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/structure/Pk.eigenvec, Pk.eigenval`

**TRY.** compare PC1 variance against the by-region split — if PC1 is
country-only and PC2 hits region, your structure is hierarchical.

---

### `pca_variance`

Convert .eigenval (one variance per line) into a tidy two-column TSV
(PC, variance_percent) for the Quarto report and Stage-3b axis labels.

**WHAT.** awk sum-and-normalise

**WHY.** Eigenvalues alone are unitless; percent-variance is what plots
and tables display.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/structure/pca_variance.tsv`

**TRY.** sum the first three rows — typically explains the bulk of
between-region structure on this scale of cohort.

---

### `distance_matrix`

Pairwise IBS-distance matrix used for the NJ tree in Stage 3b.

**WHAT.** plink (1.9) --distance square

**WHY.** plink2 alpha dropped --distance; the 1.9 implementation is the
reference. Square output is the matrix shape ape::nj() expects.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/structure/Pk.dist + Pk.dist.id`

**TRY.** convert to a phylo and inspect a quick `ape::nj` plot to sanity-
check the cohort splits before Stage-3b adds colour + labels.

---

## Stage 3 — Structure prep (WGS path)

Source: `workflow/rules/structure_prep_wgs.smk`

### `exclude_high_moi`

Drop high-MOI samples (from Stage 2) from the Stage-1 seam VCF.

**WHAT.** bcftools view -S ^exclude --force-samples

**WHY.** Polyclonal samples break single-genotype assumptions in
ADMIXTURE / PCA / NJ. Filter once, here, at the entry of Stage 3.

**TUNABLES.**

- `moi.fws_exclusion_cutoff` = `0.95`

**OUTPUT.** `{outputs}/structure/snps.moi_excluded.vcf.gz`

**TRY.** after raising fws_exclusion_cutoff, count the remaining samples
and watch Stage-3's final cohort drop.

---

### `normalise_vcf`

Split multiallelic records, left-align, set ID = CHROM:POS:REF:ALT.

**WHAT.** bcftools norm -m-any | bcftools norm --check-ref w -f <ref>
| bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT'

**WHY.** PLINK needs unique variant IDs and biallelic records; ADMIXTURE
is strict about both. Doing it here, once, keeps every downstream
tool happy without re-normalisation.

**TUNABLES.**

- `reference.fasta` = `'data/reference/strain_A1_H.1.Icor.fasta'`

**OUTPUT.** `{outputs}/structure/snps.normalised.vcf.gz`

**TRY.** head the .vcf.gz — IDs should look like
`ordered_PKNH_01_v2:12345:A:C` (contig:pos:ref:alt).

---

### `build_variant_chrom_map`

Build the variant-ID → integer-chrom-code map for PLINK by joining the
normalised VCF's CHROM column to {outputs}/setup/contig_map.tsv.

**WHAT.** bcftools query -f '%ID\t%CHROM\n' | awk join against contig_map

**WHY.** PLINK / ADMIXTURE want integer chromosome codes. Stage 0 already
derived (contig_name → integer) from the reference .fai (in fai
order); we reuse that single source of truth instead of V1's
PKA1H1-specific `gsub("ordered_PKNH_",""); gsub("_v2","")`.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/structure/chrom_update.txt   (variant_ID  int_code)`

**TRY.** diff this against V1's outputs/structure/chrom_update.txt (in
agnostic dev) — the integer codes should match the .fai order.

---

### `vcf_to_plink`

Convert the normalised VCF to PLINK (.bed/.bim/.fam) with integer
chromosome codes and sample IDs preserved verbatim.

**WHAT.** plink2 --vcf … --double-id --update-chr <map> --sort-vars
--make-pgen → plink2 --pfile … --make-bed
(two-step because plink2's --make-bed with --update-chr +
--sort-vars requires a sorted .pgen intermediate)

**WHY.** ADMIXTURE / downstream PLINK rules want a clean integer-coded
bfile with unique sample IDs.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/structure/Pk.{bed,bim,fam}`

**TRY.** head Pk.bim — chrom column should be 1..N integers, no contig
strings; sample IDs in Pk.fam should match VCF sample IDs verbatim.

---

### `plink_freq_missing`

Initial --freq + --missing on the full PLINK set — drives duplicate
detection (lowest-missingness replicate) + the variant-filter chain.

**WHAT.** plink2 --bfile Pk --freq --missing

**WHY.** Per-sample F_MISS is what dedup needs; allele frequencies feed
the final --maf filter further down.

**TUNABLES.** *(none)*

**OUTPUT.** `Pk.afreq, Pk.smiss, Pk.vmiss`

**TRY.** eyeball Pk.smiss F_MISS distribution; any extreme outliers
usually flag a sequencing failure.

---

### `find_duplicates`

Identify replicate samples (same biological isolate sequenced twice)
and write the IDs of the higher-missingness replicates to a PLINK
--remove file.

**WHAT.** scripts/R/find_duplicates.R reads .fam + .smiss, strips the
configured duplicate_id_pattern (default: don't strip), then
normalises out underscores/hyphens to derive a dedup base ID.
Keeps the lowest-missing replicate per base ID.

**WHY.** V1 stripped Illumina lane suffix `_DK.*` which is PKA1H1-cohort-
specific. The agnostic version takes the regex from config; null
= no strip (still dedups on the underscore/hyphen normalisation).

**TUNABLES.**

- `structure.duplicate_id_pattern` = `'_DK.*'`

**OUTPUT.** `{outputs}/structure/Pk.dups`

**TRY.** set duplicate_id_pattern to ".*" — every sample becomes a dup
of every other; rule will keep one. Useful sanity check that
the strip is happening at all.

---

### `final_filters`

Sequential 4-step filter chain (legacy order V→S→V→M):
1. Remove duplicate replicates + lenient variant filter (--geno 0.20).
2. Sample-missingness filter (--mind) on the cleaned variants.
3. Stricter variant filter (--geno) on the post-sample set.
4. MAF filter (--maf).

**WHAT.** four plink2 passes, intermediates cleaned at the end.

**WHY.** V1 found that a strict sample-first order dropped 73% of the
Indonesia cohort because the 1.4M-variant unfiltered set carries
many low-coverage sites that drag down per-sample missingness.
The V→S→V→M order matches HPC behaviour and recovers the cohort.

**TUNABLES.**

- `structure.max_sample_missing` = `0.1`
- `structure.max_variant_missing` = `0.1`
- `structure.min_maf` = `0.01`

**OUTPUT.** `{outputs}/structure/cleaned.{bed,bim,fam}`

**TRY.** bump max_sample_missing to 0.05 and re-run — every step shows
the cohort shrinking by tens.

---

### `ld_prune`

LD-prune the cleaned bfile before ADMIXTURE.

**WHAT.** plink2 --indep-pairwise 50 5 0.5 → plink --extract --make-bed
(plink 1.9 for the extract step; plink2 alpha segfaults on this
combination, both produce identical bfiles).

**WHY.** ADMIXTURE docs recommend unlinked SNPs. Without pruning, K=5
alone took >3 h on a single core in V1; pruning takes it to
~minutes total. The standard 50/5/0.5 window typically retains
5-10% of variants with full population-structure signal.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/structure/cleaned.ld.{bed,bim,fam}
+ {outputs}/structure/cleaned.prune.in (kept-variant list)`

**TRY.** after pruning, check `wc -l cleaned.prune.in` vs cleaned.bim —
~5-10% kept is typical for an outbred eukaryote at this scale.

---

## Stage 3 — Structure prep (microhap seam stub)

Source: `workflow/rules/structure_prep_microhap.smk`

### `microhap_structure_prep_stub`

STUB. Emits the seam paths so the downstream DAG wires up identically;
actually running it fails fast with the contract message.

**WHAT.** echo MicrohapNotImplemented >&2 ; exit 1

**WHY.** Refuses to let a half-built microhap path produce a biallelic-
shaped bfile from a non-biallelic dataset.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/structure/cleaned.{bed,bim,fam} + cleaned.ld.{bed,bim,fam}
(never actually written)`

**TRY.** swap cohort.input_type to "wgs" — this rule disappears and the
real WGS prep wires in instead.

---

## Stage 3b — Structure figures

Source: `workflow/rules/03_figures.smk`

### `plot_admixture_bars`

Single-panel stacked ADMIXTURE bar plot — samples ordered by dominant
cluster then dominant proportion, coloured by component.

**WHAT.** scripts/R/plot_admixture_bars.R

**WHY.** The slide-7 cohort-wide ancestry summary. Always renders; no
role gating because dominant-cluster sort uses admix_clusters.tsv.

**TUNABLES.** *(none)*

**OUTPUT.** `{reports}/figures/admixture_bars.png/.svg`

**TRY.** if the bar plot looks "noisy" you may want a different K — pass
it via `structure.admixture_k` and re-run.

---

### `plot_pca`

PC1-vs-PC2 scatter coloured by ADMIXTURE cluster.

**WHAT.** scripts/R/plot_pca.R

**WHY.** Always renders — uses admix_clusters.tsv's Cluster column which
is always produced (in numbered mode at worst).

**TUNABLES.** *(none)*

**OUTPUT.** `{reports}/figures/pca.png/.svg`

**TRY.** swap the colour role to geography in plot_pca_by_geography
(below) to see if PC1/PC2 maps onto region instead of cluster.

---

### `plot_nj_tree`

Unrooted Neighbour-Joining tree from PLINK's IBS distance matrix,
coloured by ADMIXTURE cluster.

**WHAT.** scripts/R/plot_nj_tree.R (ape::nj + ggtree daylight layout)

**WHY.** Visualises clade structure without assuming an outgroup.
Coloured edges show within-cluster vs between-cluster edges.

**TUNABLES.** *(none)*

**OUTPUT.** `{reports}/figures/njt.png/.svg`

**TRY.** open the svg and inspect: tight cluster cliques + sparse
connections between clusters typically mean clean structure.

---

## Stage 4 — IBD (per-cluster hmmIBD + clonal detection)

Source: `workflow/rules/04_ibd.smk`

### `cluster_membership`

Write the cluster's keep-list (FID + IID per line) for PLINK --keep.

**WHAT.** read admix_clusters.tsv, filter to this cluster, emit sample IDs.

**WHY.** Replaces V1's cluster_exclude_lists rule + its 4-source union
(lab isolates + manually-curated dups). Those are gone from
clean input; membership from admix_clusters is the only
agnostic input to per-cluster subsetting.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/ibd/{cluster}/keep.txt`

**TRY.** peek at the file — one line per sample, `sample sample`
format matching plink --keep expectations.

---

### `cluster_maf_filter`

Per-cluster bfile + VCF: keep only cluster members, apply within-
cluster --maf.

**WHAT.** plink2 --keep + --maf → make-bed + recode vcf.

**WHY.** Cluster-specific MAF is the point of per-cluster IBD; variants
polymorphic in one cluster may be fixed in another.

**TUNABLES.**

- `ibd.cluster_min_maf` = `0.01`

**OUTPUT.** `{outputs}/ibd/{cluster}/cleaned.{bed,bim,fam,vcf.gz}`

**TRY.** bump ibd.cluster_min_maf to 0.05 for a tighter within-cluster
variant set; usually shrinks hmmIBD input by ~30-40%.

---

### `genotype_table`

Per-cluster VCF → hmmIBD genotype-table format.

**WHAT.** scripts/R/genotype_table.R, using outputs/setup/contig_map.tsv
to map CHROM strings → integers.

**WHY.** hmmIBD wants integer CHROM. The contig_map is derived once
in Stage 0 from the reference .fai — the agnostic replacement
for V1's hardcoded `ordered_PKNH_/_v2` gsub.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/ibd/{cluster}/hmmIBD_input.tsv`

**TRY.** `head -3` the tsv — chrom column should be small ints in
the same order as contig_map.tsv.

---

### `run_hmmibd`

Run hmmIBD pairwise within a cluster.

**WHAT.** hmmIBD -i <genotype_table.tsv> -o <prefix>

**WHY.** Produces .hmm_fract.txt (per-pair fract_sites_IBD) which is the
basis for both connectivity plots and clonal-pair detection.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/ibd/{cluster}/{cluster}.hmm_fract.txt + .hmm.txt`

**TRY.** inspect the hmm_fract.txt — the fract_sites_IBD histogram
should have a bulk near 0 and a small clonal spike near 1.

---

### `combined_vcf`

Export the Stage-3 `cleaned` bfile straight to VCF. No exclusion step:
the cleaned bfile already has dups + lab isolates + manually-curated
samples removed (see data/INPUT_PREP.md for the Indo cohort's prep).

**WHAT.** plink2 --recode vcf bgz on cleaned.{bed,bim,fam}

**WHY.** Stage 5 (introgression) reads the combined genotype table for
all Stage-3-passing samples.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/ibd/combined/cleaned.vcf.gz`

**TRY.** n = wc -l cleaned.fam should equal the sample count in the
combined VCF header.

---

### `combined_genotype_table`

All-samples hmmIBD-format table (Stage 5 introgression input).

**WHAT.** scripts/R/genotype_table.R on the combined VCF using contig_map.

**WHY.** Introgression needs every sample under the same genotype-encoded
matrix, not per-cluster subsets.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/ibd/combined/hmmIBD_input.tsv`

**TRY.** Stage 5's rules read from this path; touch it to force a
downstream re-run without redoing the plink export.

---

### `clonal_clusters`

Aggregate every cluster's hmm_fract.txt, flag pairs at or above
ibd.clonal_ibd_threshold, derive clonal groups as connected components
of the pair graph. Join to canonical geography/date roles for a tidy
per-sample table.

**WHAT.** scripts/R/clonal_clusters.R over every cluster's fract file.

**WHY.** V1 hardcoded a Peninsular-only pipeline plus a manual 8-sample
tribble for downstream. The agnostic version computes clonal
groups from the IBD graph itself; the cohort-specific tribble
never returns.

**TUNABLES.**

- `ibd.clonal_ibd_threshold` = `0.95`
- `ibd.focal_cluster` = `'Peninsular'`

**OUTPUT.** `{outputs}/ibd/clonal_clusters.tsv       (always)
{outputs}/ibd/focal_<name>_clones.tsv    (only if
ibd.focal_cluster set)`

**TRY.** drop ibd.clonal_ibd_threshold to 0.5 to see near-clonal
structure; the connected-component count will spike.

---

## Stage 6 — Selection (rehh iHS)

Source: `workflow/rules/06_selection.smk`

### `run_rehh_ihs`

Per-model whole-genome iHS scan via rehh. Single-threaded — rehh's
`data2haplohh` / `scan_hh` do not parallelise usefully on this scale.

**WHAT.** scripts/R/run_rehh_ihs.R — subsets the combined Stage-4 VCF to
case + control per the model's role-expression filters, splits
by contig, accumulates a whole-genome scan, runs
`calc_candidate_regions`. Emits empty (header-only) TSVs +
a one-line summary on a negative result.

**WHY.** iHS detects loci under recent positive selection. Splitting by
case/control lets the same machinery run any comparison the
user cares about — no cohort hardcoding beyond the filter
expression the user writes in config.

**TUNABLES.**

- `selection.ihs_threshold` = `4`
- `selection.ihs_p_threshold` = `0.0001`
- `selection.window_size_bp` = `10000`
- `selection.window_overlap_bp` = `1000`
- `selection.min_extr_markers` = `3`

**OUTPUT.** `{outputs}/selection/{model}/candidate_regions_iHS.tsv,
{outputs}/selection/{model}/ihs_table.tsv,
{outputs}/selection/{model}/ihs_summary.tsv`

**TRY.** add a second model to selection.models (e.g. a
North-Kalimantan-vs-Sabah scan) and re-run — the rule
expands over the added name automatically.

---

### `plot_ihs_scan`

Genome-wide iHS scan + -log10(p) plot for a model.

**WHAT.** scripts/R/plot_ihs_scan.R

**WHY.** The scan is the standard visualisation; alternating chromosome
fills + 0.1% / 99.9% quantile dotted lines flag outliers.

**TUNABLES.** *(none)*

**OUTPUT.** `{outputs}/selection/{model}/ihs_scan.{png,svg} +
{outputs}/selection/{model}/ihs_pvalue.{png,svg}`

**TRY.** inspect the p-value plot for peaks above the 0.001 quantile
line; the candidate_regions_iHS.tsv should hit those peaks.

---

## workflow/rules/04_figures.smk

---
