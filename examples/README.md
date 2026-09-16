# Pull and play — the public example

A complete, runnable cohort that ships with the repository. Clone ZOOMAL-Flow,
fetch one reference genome from PlasmoDB, and run every stage of the pipeline
on real *Plasmodium knowlesi* data without needing access to anything
embargoed.

**71 publicly-archived ENA isolates from Malaysia · 18,279 SNPs · ~4 MB.**

It exists for two reasons: to give a new student something that works on the
first try, and to let anyone check that their environment is correct before
pointing the pipeline at data that matters.

---

## Run it

From the root of your clone (the repository root is the workspace):

```bash
pixi install && pixi run setup          # if you have not built the env yet
```

**1. Get the reference** (~24 MB; not shipped — see "Why no reference?" below).

Go to **<https://plasmodb.org>**. PlasmoDB requires a **free account, and you
must be logged in to download** — register once, then:

> **PlasmoDB** → search for *Plasmodium knowlesi* strain **A1-H.1** →
> **Download** → **Genome** → the genomic **FASTA** file,
> `PlasmoDB-<release>_PknowlesiA1H1_Genome.fasta`

This example is built against **PlasmoDB release 67**, which is the release the
example VCF's and mask's contig names were reconciled to. Take release 67 if it
is offered. A later release is also fine: the A1-H.1 assembly itself has not
changed — the FASTA headers still record `version=2017-09-19` — so the contigs
and coordinates are the same. Step 2 tells you either way, so you do not have to
take that on trust.

Save it here:

```
examples/data/reference/PlasmoDB-67_PknowlesiA1H1_Genome.fasta
```

(or put it anywhere and point `reference.fasta` in `examples/config.yaml` at it).

No direct download URL is given on purpose. PlasmoDB has reorganised its
download paths before — the ones this repository previously documented now
return 404 — and downloads sit behind a login anyway, so a copy-pasteable link
would rot and could not be fetched non-interactively regardless. The site plus a
verification step is the arrangement that keeps working.

**2. Check it and index it** — one command, a few seconds, and it saves you
finding out from a confusing error several minutes into the run:

```bash
pixi run bash examples/check_reference.sh
```

It builds the `.fai` if needed, then verifies every contig name and length
against the assembly this example was built for. If your download has different
FASTA headers, it tells you how to fix them.

**3. Run:**

```bash
pixi run snakemake --configfile examples/config.yaml --cores 4
```

Takes a few minutes on a laptop. Results land in `examples/outputs/`, figures
in `examples/reports/figures/`, per-rule logs in `examples/logs/` — all
gitignored.

### What you should see

```
outputs/structure/full/best_k.txt            → 3
outputs/structure/full/admix_clusters.tsv    → Cluster_1 24 / Cluster_2 14 / Cluster_3 33
reports/figures/full/pca.png                 → three tight, well-separated clouds
```

ADMIXTURE picks K = 3 by cross-validation error, and the three components it
finds correspond exactly — all 71 samples, no exceptions — to the three
clusters these isolates fall into in the full 981-sample cohort this example
was drawn from. That is the acceptance test: if you get K = 3 and a clean
three-way PCA, your install is working.

Other stages run too and write real output: Fws/MOI, per-cluster hmmIBD,
clonal-group detection, introgression windows for all three cluster pairs, and
the full/de-clonalized comparison. The iHS selection scan is the one stage that
produces nothing, deliberately — see "What is switched off" below.

---

## The data

Every sample is a public ENA accession. Nothing here is under embargo, and the
metadata carries only published fields: the accession, country, region, and
host.

| Source-cohort group | Region | n | Label with the shipped seed | Accessions |
|---|---|---|---|---|
| Peninsular | Peninsular | 24 | `Cluster_1` | ERR3374031–3374058 (24 of the 28 in that range) |
| Mn | Sarawak | 14 | `Cluster_2` | ERR2214838, ERR2214849, ERR2214852, ERR274224–225, ERR366425, ERR985411–416, ERR985418–419 |
| Mf | Sarawak | 33 | `Cluster_3` | ERR274221–222, ERR366426, ERR985372–374, ERR985377–385, ERR985387–394, ERR985398–404, ERR985406–408 |

The cluster *numbers* carry no meaning — they are an artefact of the seed, not
a property of the data. ADMIXTURE labels its components in whatever order it
converges on, so which group is `Cluster_1` depends on
`structure.admixture_seed`. With the shipped seed you will reproduce the
right-hand column exactly; change the seed and the same three groups can come
back under permuted names. Read the partition, not the labels: the stable
handles are the group sizes (24 / 14 / 33) and the accessions in each.

`ACCESSIONS.tsv` in this directory is the authoritative per-sample list, with
each sample's country, region, host, and the cluster it occupies in the source
cohort.

### Provenance and citation

This example is a **public subset derived from the dataset published in**:

> Westaway JAF, Diez Benavente E, Auburn S, Kucharski M, Aranciaga N, Nayak S,
> William T, Rajahram GS, Piera KA, Braima K, Tan AF, Alaza DA, Barber BE,
> Drakeley C, Amato R, Sutanto E, Trimarsanto H, Jelip J, Anstey NM, Bozdech Z,
> Field M, Grigg MJ. **Genomic epidemiology of *Plasmodium knowlesi* reveals
> putative genetic drivers of adaptation in Malaysia.** *PLOS Neglected
> Tropical Diseases* 19(3): e0012885 (2025).
> doi:[10.1371/journal.pntd.0012885](https://doi.org/10.1371/journal.pntd.0012885)

**If you use this example in anything you publish or present, cite that paper.**

That study combined newly-sequenced Malaysian genomes with publicly archived
ones, and references the original public ENA studies the archived samples came
from. The 71 samples here are exclusively the **publicly archived** portion —
every one is an ENA accession that was already public — so citing the paper
above covers attribution, and its reference list is where to follow the
accessions back to their original studies.

Individual accessions also resolve directly at
<https://www.ebi.ac.uk/ena/browser/text-search>, which gives the study
accession and submitting centre for any one of them.

Geography does **not** separate these three groups. One of them is from
Peninsular Malaysia, but the other two are *both* from Sarawak — they are
host-associated lineages (the `Mf` and `Mn` groups of the source cohort), not
geographic ones. That is why the config sets `cluster_labelling: "numbered"`:
`auto` names each component after its members' majority region, which here
would hand the same name, "Sarawak", to two different clusters. Numbered labels
are the honest option for a cohort whose structure is not geographic — at the
cost, noted above, of the numbers themselves meaning nothing.

### What was done to the VCF

Built by `build_example.sh` from the full cohort VCF, deterministically:

1. **Samples** — every public accession in the cohort (`^(ERR|ERS|ERX|SRR|SRS)[0-9]+$`),
   which is all 71 of them. Not a balanced subset: PLINK2 refuses to estimate
   linkage disequilibrium from fewer than 50 samples, and Stage 3 LD-prunes, so
   a tidy 8-per-cluster example builds fine and then dies at `ld_prune`.
2. **Sites** — biallelic PASS SNPs on the 14 nuclear chromosomes, then
   re-filtered against these 71 samples (called in ≥ 90% of them, MAC ≥ 2).
   613,518 SNPs survive that.
3. **Thinned** to ≥ 1100 bp spacing → **18,279 SNPs**, which is what keeps the
   committed file near 4 MB. Even after thinning there are ~9 SNPs in each
   10 kb introgression window, comfortably above the configured minimum of 5.
4. **FORMAT trimmed** to `GT`, `AD`, `DP`. `AD` stays because moimix needs
   within-sample allele balance to estimate Fws; `PL`/`SB`/`GQ` and friends are
   bulk with nothing downstream reading them.
5. **Contigs renamed** to the PlasmoDB assembly's names, and the header contig
   lengths rewritten from its `.fai` — see below.

Numbers through the run: 18,279 sites in → 17,842 after QC masking and
filtering → 15,336 after LD pruning → ADMIXTURE.

---

## The reference must match

This is the one place it is easy to go wrong, so it is worth a paragraph.

The cohort VCF was called against a local copy of the *P. knowlesi* A1-H.1
assembly whose contigs are named `ordered_PKNH_01_v2` … `ordered_PKNH_14_v2`.
PlasmoDB distributes the same assembly with its contigs named by their ENA
accession, `LT727648` … `LT727663`. Different names, same sequence.

That was verified rather than assumed, per contig, at 10/25/50/75/95% of its
length: the two are identical base-for-base at identical coordinates. The local
copy is slightly *longer* on most contigs (+99 bp, and +1208 on chr4, +7138 on
chr11), and every one of those extra bases is an `N` pad on the end. So
reconciling the two is a pure rename, with no coordinate shift to worry about —
and `build_example.sh` still asserts that no site in the example falls beyond a
PlasmoDB contig's end before it will write the file.

The example VCF is therefore shipped **already renamed to PlasmoDB's contig
names, with PlasmoDB's contig lengths in its header**. It is built to run
against the reference you download, not against the one on the maintainer's
machine.

To check you fetched the right assembly, compare your `.fai` against this:

| Contig | Length | | Contig | Length |
|---|---|---|---|---|
| LT727648 | 892,471   | | LT727656 | 2,195,390 |
| LT727649 | 787,547   | | LT727657 | 1,516,883 |
| LT727650 | 1,050,471 | | LT727658 | 2,361,529 |
| LT727651 | 1,130,816 | | LT727659 | 3,161,709 |
| LT727652 | 782,346   | | LT727660 | 2,596,912 |
| LT727653 | 1,080,557 | | LT727661 | 3,267,653 |
| LT727654 | 1,522,076 | | LT727662 | 5,957 (mitochondrion) |
| LT727655 | 1,924,455 | | LT727663 | 34,438 (apicoplast) |

`examples/check_reference.sh` checks exactly this table for you, which is why
step 2 above is worth the few seconds. Failing that, a mismatched reference
shows up as a contig-name error in Stage 0 — loudly, not as quietly wrong
results.

Note that the organelles are excluded by **accession**
(`exclude_contigs: ["LT727662", "LT727663"]`), not by the `["MIT", "API"]`
patterns a cohort using the other naming would use. Exclusion patterns are
regexes matched against contig names, and there is no "MIT" substring in
`LT727662` to match.

### Why no reference?

It is ~24 MB of sequence that PlasmoDB already distributes and versions
properly, behind a login. Committing a copy would bloat the repository, a stale
copy would be worse than none, and redistributing it is not ours to do.
Pointing at the source is the correct arrangement.

---

## What is switched off, and why

The example runs every stage the data can support. Three things are
deliberately null, and each one demonstrates the pipeline's
graceful-degradation behaviour rather than a gap:

| Off | Config | Why |
|---|---|---|
| Temporal clonality plot | `roles.date: null` | Collection dates are not public for these samples. The stage logs a note and skips. |
| Province map | `metadata.gis: null` | No public per-sample coordinates. |
| Gene-family masking | `gene_family_filters: []`, `gff: null` | Needs the PlasmoDB GFF, a second download. **On a real *P. knowlesi* cohort, turn this on** — SICAvar and KIR are hypervariable and generate false-positive introgression calls. |
| iHS selection scan | `selection.models: []` | An iHS scan needs a case/control contrast, which is a cohort-specific hypothesis. The example has none to make. |
| Focal enrichment test | `focal_group: null` | Same reason: asking whether introgression is enriched in a named subgroup is a claim about a particular cohort. |

Some thresholds are also relaxed from what a real cohort would use, and the
reasons are commented inline in `config.yaml`: `max_sample_missing` and
`max_variant_missing` at 0.25 (the archival Peninsular batch is
lower-coverage, and 0.10 would drop it and take the third cluster with it), and
`fws_exclusion_cutoff` at 0.90.

**Read the introgression output as a demonstration, not a finding.** On a real
cohort the per-cluster support floor is *derived* by permutation against an
explicit FDR target, and is cohort-specific — inheriting a number from another
study is wrong. With 14–33 samples per cluster a derived floor is not
meaningful, so `per_cluster_min_pct: 0.25` is just a value that keeps the stage
running so you can see the shape of its outputs.

---

## Rebuilding it

`build_example.sh` regenerates every committed file here. You do not need it to
*run* the example — it exists so the example's provenance is code rather than
prose, and so the selection rule can be changed and re-applied.

It reads the full cohort VCF, so it only runs on a machine that has it. Its
tunables are `PER_CLUSTER` and `MIN_SPACING` at the top.

Everything it writes is public by construction, and it proves that before
finishing: it selects only from public-accession IDs, carries only public
metadata fields, and asserts that the sample list in the example VCF's own
header shares **zero** IDs with the non-public cohort. If that assertion fails
it deletes the VCF it just wrote and exits non-zero.

---

## A note on what this is

This is a demonstration cohort: small, heavily thinned, and chosen to exercise
the pipeline rather than to answer a question. It shows that the pipeline runs,
that your environment is correct, and that the structure it recovers is the
structure that is there. It is not a result about *P. knowlesi*, and no figure
from it should be presented as one.
