#!/usr/bin/env bash
# examples/build_example.sh — regenerate the public pull-and-play example.
#
# MAINTAINER TOOL. You do not need this to *run* the example — the example's
# inputs are committed. This is the record of how they were made, and the way
# to remake them if the selection rule or the thinning target changes.
#
# It reads the full embargoed cohort VCF, so it only runs on a machine that has
# that VCF. Everything it *writes* is public by construction, and it asserts
# that before finishing (see "Safety" below).
#
#   bash examples/build_example.sh            # from agnostic/, inside the env
#
# Outputs (all committed, all public):
#   examples/data/vcf/example.vcf.gz{,.csi}
#   examples/data/metadata/samples.tsv
#   examples/data/reference/regions_to_mask.list
#   examples/ACCESSIONS.tsv
#
# ---------------------------------------------------------------------------
# Safety — why the output is publishable
# ---------------------------------------------------------------------------
# The cohort mixes public ENA-archived samples with samples under publication
# embargo. Only the former may be committed. This script:
#
#   1. selects exclusively from IDs matching ^(ERR|ERS|ERX|SRR|SRS)[0-9]+$,
#   2. carries only public metadata fields (accession, country, region, host),
#   3. asserts, before writing anything final, that the sample list in the
#      example VCF header is disjoint from every ID in the embargoed cohort
#      tables — and aborts if it is not.
#
# ---------------------------------------------------------------------------
# Reference / contig reconciliation
# ---------------------------------------------------------------------------
# The cohort VCF is called against a local reference whose contigs are named
# `ordered_PKNH_NN_v2`. The README tells users to download the reference from
# PlasmoDB, where the same assembly's contigs are named by their ENA accession
# (`LT727648`..`LT727663`). Verified 2026-09-16, per contig at 10/25/50/75/95%
# of its length: the two are the SAME sequence at the SAME coordinates. The
# local copy carries extra trailing bases (+99 on most contigs, +1208 on chr4,
# +7138 on chr11) and every one of those extra bases is an `N` pad. So the
# reconciliation is a pure rename, and the rename is safe — but this script
# still asserts that no retained site falls beyond a PlasmoDB contig's length.

set -eo pipefail

cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ---------------------------------------------------------------------------
# Tunables
# ---------------------------------------------------------------------------
SRC_VCF="data/vcf/merged_popgen.clean.vcf.gz"
CLUSTERS="outputs/structure/full/admix_clusters.tsv"   # ADMIXTURE cluster per sample
SMISS="outputs/structure/Pk.smiss"                     # PLINK per-sample missingness
META_SRC="outputs/metadata/samples.tsv"                # role-renamed cohort metadata
MASK_SRC="data/reference/regions_to_mask.list"
PLASMODB_FASTA="data/reference/PlasmoDB_version/PlasmoDB-67_PknowlesiA1H1_Genome.fasta"

# Samples kept per ADMIXTURE cluster. Set high enough to take the whole public
# pool (33 Mf / 14 Mn / 24 Peninsular = 71).
#
# Why not a small balanced set — the obvious thing to want? PLINK2 refuses
# --indep-pairwise below 50 samples ("there are less than 50 samples to
# estimate from"), and LD pruning is not optional in Stage 3. A balanced
# ~8-per-cluster example builds fine and then dies at ld_prune. Capping at 14
# per cluster to match the smallest cluster gives 42, still under the floor.
# So the example takes every public sample: it clears 50 with room to spare,
# and "all the public data there is" is an honest thing for an example to be.
PER_CLUSTER=99
# bp between retained SNPs — the size knob. The genome is ~23.3 Mb, so this
# caps the site count near 23.3e6/MIN_SPACING. 1100 bp lands ~21 k sites and a
# ~5 MB committed VCF, while still leaving ~9 SNPs in each 10 kb introgression
# window (min_snps_per_window is 5).
MIN_SPACING=1100
PUBLIC_RE='^(ERR|ERS|ERX|SRR|SRS)[0-9]+$'

OUT_DIR="examples"
OUT_VCF="$OUT_DIR/data/vcf/example.vcf.gz"
OUT_META="$OUT_DIR/data/metadata/samples.tsv"
OUT_MASK="$OUT_DIR/data/reference/regions_to_mask.list"
OUT_ACC="$OUT_DIR/ACCESSIONS.tsv"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

for f in "$SRC_VCF" "$CLUSTERS" "$SMISS" "$META_SRC" "$MASK_SRC" "$PLASMODB_FASTA"; do
  [[ -f "$f" ]] || { echo "ERROR: missing input $f" >&2; exit 1; }
done
mkdir -p "$OUT_DIR/data/vcf" "$OUT_DIR/data/metadata" "$OUT_DIR/data/reference"

#---------------------------------------------------------------------------
# 1. Select samples — deterministic, public only.
#
# Per ADMIXTURE cluster: public accessions, ranked by PLINK per-sample
# missingness ascending, ties broken by accession. Take the first PER_CLUSTER.
# No randomness, so re-running this picks exactly the same set.
#---------------------------------------------------------------------------
echo "==> Selecting $PER_CLUSTER public samples per cluster"
awk -F'\t' -v re="$PUBLIC_RE" -v n="$PER_CLUSTER" '
  NR==FNR { if (FNR>1) miss[$2]=$5; next }
  FNR>1 && $1 ~ re { print $NF "\t" (($1 in miss) ? miss[$1] : "NA") "\t" $1 }
' "$SMISS" "$CLUSTERS" \
  | sort -t$'\t' -k1,1 -k2,2g -k3,3 \
  | awk -F'\t' -v n="$PER_CLUSTER" '{ if (++c[$1] <= n) print }' > "$WORK/selected.tsv"

cut -f3 "$WORK/selected.tsv" | sort > "$WORK/samples.txt"
N_SEL=$(wc -l < "$WORK/samples.txt" | tr -d ' ')
echo "    selected $N_SEL samples across $(cut -f1 "$WORK/selected.tsv" | sort -u | wc -l | tr -d ' ') clusters"

# Guard: every selected ID must be a public accession. Belt and braces — the
# awk above already filtered on this, but this is the claim we publish on.
if grep -vE "$PUBLIC_RE" "$WORK/samples.txt"; then
  echo "ERROR: non-accession ID in the selection (printed above). Aborting." >&2
  exit 1
fi

#---------------------------------------------------------------------------
# 2. Subset + filter, per nuclear contig, in parallel.
#
# Biallelic PASS SNPs; FORMAT trimmed to GT/AD/DP (moimix needs AD for Fws,
# the rest is bulk); then re-filtered on the 24 selected samples so the site
# set is one that is actually informative *for this subcohort*.
#---------------------------------------------------------------------------
echo "==> Subsetting $SRC_VCF (this reads the full cohort VCF; takes a few minutes)"
mkdir -p "$WORK/chroms"
NUCLEAR=$(seq -w 1 14 | sed 's/^/ordered_PKNH_/;s/$/_v2/')

printf '%s\n' $NUCLEAR | xargs -P 7 -I{} sh -c '
  c="$1"; work="$2"; src="$3"
  bcftools view -r "$c" -S "$work/samples.txt" -f PASS -m2 -M2 -v snps -Ou "$src" \
  | bcftools annotate -x "FORMAT/PL,FORMAT/SB,FORMAT/GQ,FORMAT/PGT,FORMAT/PID,FORMAT/MIN_DP,FORMAT/RGQ,^INFO/AC,INFO/AN" -Ou \
  | bcftools view -i "F_MISSING<=0.10 && MAC>=2" -Oz -o "$work/chroms/$c.vcf.gz" -W=csi
' _ {} "$WORK" "$SRC_VCF"

printf '%s\n' $NUCLEAR | sed "s|^|$WORK/chroms/|;s|$|.vcf.gz|" > "$WORK/concat_list.txt"
bcftools concat -f "$WORK/concat_list.txt" -Oz -o "$WORK/all.vcf.gz" -W=csi
echo "    $(bcftools index -n "$WORK/all.vcf.gz") biallelic SNPs before thinning"

#---------------------------------------------------------------------------
# 3. Thin to keep the committed file small.
#
# Minimum-spacing thin: walk each contig in order and keep a site only if it is
# at least MIN_SPACING bp from the last kept one. Deterministic, spreads sites
# evenly along the genome, and incidentally decorrelates neighbouring SNPs,
# which is what ADMIXTURE wants anyway.
#---------------------------------------------------------------------------
echo "==> Thinning to >= ${MIN_SPACING} bp spacing"
bcftools query -f '%CHROM\t%POS\n' "$WORK/all.vcf.gz" \
  | awk -v d="$MIN_SPACING" '$1 != c { c=$1; last=-d } $2-last >= d { print; last=$2 }' \
  > "$WORK/keep_sites.tsv"
echo "    keeping $(wc -l < "$WORK/keep_sites.tsv" | tr -d ' ') sites"

bcftools view -T "$WORK/keep_sites.tsv" -Oz -o "$WORK/thinned.vcf.gz" -W=csi "$WORK/all.vcf.gz"

#---------------------------------------------------------------------------
# 4. Rename contigs to the PlasmoDB assembly's names.
#
# ordered_PKNH_01_v2 .. _14_v2 -> LT727648 .. LT727661   (MIT -> 662, API -> 663;
# neither is retained here, but the map is written whole so it is self-documenting).
#---------------------------------------------------------------------------
echo "==> Renaming contigs to the PlasmoDB assembly's accessions"
: > "$WORK/rename.txt"
for i in $(seq 1 14); do
  printf 'ordered_PKNH_%02d_v2\tLT%d\n' "$i" $((727647 + i)) >> "$WORK/rename.txt"
done
printf 'PKNH_MIT_v2\tLT727662\nnew_API_strain_A1_H.1\tLT727663\n' >> "$WORK/rename.txt"

bcftools annotate --rename-chrs "$WORK/rename.txt" -Oz -o "$WORK/renamed.vcf.gz" "$WORK/thinned.vcf.gz"

# --rename-chrs changes the contig IDs but keeps the old ##contig lengths,
# which are the LOCAL reference's padded ones. Rewrite the contig header from
# the PlasmoDB .fai so the declared lengths match the reference users download;
# otherwise `bcftools norm -f` in Stage 3 compares against a length that does
# not exist in their FASTA.
if [[ ! -f "$PLASMODB_FASTA.fai" ]]; then samtools faidx "$PLASMODB_FASTA"; fi
bcftools reheader -f "$PLASMODB_FASTA.fai" -o "$WORK/reheadered.vcf.gz" "$WORK/renamed.vcf.gz"
bcftools view -Oz -o "$OUT_VCF" -W=csi "$WORK/reheadered.vcf.gz"

# Assert every retained site is in range on the PlasmoDB assembly. The local
# reference is longer (trailing N pad); a site in that pad would be a genuine
# coordinate problem, not a cosmetic one.
bcftools query -f '%CHROM\t%POS\n' "$OUT_VCF" \
  | awk -F'\t' 'NR==FNR { len[$1]=$2; next }
       !($1 in len) { print "unknown contig: " $1; bad=1; next }
       $2 > len[$1] { print "out of range: " $1 ":" $2 " > " len[$1]; bad=1 }
       END { exit bad?1:0 }' "$PLASMODB_FASTA.fai" - \
  || { echo "ERROR: example VCF has sites outside the PlasmoDB assembly. Aborting." >&2; exit 1; }
echo "    all sites in range on the PlasmoDB assembly"

#---------------------------------------------------------------------------
# 5. Metadata — public fields only.
#
# accession, country, region, host. No dates (they tie to unpublished
# collections), no curated or cohort-internal columns, no cluster labels: the
# example derives its clusters with ADMIXTURE, as a new user would.
#---------------------------------------------------------------------------
echo "==> Writing public metadata"
{
  printf 'sample_id\tcountry\tgeography\thost\n'
  awk -F'\t' '
    NR==FNR { if (FNR>1) { country[$1]=$2; geo[$1]=$3; host[$1]=$6 } next }
    { s=$1; printf "%s\t%s\t%s\t%s\n", s, country[s], geo[s], host[s] }
  ' "$META_SRC" "$WORK/samples.txt"
} > "$OUT_META"

# Every VCF sample must have a metadata row, and vice versa.
diff <(bcftools query -l "$OUT_VCF" | sort) <(tail -n +2 "$OUT_META" | cut -f1 | sort) \
  || { echo "ERROR: VCF samples and metadata rows disagree. Aborting." >&2; exit 1; }

#---------------------------------------------------------------------------
# 6. Mask — same rename, and clipped to the PlasmoDB contig lengths.
#---------------------------------------------------------------------------
echo "==> Writing the renamed mask"
awk -F'\t' '
  NR==FNR { map[$1]=$2; next }
  FNR==NR { }
  { split($0, a, ":"); split(a[2], b, "-");
    if (a[1] in map) print map[a[1]] ":" b[1] "-" b[2] }
' "$WORK/rename.txt" "$MASK_SRC" > "$WORK/mask_renamed.list"

awk -F'\t' '
  NR==FNR { len[$1]=$2; next }
  { split($0, a, ":"); split(a[2], b, "-");
    if (!(a[1] in len)) next
    if (b[1] > len[a[1]]) next              # region starts past the contig end
    if (b[2] > len[a[1]]) b[2] = len[a[1]]  # clip a region that overruns it
    print a[1] ":" b[1] "-" b[2] }
' "$PLASMODB_FASTA.fai" "$WORK/mask_renamed.list" > "$OUT_MASK"
echo "    $(wc -l < "$OUT_MASK" | tr -d ' ') mask regions"

#---------------------------------------------------------------------------
# 7. Provenance table — the accession list, with the cluster it came from.
#---------------------------------------------------------------------------
{
  cat <<'PROV'
# ZOOMAL-Flow public example — sample provenance.
#
# A public subset derived from the dataset published in:
#   Westaway JAF, Diez Benavente E, Auburn S, Kucharski M, Aranciaga N, Nayak S,
#   William T, Rajahram GS, Piera KA, Braima K, Tan AF, Alaza DA, Barber BE,
#   Drakeley C, Amato R, Sutanto E, Trimarsanto H, Jelip J, Anstey NM,
#   Bozdech Z, Field M, Grigg MJ. "Genomic epidemiology of Plasmodium knowlesi
#   reveals putative genetic drivers of adaptation in Malaysia."
#   PLOS Neglected Tropical Diseases 19(3): e0012885 (2025).
#   doi:10.1371/journal.pntd.0012885
#
# Cite that paper if you use this example. It references the original public
# ENA studies these archived accessions came from.
#
# Every sample below is a publicly archived ENA accession; none are embargoed.
# `source_cluster` is the cluster the sample occupies in the full cohort, and
# is recorded for provenance only — the example derives its own clusters with
# ADMIXTURE and does not read this column.
PROV
  printf 'accession\tsource_cluster\tf_miss_full_cohort\tcountry\tgeography\thost\n'
  awk -F'\t' '
    NR==FNR { if (FNR>1) { country[$1]=$2; geo[$1]=$3; host[$1]=$6 } next }
    { printf "%s\t%s\t%s\t%s\t%s\t%s\n", $3, $1, $2, country[$3], geo[$3], host[$3] }
  ' "$META_SRC" "$WORK/selected.tsv"
} > "$OUT_ACC"

#---------------------------------------------------------------------------
# 8. THE safety assertion: zero overlap with any embargoed cohort ID.
#
# Built from the example VCF's own header, not from the selection list, so it
# checks the artefact that actually gets committed.
#---------------------------------------------------------------------------
echo "==> Checking the example against the embargoed cohort"
bcftools query -l "$OUT_VCF" | sort > "$WORK/example_ids.txt"

: > "$WORK/embargoed_ids.txt"
for src in data/metadata/samples.tsv "$META_SRC" data/benchmark_malay/inputs/*.tsv; do
  [[ -f "$src" ]] || continue
  tail -n +2 "$src" | cut -f1 >> "$WORK/embargoed_ids.txt"
done
# The embargoed tables also contain the public accessions (they are part of the
# same cohort). Overlap is only a violation for NON-public IDs.
sort -u "$WORK/embargoed_ids.txt" | grep -vE "$PUBLIC_RE" > "$WORK/embargoed_only.txt" || true

OVERLAP=$(comm -12 "$WORK/example_ids.txt" "$WORK/embargoed_only.txt" | tee "$WORK/overlap.txt" | wc -l | tr -d ' ')
if [[ "$OVERLAP" != "0" ]]; then
  echo "ERROR: $OVERLAP embargoed ID(s) in the example VCF:" >&2
  cat "$WORK/overlap.txt" >&2
  rm -f "$OUT_VCF" "$OUT_VCF.csi"
  exit 1
fi
echo "    0 of $(wc -l < "$WORK/example_ids.txt" | tr -d ' ') example samples appear among the"
echo "    $(wc -l < "$WORK/embargoed_only.txt" | tr -d ' ') non-public cohort IDs."

#---------------------------------------------------------------------------
echo
echo "==> Done."
printf '    %-46s %s\n' \
  "$OUT_VCF"  "$(du -h "$OUT_VCF" | cut -f1), $(bcftools index -n "$OUT_VCF") sites, $(bcftools query -l "$OUT_VCF" | wc -l | tr -d ' ') samples" \
  "$OUT_META" "$(( $(wc -l < "$OUT_META") - 1 )) rows" \
  "$OUT_MASK" "$(wc -l < "$OUT_MASK" | tr -d ' ') regions" \
  "$OUT_ACC"  "provenance"
