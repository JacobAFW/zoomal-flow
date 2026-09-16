#!/usr/bin/env bash
# examples/check_reference.sh — confirm the reference you downloaded is the one
# the example VCF was built against, BEFORE you spend a few minutes running the
# pipeline against the wrong assembly.
#
#   bash examples/check_reference.sh [path/to/reference.fasta]
#
# Default path is the one examples/config.yaml expects. The script builds the
# .fai if it is missing, then checks every contig name and length against the
# P. knowlesi A1-H.1 assembly as PlasmoDB distributes it.
#
# A pass means: right assembly, right contig naming, right lengths — the
# example will run. A fail tells you exactly which contig is wrong.

set -o pipefail

REF="${1:-examples/data/reference/PlasmoDB-67_PknowlesiA1H1_Genome.fasta}"

if [[ ! -f "$REF" ]]; then
  cat >&2 <<EOF
ERROR: no reference at $REF

Download the P. knowlesi A1-H.1 genome FASTA from PlasmoDB:
  https://plasmodb.org/plasmo/app/downloads/
  → release-NN → PknowlesiA1H1 → fasta/data/
  → PlasmoDB-NN_PknowlesiA1H1_Genome.fasta

Any release works: the A1-H.1 assembly has not changed (the FASTA headers
record version=2017-09-19), so the contigs below are the same in every one.
Save it to $REF — or pass your path as an argument and point
examples/config.yaml's reference.fasta at it.
EOF
  exit 1
fi

if [[ ! -f "$REF.fai" ]]; then
  echo "==> No .fai beside the FASTA; building one"
  samtools faidx "$REF" || { echo "ERROR: samtools faidx failed." >&2; exit 1; }
fi

# The expected assembly: 14 nuclear chromosomes + mitochondrion + apicoplast,
# named by ENA accession, with the lengths PlasmoDB distributes.
read -r -d '' EXPECTED <<'EOF'
LT727648	892471
LT727649	787547
LT727650	1050471
LT727651	1130816
LT727652	782346
LT727653	1080557
LT727654	1522076
LT727655	1924455
LT727656	2195390
LT727657	1516883
LT727658	2361529
LT727659	3161709
LT727660	2596912
LT727661	3267653
LT727662	5957
LT727663	34438
EOF

echo "==> Checking $REF.fai"
if diff <(printf '%s\n' "$EXPECTED") <(cut -f1,2 "$REF.fai"); then
  echo "OK — this is the assembly the example VCF was built against."
  exit 0
fi

cat >&2 <<'EOF'

FAIL — the reference does not match (differences above; "<" is expected,
">" is what you have).

Two things usually cause this:

1. Wrong assembly. The example needs P. knowlesi strain A1-H.1. If your
   lengths are close but not equal, you may have a copy with trailing N
   padding (some pipelines' references carry +99 bp on most chromosomes);
   that copy has the same coordinates but will not pass this check, and
   bcftools will warn about contig-length mismatches.

2. Different FASTA headers. samtools takes the contig name from the first
   whitespace-delimited token, so a header like

       >LT727648 | organism=Plasmodium_knowlesi_strain_A1H1 | ...

   gives "LT727648" and is fine, but one like ">ENA|LT727648|LT727648.1"
   gives "ENA|LT727648|LT727648.1" and is not. Rewrite the headers to bare
   accessions, re-index, and re-run this check:

       sed -E 's/^>.*(LT7276[0-9]{2}).*/>\1/' in.fasta > fixed.fasta
       samtools faidx fixed.fasta
EOF
exit 1
