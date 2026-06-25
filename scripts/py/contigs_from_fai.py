#!/usr/bin/env python3
"""
Derive the nuclear-contig list and a contig→integer map from a FASTA index.

Both a CLI (called from a Snakemake rule to write `contig_map.tsv`) and an
importable function (`derive_contigs`, used at Snakefile parse time to expose
the nuclear-contig list to other rules without a checkpoint dance).

Replaces the hand-listed `nuclear_contigs` array in V1's config.yaml. Each
pattern in `--exclude` is a Python regex tested with `re.search` against the
contig name (column 1 of the .fai); any hit drops the contig.

CLI:
    python contigs_from_fai.py <fasta.fai> \\
        --exclude MIT --exclude API \\
        --out-list outputs/setup/nuclear_contigs.txt \\
        --out-map  outputs/setup/contig_map.tsv

Importable:
    from contigs_from_fai import derive_contigs
    contigs = derive_contigs("ref.fasta.fai", exclude=["MIT", "API"])
    # → ["ordered_PKNH_01_v2", "ordered_PKNH_02_v2", ...]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import List, Sequence, Tuple


def derive_contigs(fai_path: str | Path, exclude: Sequence[str] = ()) -> List[str]:
    """
    Read a .fai and return its contig names in file order, minus any
    matching one of `exclude` (regex, re.search).
    """
    fai = Path(fai_path)
    if not fai.exists():
        raise FileNotFoundError(f".fai not found: {fai}")
    patterns = [re.compile(p) for p in exclude]
    kept: List[str] = []
    for line in fai.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        name = line.split("\t", 1)[0]
        if any(p.search(name) for p in patterns):
            continue
        kept.append(name)
    return kept


def contig_map(contigs: Sequence[str]) -> List[Tuple[str, int]]:
    """Return [(contig, 1-based integer code), ...] in input order."""
    return [(c, i + 1) for i, c in enumerate(contigs)]


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("fai", help="Path to <fasta>.fai")
    p.add_argument(
        "--exclude",
        action="append",
        default=[],
        metavar="REGEX",
        help="Regex (re.search) matched against contig names; matches dropped. Repeatable.",
    )
    p.add_argument("--out-list", help="Write one contig name per line.")
    p.add_argument("--out-map",  help="Write tab-separated contig\\tinteger map.")
    args = p.parse_args(argv)

    contigs = derive_contigs(args.fai, args.exclude)
    if not contigs:
        print(
            f"ERROR: no contigs left after applying exclude patterns {args.exclude}",
            file=sys.stderr,
        )
        return 2

    if args.out_list:
        Path(args.out_list).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out_list).write_text("\n".join(contigs) + "\n")
    if args.out_map:
        Path(args.out_map).parent.mkdir(parents=True, exist_ok=True)
        lines = [f"{name}\t{code}" for name, code in contig_map(contigs)]
        Path(args.out_map).write_text("\n".join(lines) + "\n")
    if not args.out_list and not args.out_map:
        for c in contigs:
            print(c)
    return 0


if __name__ == "__main__":
    sys.exit(main())
