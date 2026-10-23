#!/usr/bin/env python3
"""Generate (or verify) the pinned-output checksum manifest for the tiny cohort.

The manifest is the byte-identity guard for Stages 0-5 plus the clonality
pass. It exists so that work on Stage 6 is bounded by evidence rather than by
assertion: if a Stage-6 change perturbs an upstream stage, `pytest tests/` says
so instead of nobody noticing.

What is pinned and what is not: see `CONTEXT-pipeline.md`. The short version is
that logs, bgzipped VCFs, their indexes, and the GDS container are excluded
because they are NOT byte-stable between runs — they embed timestamps,
hostnames, the configured output paths, or compressor state. Everything with
scientific content is pinned.

Regenerate after a DELIBERATE change to an upstream stage:

    pixi run snakemake --configfile tests/tiny_cohort/config.yaml --cores 4
    python3 tests/tiny_cohort/pin_outputs.py --write

Verify by hand (the pytest does this automatically):

    python3 tests/tiny_cohort/pin_outputs.py --check
"""

import argparse
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(HERE, "pinned_outputs.sha256")
DEFAULT_OUTPUTS = os.path.join(HERE, "outputs")

# Suffixes that are not byte-stable between runs of identical code, verified by
# running the tiny cohort three times (the committed September run plus two
# runs a month later) and diffing:
#   .log            plink/bcftools logs embed a timestamp, a hostname and the
#                   configured paths
#   .vcf.gz / .csi  BGZF output embeds an mtime, and bcftools writes a
#                   ##bcftools_command header carrying the full paths
#   .gds            binary container with its own internal metadata
UNSTABLE_SUFFIXES = (".log", ".vcf.gz", ".vcf.gz.csi", ".csi", ".gds")


def is_pinnable(rel):
    return not any(rel.endswith(s) for s in UNSTABLE_SUFFIXES)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def walk(root):
    out = {}
    for dirpath, _, filenames in os.walk(root):
        for name in filenames:
            full = os.path.join(dirpath, name)
            rel = os.path.relpath(full, root)
            if is_pinnable(rel):
                out[rel] = full
    return out


def read_manifest():
    if not os.path.exists(MANIFEST):
        return None
    pinned = {}
    with open(MANIFEST) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            digest, rel = line.split("  ", 1)
            pinned[rel] = digest
    return pinned


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs", default=DEFAULT_OUTPUTS,
                    help="tiny-cohort outputs directory to read")
    ap.add_argument("--write", action="store_true",
                    help="(re)write the manifest from --outputs")
    ap.add_argument("--check", action="store_true",
                    help="verify --outputs against the manifest")
    args = ap.parse_args()

    if not args.write and not args.check:
        ap.error("pass --write or --check")

    found = walk(args.outputs)
    if not found:
        print(f"no pinnable files under {args.outputs} — run the tiny cohort "
              f"first:\n  pixi run snakemake --configfile "
              f"tests/tiny_cohort/config.yaml --cores 4", file=sys.stderr)
        return 2

    if args.write:
        lines = [
            "# Pinned tiny-cohort outputs — the byte-identity guard for",
            "# Stages 0-5 plus the clonality pass. Regenerate with",
            "#   python3 tests/tiny_cohort/pin_outputs.py --write",
            "# only after a DELIBERATE upstream change, and say so in",
            "# CONTEXT-pipeline.md. Excluded by design: .log, .vcf.gz, .csi,",
            "# .gds — not byte-stable between runs (see the module docstring).",
        ]
        for rel in sorted(found):
            lines.append(f"{sha256(found[rel])}  {rel}")
        with open(MANIFEST, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"[pin] wrote {MANIFEST} ({len(found)} files)")
        return 0

    pinned = read_manifest()
    if pinned is None:
        print(f"no manifest at {MANIFEST}", file=sys.stderr)
        return 2
    missing = sorted(set(pinned) - set(found))
    extra = sorted(set(found) - set(pinned))
    changed = sorted(rel for rel in set(pinned) & set(found)
                     if sha256(found[rel]) != pinned[rel])
    for rel in missing:
        print(f"MISSING  {rel}")
    for rel in changed:
        print(f"CHANGED  {rel}")
    for rel in extra:
        print(f"NEW      {rel}")
    ok = not (missing or changed or extra)
    print(f"[pin] {len(pinned)} pinned, {len(changed)} changed, "
          f"{len(missing)} missing, {len(extra)} new -> "
          f"{'OK' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
