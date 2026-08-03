#!/usr/bin/env python3
"""
render_walkthrough.py — regenerate docs/WALKTHROUGH.md from rule docstrings.

Parses every workflow/rules/*.smk file in stage order, extracts each rule's
one-line summary + the WHAT / WHY / TUNABLES / OUTPUT / TRY block from its
docstring, and renders `docs/WALKTHROUGH.md`. The Markdown is *raw* — this
script is the single source of truth (DESIGN §6); the doc is committed as
Markdown, not knitted.

Two modes:

  --write   Regenerate WALKTHROUGH.md in place.
  --check   Regenerate to a temp file and diff against the committed copy;
            exit non-zero if they differ (for the `walkthrough` staleness
            check Snakemake target).

CLI:
    python scripts/py/render_walkthrough.py \\
        --config config/cohort.example.yaml \\
        --write

The header stamps the git commit + the config file it was rendered from
so a reader knows exactly which state the resolved TUNABLES came from.
"""

from __future__ import annotations

import argparse
import difflib
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import yaml


HERE     = Path(__file__).resolve().parent
AGNOSTIC = HERE.parent.parent
RULES_DIR = AGNOSTIC / "workflow" / "rules"
WALKTHROUGH_PATH = AGNOSTIC / "docs" / "WALKTHROUGH.md"


# --------------------------------------------------------------------------
# Data structures
# --------------------------------------------------------------------------

@dataclass
class Rule:
    name: str
    file: str
    summary: str = ""
    what: str = ""
    why: str = ""
    tunables: List[str] = field(default_factory=list)
    output: str = ""
    try_: str = ""


# --------------------------------------------------------------------------
# Docstring parser
# --------------------------------------------------------------------------

RULE_HEADER_RE = re.compile(
    r'^(?:rule|checkpoint)\s+([A-Za-z_][A-Za-z0-9_]*)\s*:\s*$'
)


def parse_docstring(body: str) -> Dict[str, str]:
    """
    Parse a rule's docstring body (the string content between the triple
    quotes) into a dict of sections. Recognises the fields WHAT, WHY,
    TUNABLES, OUTPUT, TRY. The first non-blank paragraph before any
    section header is the summary.

    Section headers are case-sensitive lines of the shape "FIELD:" at the
    start of a line (after leading whitespace strip).
    """
    lines = [ln.rstrip() for ln in body.strip("\n").splitlines()]
    sections: Dict[str, List[str]] = {"summary": []}
    current = "summary"
    section_re = re.compile(r'^\s*(WHAT|WHY|TUNABLES|OUTPUT|TRY)\s*:\s*(.*)$')
    for ln in lines:
        stripped = ln.strip()
        m = section_re.match(ln)
        if m:
            current = m.group(1).lower()
            first_val = m.group(2).strip()
            sections.setdefault(current, [])
            if first_val:
                sections[current].append(first_val)
            continue
        # Continuation line — strip a common indent (crudely: everything
        # after the first non-space, so multi-line WHY blocks flow).
        sections.setdefault(current, [])
        sections[current].append(stripped)
    # Join and collapse leading/trailing blanks in each section
    return {k: "\n".join(v).strip() for k, v in sections.items()}


def parse_rules_from_smk(smk_path: Path) -> List[Rule]:
    """
    Walk the .smk file line by line, find rule/checkpoint headers, then
    look for a triple-quoted docstring immediately inside the block.
    """
    text = smk_path.read_text()
    lines = text.splitlines()
    rules: List[Rule] = []
    i = 0
    while i < len(lines):
        m = RULE_HEADER_RE.match(lines[i])
        if not m:
            i += 1
            continue
        name = m.group(1)
        # scan for opening """ on the next few lines
        j = i + 1
        while j < len(lines) and not lines[j].strip().startswith(('"""', "'''")):
            # If we hit the next rule/checkpoint or top-level statement
            # without seeing a docstring, give up on this rule.
            if RULE_HEADER_RE.match(lines[j]):
                break
            j += 1
        if j >= len(lines) or not lines[j].strip().startswith(('"""', "'''")):
            i = j
            continue
        opener = lines[j].strip()[:3]
        # Collect from after opener until closing triple-quote
        doc_lines: List[str] = []
        # rest of the opening line, if any
        rest = lines[j].strip()[3:]
        if rest.endswith(opener) and len(rest) >= 3:
            doc_lines.append(rest[: -3])
            k = j + 1
        else:
            if rest:
                doc_lines.append(rest)
            k = j + 1
            while k < len(lines):
                if opener in lines[k]:
                    doc_lines.append(lines[k].split(opener, 1)[0])
                    break
                doc_lines.append(lines[k])
                k += 1
        body = "\n".join(doc_lines)
        sec = parse_docstring(body)

        rule = Rule(name=name, file=str(smk_path.relative_to(AGNOSTIC)))
        rule.summary = sec.get("summary", "").strip()
        rule.what    = sec.get("what", "").strip()
        rule.why     = sec.get("why", "").strip()
        rule.output  = sec.get("output", "").strip()
        rule.try_    = sec.get("try", "").strip()
        # TUNABLES: split on comma / whitespace; keep any dotted config
        # paths intact.
        tun_block = sec.get("tunables", "").strip()
        if tun_block and tun_block.lower() not in {"(none)", "none"}:
            # Strip parenthesised commentary that would otherwise get parsed
            # as bogus dotted keys (e.g. "cohort.vcf (bgzipped + indexed)"
            # → drop the "(bgzipped + indexed)" part).
            tun_block_clean = re.sub(r'\([^)]*\)', '', tun_block)
            tokens = re.split(r'[,\s]+', tun_block_clean)
            rule.tunables = []
            for tok in tokens:
                tok = tok.strip().rstrip(",.:")
                if not tok:
                    continue
                # Only accept dotted or lowercase-underscore identifiers
                if re.fullmatch(r'[A-Za-z_][A-Za-z0-9_]*(\.[A-Za-z0-9_]+)*', tok):
                    rule.tunables.append(tok)
        rules.append(rule)
        i = k + 1
    return rules


# --------------------------------------------------------------------------
# Config resolver
# --------------------------------------------------------------------------

def _merge_dicts(base: dict, override: dict) -> dict:
    """Shallow-recursive merge (override wins)."""
    out = dict(base)
    for k, v in override.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = _merge_dicts(out[k], v)
        else:
            out[k] = v
    return out


def load_config(config_path: Path) -> dict:
    """
    Load a cohort or template config, layered on top of the default
    config.yaml (Snakemake's --configfile semantics do the same).
    """
    default_cfg = AGNOSTIC / "config" / "config.yaml"
    with default_cfg.open() as fh:
        cfg = yaml.safe_load(fh) or {}
    if config_path.resolve() != default_cfg.resolve():
        with config_path.open() as fh:
            override = yaml.safe_load(fh) or {}
        cfg = _merge_dicts(cfg, override)
    return cfg


def resolve_dotted(cfg: dict, path: str):
    """
    Resolve `structure.min_maf` → cfg['structure']['min_maf']. Returns
    "(unset)" when the path doesn't resolve.
    """
    ref = cfg
    for part in path.split("."):
        if isinstance(ref, dict) and part in ref:
            ref = ref[part]
        else:
            return None
    return ref


# --------------------------------------------------------------------------
# Renderer
# --------------------------------------------------------------------------

STAGE_ORDER = [
    ("00_setup",            "Stage 0 — Setup (VCF sample list, contig map, metadata validation)"),
    ("01_qc",               "Stage 1 — QC (common)"),
    ("qc_wgs",              "Stage 1 — QC (WGS path)"),
    ("qc_microhap",         "Stage 1 — QC (microhap seam stub)"),
    ("02_moi",              "Stage 2 — MOI / Fws (common)"),
    ("moi_wgs",             "Stage 2 — MOI / Fws (WGS path)"),
    ("moi_microhap",        "Stage 2 — MOI / Fws (microhap seam stub)"),
    ("03_structure",        "Stage 3 — Structure (analysis core)"),
    ("structure_prep_wgs",  "Stage 3 — Structure prep (WGS path)"),
    ("structure_prep_microhap", "Stage 3 — Structure prep (microhap seam stub)"),
    ("03_figures",          "Stage 3b — Structure figures"),
    ("04_ibd",              "Stage 4 — IBD (per-cluster hmmIBD + clonal detection)"),
    ("04_figures",          "Stage 4b — IBD figures"),
    ("06_selection",        "Stage 6 — Selection (rehh iHS)"),
]


def render_walkthrough(cfg: dict, config_path: Path, commit_hash: str,
                       rules_by_file: Dict[str, List[Rule]]) -> str:
    out: List[str] = []
    out.append("# ZOOMAL-Flow — Walkthrough")
    out.append("")
    out.append(f"*Generated from `workflow/rules/*.smk` docstrings by "
               f"`scripts/py/render_walkthrough.py`.*  ")
    out.append("*Do NOT hand-edit — regenerate via `snakemake walkthrough` or "
               "`python scripts/py/render_walkthrough.py --write`.*")
    out.append("")
    out.append(f"- Config: `{config_path.relative_to(AGNOSTIC) if config_path.is_relative_to(AGNOSTIC) else config_path}`")
    out.append(f"- Commit: `{commit_hash}`")
    out.append("")
    out.append("Each rule below carries a WHAT/WHY block, its resolved TUNABLES "
               "(current values from the config above), its OUTPUT path(s), and "
               "a TRY suggestion — a concrete experiment you can run by editing "
               "the config and re-invoking that stage's target.")
    out.append("")
    out.append("---")
    out.append("")

    seen: set = set()
    for basename, heading in STAGE_ORDER:
        smk_key = f"workflow/rules/{basename}.smk"
        rules = rules_by_file.get(smk_key)
        if not rules:
            continue
        seen.add(smk_key)
        out.append(f"## {heading}")
        out.append("")
        out.append(f"Source: `{smk_key}`")
        out.append("")
        stage_targets = [r.output for r in rules if r.output and "outputs/" in r.output]
        if stage_targets:
            out.append(f"Run this stage alone:")
            out.append("")
            out.append("```bash")
            out.append(f"# any rule below can be invoked as `snakemake <rule_name>`,")
            out.append(f"# or reach the stage's core outputs directly:")
            for r in rules[:3]:
                if r.output and r.output.startswith("{outputs}"):
                    out.append(f"snakemake {r.output.splitlines()[0]}")
            out.append("```")
            out.append("")
        for r in rules:
            out.append(f"### `{r.name}`")
            out.append("")
            if r.summary:
                out.append(r.summary)
                out.append("")
            if r.what:
                out.append(f"**WHAT.** {r.what}")
                out.append("")
            if r.why:
                out.append(f"**WHY.** {r.why}")
                out.append("")
            if r.tunables:
                lines = ["**TUNABLES.**", ""]
                for t in r.tunables:
                    val = resolve_dotted(cfg, t)
                    if val is None:
                        lines.append(f"- `{t}`: *(not set in this config)*")
                    else:
                        lines.append(f"- `{t}` = `{val!r}`")
                out.extend(lines)
                out.append("")
            else:
                out.append("**TUNABLES.** *(none)*")
                out.append("")
            if r.output:
                out.append(f"**OUTPUT.** `{r.output}`")
                out.append("")
            if r.try_:
                out.append(f"**TRY.** {r.try_}")
                out.append("")
            out.append("---")
            out.append("")

    # Any rules files we didn't include in STAGE_ORDER
    for key, rules in rules_by_file.items():
        if key in seen:
            continue
        out.append(f"## {key}")
        out.append("")
        for r in rules:
            out.append(f"### `{r.name}`")
            out.append("")
            if r.summary:
                out.append(r.summary)
                out.append("")
        out.append("---")
        out.append("")

    return "\n".join(out).rstrip() + "\n"


def _git_commit() -> str:
    try:
        h = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=AGNOSTIC, text=True
        ).strip()
        dirty = subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=AGNOSTIC, text=True
        ).strip()
        return f"{h}{'-dirty' if dirty else ''}"
    except Exception:
        return "(unknown)"


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--config", type=Path,
                   default=AGNOSTIC / "config" / "cohort.example.yaml",
                   help="Config file for TUNABLES resolution.")
    p.add_argument("--write",  action="store_true", help="Regenerate WALKTHROUGH.md in place.")
    p.add_argument("--check",  action="store_true", help="Fail if committed WALKTHROUGH.md is stale.")
    args = p.parse_args(argv)

    if not args.write and not args.check:
        args.write = True  # default

    cfg = load_config(args.config)
    commit = _git_commit()

    rules_by_file: Dict[str, List[Rule]] = {}
    for smk in sorted(RULES_DIR.glob("*.smk")):
        rules_by_file[str(smk.relative_to(AGNOSTIC))] = parse_rules_from_smk(smk)

    rendered = render_walkthrough(cfg, args.config, commit, rules_by_file)

    WALKTHROUGH_PATH.parent.mkdir(parents=True, exist_ok=True)

    if args.check:
        if not WALKTHROUGH_PATH.exists():
            print(f"WALKTHROUGH.md missing: run --write first.", file=sys.stderr)
            return 2
        current = WALKTHROUGH_PATH.read_text()
        # Strip the header stamp lines (config + commit) from both sides
        # before diffing — those depend on which config was in play and
        # the working-tree state, so they'd trip the check even when the
        # content matches.
        def _strip_stamp(t: str) -> str:
            return "\n".join(
                ln for ln in t.splitlines()
                if not ln.startswith("- Config:") and not ln.startswith("- Commit:")
            )
        if _strip_stamp(current) == _strip_stamp(rendered):
            print("WALKTHROUGH.md up to date.")
            return 0
        # Show a compact diff on the stamp-stripped content
        diff = list(difflib.unified_diff(
            _strip_stamp(current).splitlines(keepends=True),
            _strip_stamp(rendered).splitlines(keepends=True),
            fromfile="committed",
            tofile="regenerated",
            n=2,
        ))
        sys.stderr.write("".join(diff[:80]))
        if len(diff) > 80:
            sys.stderr.write(f"... ({len(diff) - 80} more diff lines suppressed) ...\n")
        print("WALKTHROUGH.md is stale — regenerate with `snakemake walkthrough`.",
              file=sys.stderr)
        return 1

    WALKTHROUGH_PATH.write_text(rendered)
    print(f"Wrote {WALKTHROUGH_PATH.relative_to(AGNOSTIC)} ({len(rendered)} bytes).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
