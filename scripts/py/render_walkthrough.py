#!/usr/bin/env python3
"""
render_walkthrough.py — regenerate docs/WALKTHROUGH.md from rule docstrings.

Parses every workflow/rules/*.smk file in stage order, extracts each rule's
one-line summary + the WHAT / WHY / TUNABLES / OUTPUT / TRY block from its
docstring, and renders two things from that one pass:

  docs/WALKTHROUGH.md   the whole pipeline as one linear read
  docs/steps/           one page per stage, with the run mechanics for
                        stepping through the pipeline by hand

The Markdown is *raw* — this script is the single source of truth (DESIGN §6);
the docs are committed as Markdown, not knitted. Generating both here is what
guarantees they cannot drift apart.

Two modes:

  --write   Regenerate docs/WALKTHROUGH.md and docs/steps/ in place.
  --check   Re-render both and diff against the committed copies; exit
            non-zero if any differ (for the `walkthrough` staleness check
            Snakemake target).

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
STEPS_DIR        = AGNOSTIC / "docs" / "steps"


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
    # Structural fields — used to derive each stage's terminal rules, so the
    # `--until` command on a step page is read off the rules themselves rather
    # than hand-maintained.
    dep_rules: List[str] = field(default_factory=list)      # rules.X.output / checkpoints.X
    input_paths: List[str] = field(default_factory=list)    # literal path templates
    output_paths: List[str] = field(default_factory=list)


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


DIRECTIVE_RE = re.compile(
    r'^\s{4}(input|output|log|params|threads|message|shell|run|resources|'
    r'wildcard_constraints|benchmark|priority|group|conda|container|script|'
    r'notebook|cache|default_target|retries|localrule|envmodules|shadow|'
    r'template_engine|handover)\s*:'
)

STRING_LITERAL_RE = re.compile(r'(?:f|rf|fr)?"([^"\n]*)"|(?:f|rf|fr)?\'([^\'\n]*)\'')
RULE_REF_RE = re.compile(r'\b(?:rules|checkpoints)\.([A-Za-z_][A-Za-z0-9_]*)\.')


def _directive_block(body_lines: List[str], name: str) -> str:
    """
    Return the text of one directive block (e.g. `input:`) from a rule body,
    i.e. every line after `    <name>:` up to the next directive at the same
    indent. Returns "" when the rule has no such directive.
    """
    out: List[str] = []
    collecting = False
    for ln in body_lines:
        m = DIRECTIVE_RE.match(ln)
        if m:
            if collecting:
                break
            collecting = (m.group(1) == name)
            # keep anything on the same line after the colon
            rest = ln.split(":", 1)[1] if collecting else ""
            if rest.strip():
                out.append(rest)
            continue
        if collecting:
            out.append(ln)
    return "\n".join(out)


def _string_literals(block: str) -> List[str]:
    """Every single- or double-quoted literal in a block, in order."""
    vals: List[str] = []
    for m in STRING_LITERAL_RE.finditer(block):
        vals.append(m.group(1) if m.group(1) is not None else m.group(2))
    return vals


def _path_literals(block: str) -> List[str]:
    """
    Literals that look like file paths — the ones worth comparing between one
    rule's output and another's input. Drops bare dict keys and flags.
    """
    return [v for v in _string_literals(block) if "/" in v]


def _rule_refs(block: str) -> List[str]:
    """Rule names referenced as `rules.X.output` / `checkpoints.X.get(...)`."""
    seen: List[str] = []
    for m in RULE_REF_RE.finditer(block):
        if m.group(1) not in seen:
            seen.append(m.group(1))
    return seen


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

        # Structural scan: everything from just after the docstring to the
        # next rule header is this rule's body. We read `input:`/`output:` out
        # of it so the stage pages can work out which rules a stage ends on.
        b = k + 1
        body_lines: List[str] = []
        while b < len(lines) and not RULE_HEADER_RE.match(lines[b]):
            body_lines.append(lines[b])
            b += 1
        in_block  = _directive_block(body_lines, "input")
        out_block = _directive_block(body_lines, "output")
        rule.dep_rules    = _rule_refs(in_block)
        rule.input_paths  = _path_literals(in_block)
        rule.output_paths = _path_literals(out_block)

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
    ("05_introgression",    "Stage 5 — Introgression (pairwise density-cloud detection)"),
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


# --------------------------------------------------------------------------
# Per-step pages (docs/steps/)
# --------------------------------------------------------------------------
#
# Same parse as WALKTHROUGH.md, split one file per stage and given the run
# mechanics WALKTHROUGH.md lacks: how to run *up to* a stage, stop, look at
# what came out, change something, and re-run just that step.
#
# Rule and target names are read off the rules, never typed here, so a renamed
# rule changes these pages on the next `snakemake walkthrough`.

# The config shown in the example commands. A placeholder for the reader's own
# cohort config — not the config TUNABLES were resolved from.
STEP_CMD_CONFIG = "config/config.yaml"
STEP_CMD_CORES = 4


@dataclass
class Stage:
    slug: str                      # file stem, e.g. "03_structure"
    title: str
    files: List[str]               # .smk basenames, in order
    purpose: str                   # one line for the index table
    does: str                      # what this stage does
    needs: str                     # what it needs from the stage before
    deep_doc: Optional[str] = None  # path to a deeper design doc, if any
    # Some stages are reached by a named target rather than by --until: the
    # report is deliberately NOT part of `rule all`, so --until render_report
    # cannot resolve it. Where this is set, the run command uses it directly.
    target: Optional[str] = None


STAGES: List[Stage] = [
    Stage(
        slug="00_setup", purpose="Preflight: sample roster, contig map, metadata validated against the role map.", title="Stage 0 — Setup",
        files=["00_setup"],
        does="Reads the cohort VCF's sample roster, works out which contigs are "
             "nuclear, and validates your `samples.tsv` against the role map in "
             "the config. Nothing here changes any genotype — it is the preflight "
             "that lets every later stage assume its inputs are sane.",
        needs="Nothing but your three inputs: the bgzipped VCF, the reference "
              "FASTA, and the metadata table the config points at.",
    ),
    Stage(
        slug="01_qc", purpose="Filter to biallelic PASS SNPs on nuclear contigs, controls and masked regions removed.", title="Stage 1 — QC",
        files=["01_qc", "qc_wgs", "qc_microhap"],
        does="Drops controls and non-nuclear contigs, masks the regions listed in "
             "the config, and keeps biallelic PASS SNPs above the MAC floor. Ends "
             "at the QC seam — one filtered VCF that both the WGS and microhap "
             "paths converge on.",
        needs="Stage 0's sample list, contig map, and validated metadata.",
    ),
    Stage(
        slug="02_moi", purpose="Fws per sample — how mixed each infection is — and the high-MOI exclusion list.", title="Stage 2 — MOI / Fws",
        files=["02_moi", "moi_wgs", "moi_microhap"],
        does="Estimates within-host diversity per sample. Fws near 1 means a "
             "single dominant genotype; lower values mean a mixed (polyclonal) "
             "infection. Produces the Fws table, the density figures, and the "
             "exclusion list that keeps highly-mixed samples out of Stage 3.",
        needs="The QC seam VCF from Stage 1.",
    ),
    Stage(
        slug="03_structure", purpose="ADMIXTURE, PCA and NJ tree; picks K and assigns every sample to a cluster.", title="Stage 3 — Population structure",
        files=["03_structure", "structure_prep_wgs", "structure_prep_microhap",
               "03_figures"],
        does="Turns the filtered VCF into a PLINK bfile, LD-prunes it, then runs "
             "ADMIXTURE across a range of K, picks K by cross-validation error, "
             "assigns each sample to a cluster, and draws the PCA, NJ tree and "
             "ancestry bars. This is the stage whose cluster assignments Stages 4 "
             "and 5 are built on.",
        needs="Stage 1's QC seam VCF and Stage 2's high-MOI exclusion list.",
        deep_doc="docs/clonality.md",
    ),
    Stage(
        slug="04_ibd", purpose="Pairwise IBD within each cluster, clonal groups, and the two sample-set keep-lists.", title="Stage 4 — IBD and clonality",
        files=["04_ibd", "04b_clonality", "04_figures"],
        does="Runs hmmIBD within each cluster to get pairwise identity-by-descent, "
             "calls clonal groups above the IBD threshold, and writes the two "
             "sample-set keep-lists (`full` and `unique`) that every frequency-based "
             "stage is then built twice against.",
        needs="Stage 3's cluster assignments and cleaned bfile.",
        deep_doc="docs/clonality.md",
    ),
    Stage(
        slug="05_introgression", purpose="Windows that look introgressed between each pair of clusters.", title="Stage 5 — Introgression",
        files=["05_introgression"],
        does="Scans windows for genotype patterns that look introgressed between "
             "each pair of clusters, applies the artifact masks and the support "
             "floor, and aggregates the surviving windows.",
        needs="Stage 3's clusters and Stage 4's sample-set keep-lists.",
        deep_doc="docs/introgression_analysis_spec.md",
    ),
    Stage(
        slug="06_selection", purpose="rehh iHS scan, per case/control model declared in the config.", title="Stage 6 — Selection (iHS)",
        files=["06_selection"],
        does="Runs the rehh iHS scan for each case/control model declared in the "
             "config. With no models declared the stage produces nothing, by "
             "design — an iHS scan needs a cohort-specific hypothesis.",
        needs="Stage 3's phased/cleaned genotypes and cluster labels.",
    ),
    Stage(
        slug="07_report", purpose="Full vs de-clonalized comparison tables, then the self-contained HTML report.", title="Stage 7 — Comparison tables and the report",
        files=["07_declonalization", "99_report"],
        target="report",
        does="Builds the full-vs-de-clonalized comparison tables, then renders the "
             "single self-contained HTML report from whatever the run actually "
             "produced. Stages that did not run are reported as such rather than "
             "omitted.",
        needs="Everything above — the report target pulls the whole pipeline.",
    ),
]


def rule_seam(r: Rule) -> Optional[str]:
    """
    Which seam path a rule belongs to, or None for the shared ones. The WGS and
    microhap paths are alternatives: `cohort.input_type` decides which set of
    rule files the Snakefile includes, so the other path's rules are not in the
    DAG at all and must stay out of any command we print.
    """
    stem = Path(r.file).stem
    if stem.endswith("_wgs"):
        return "wgs"
    if stem.endswith("_microhap"):
        return "microhap"
    return None


def active_rules(rules: List[Rule], cfg: dict) -> List[Rule]:
    """Rules that are actually in the DAG for this config's input_type."""
    want = (resolve_dotted(cfg, "cohort.input_type") or "wgs")
    return [r for r in rules if rule_seam(r) in (None, want)]

def stage_rules(stage: Stage, rules_by_file: Dict[str, List[Rule]]) -> List[Rule]:
    """Every parsed rule belonging to a stage, in file then declaration order."""
    out: List[Rule] = []
    for base in stage.files:
        out.extend(rules_by_file.get(f"workflow/rules/{base}.smk", []))
    return out


def terminal_rules(rules: List[Rule]) -> List[Rule]:
    """
    The rules a stage *ends* on: those whose output no other rule in the same
    stage consumes. `--until` these and Snakemake runs the whole stage and
    stops.

    Consumption is detected two ways: an explicit `rules.X.output` reference in
    another rule's input, and a literal path template appearing in both one
    rule's output and another's input. Inputs built by helper functions are
    invisible to both, so a rule can be reported terminal when something later
    in the stage does in fact consume it. That over-reports rather than
    under-reports, and `--until` with an extra rule still runs the whole stage
    — the failure mode is a slightly longer command, not a wrong one.
    """
    names = {r.name for r in rules}
    consumed: set = set()
    # path -> every rule declaring it. A list, not a single name: the WGS and
    # microhap seam paths deliberately declare the same output paths, so
    # last-writer-wins would hide whichever producer was parsed first and
    # report it as terminal.
    produces: Dict[str, List[str]] = {}
    for r in rules:
        for pth in r.output_paths:
            produces.setdefault(pth, []).append(r.name)
    for r in rules:
        for dep in r.dep_rules:
            if dep in names and dep != r.name:
                consumed.add(dep)
        for pth in r.input_paths:
            for owner in produces.get(pth, []):
                if owner != r.name:
                    consumed.add(owner)
    term = [r for r in rules if r.name not in consumed]
    return term or rules


SMK_CONST_RE = re.compile(
    r'^([A-Z][A-Z0-9_]*)\s*=\s*(?:f|rf|fr)?"([^"\n]*)"', re.M
)


def smk_constants() -> Dict[str, str]:
    """
    Module-level string constants defined in the rule files, e.g.
    `CLONALITY_DIR = f"{PATHS['outputs']}/clonality"`. Rules interpolate these
    into their output paths, so resolving a path template means resolving
    these first. Read from the sources, so a renamed or repointed constant
    follows automatically.
    """
    consts: Dict[str, str] = {}
    for smk in sorted(RULES_DIR.glob("*.smk")):
        for name, val in SMK_CONST_RE.findall(smk.read_text()):
            consts.setdefault(name, val)
    return consts

def resolve_output_template(tpl: str, cfg: dict) -> Optional[str]:
    """
    Turn a raw `output:` template from the source into a real path.

    The source text is an f-string body, so it looks like
    `{PATHS['outputs']}/ibd/clonal_clusters.tsv` with wildcards doubled as
    `{{cluster}}`. Substitute the path roots from the config; return None if
    any wildcard survives, because a wildcard path cannot be named as a target.
    """
    paths = cfg.get("paths") or {}
    out = tpl
    # Expand any module-level constant first — it may itself contain a
    # {PATHS[...]} reference, which the loop below then resolves.
    for name, val in smk_constants().items():
        out = out.replace("{%s}" % name, val)
    for key in ("outputs", "logs", "reports"):
        val = paths.get(key)
        if val is None:
            continue
        out = out.replace("{PATHS['%s']}" % key, str(val))
        out = out.replace('{PATHS["%s"]}' % key, str(val))
    # `{sampleset}` is the one wildcard whose values the config already fixes
    # (`full`, plus `unique` when clonality.declonalize is on), so it can be
    # expanded here. `{cluster}` cannot — it comes out of the assign_clusters
    # checkpoint and is not knowable before the run — and `{K}`/`{model}` are
    # per-item outputs rather than a stage endpoint. Those are left to fail the
    # brace test below and are skipped.
    if "{{sampleset}}" in out:
        out = out.replace("{{sampleset}}", "full")
    if "{{model}}" in out:
        models = resolve_dotted(cfg, "selection.models") or []
        names = [m.get("name") for m in models if isinstance(m, dict) and m.get("name")]
        if not names:
            return None          # no models declared → the stage runs nothing
        out = out.replace("{{model}}", str(names[0]))
    if "{" in out or "}" in out:
        return None
    return out


def stage_targets(rules: List[Rule], cfg: dict) -> List[str]:
    """
    Concrete output files that stand for "this stage is done".

    Named file targets rather than `--until <rule>`: `--until` cannot be used
    once a checkpoint is involved. Stage 4's clusters come from the
    `assign_clusters` checkpoint, and `--until` a rule downstream of it
    deadlocks with "Out of jobs ready to be started, but not all files built
    yet" (snakemake#823). Asking for the files instead resolves the checkpoint
    normally. Wildcard outputs are skipped — there is no single path to name.
    """
    targets: List[str] = []
    for r in terminal_rules(rules):
        for tpl in r.output_paths:
            got = resolve_output_template(tpl, cfg)
            if got and got not in targets:
                targets.append(got)
                break           # one representative file per terminal rule
    return targets

def _wrap_targets(names: List[str], indent: str = "    ") -> str:
    """Lay a target list out over continuation lines so it stays readable."""
    if len(names) <= 3:
        return " ".join(names)
    lines = []
    for i in range(0, len(names), 3):
        lines.append(" ".join(names[i:i + 3]))
    return (" \\\n" + indent).join(lines)


def render_rule_block(r: Rule, cfg: dict, out: List[str]) -> None:
    """The WHAT/WHY/TUNABLES/OUTPUT/TRY block — identical to WALKTHROUGH.md."""
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
        out.append("**TUNABLES.**")
        out.append("")
        for t in r.tunables:
            val = resolve_dotted(cfg, t)
            if val is None:
                out.append(f"- `{t}`: *(not set in this config)*")
            else:
                out.append(f"- `{t}` = `{val!r}`")
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


def render_step_page(stage: Stage, cfg: dict, config_path: Path,
                     commit_hash: str,
                     rules_by_file: Dict[str, List[Rule]],
                     position: int, total: int,
                     prev_stage: Optional[Stage],
                     next_stage: Optional[Stage]) -> str:
    rules = stage_rules(stage, rules_by_file)
    # Commands are built from the seam path this config actually runs; the
    # other path's rules are documented below but are not in the DAG, and
    # naming one in --until fails with MissingRuleException.
    runnable = active_rules(rules, cfg)
    term = terminal_rules(runnable)
    file_targets = stage_targets(runnable, cfg)
    input_type = resolve_dotted(cfg, "cohort.input_type") or "wgs"
    inactive = [r for r in rules if r not in runnable]
    cfgp = (config_path.relative_to(AGNOSTIC)
            if config_path.is_relative_to(AGNOSTIC) else config_path)

    o: List[str] = []
    o.append(f"# {stage.title}")
    o.append("")
    o.append("*Generated by `scripts/py/render_walkthrough.py` — do NOT hand-edit. "
             "Regenerate with `snakemake walkthrough`.*")
    o.append("")
    o.append(f"- Step {position} of {total} · "
             f"[index](README.md)"
             + (f" · previous: [{prev_stage.title}]({prev_stage.slug}.md)" if prev_stage else "")
             + (f" · next: [{next_stage.title}]({next_stage.slug}.md)" if next_stage else ""))
    o.append(f"- Config the TUNABLES below were resolved from: `{cfgp}`")
    o.append(f"- Commit: `{commit_hash}`")
    o.append("")
    o.append("## What this stage does")
    o.append("")
    o.append(stage.does)
    o.append("")
    o.append(f"**Needs from before it.** {stage.needs}")
    o.append("")
    if stage.deep_doc:
        o.append(f"**Deeper reading.** [`{stage.deep_doc}`](../{Path(stage.deep_doc).name}) "
                 "covers the design behind this stage in detail.")
        o.append("")
    o.append(f"Rules: " + ", ".join(f"`{r.name}`" for r in rules) if rules
             else "Rules: *(none parsed for this stage)*")
    o.append("")
    o.append("## Running it")
    o.append("")
    o.append("Replace the config with your own, and `--cores` with what your "
             "machine has. See [the index](README.md#stepping-through-manually) "
             "for how stopping and resuming works.")
    o.append("")
    if not stage.target and not file_targets:
        o.append("**This stage produces nothing with the config above**, so there "
                 "is no target to ask for. That is the stage's configured "
                 "behaviour, not a failure — see the rules below for what "
                 "switches it on. Once it is configured, its outputs become "
                 "nameable targets and this page will show them.")
        o.append("")
        o.append("---")
        o.append("")
        o.append("## The steps")
        o.append("")
        for r in rules:
            o.append(f"### `{r.name}`")
            o.append("")
            o.append(f"<sub>`{r.file}`</sub>")
            o.append("")
            render_rule_block(r, cfg, o)
            o.append("---")
            o.append("")
        if next_stage:
            o.append(f"Next: [{next_stage.title}]({next_stage.slug}.md)")
        else:
            o.append("That is the last stage. Back to [the index](README.md).")
        o.append("")
        return "\n".join(o).rstrip() + "\n"

    o.append("Run this stage and everything it depends on, then stop:")
    o.append("")
    o.append("```bash")
    if stage.target:
        o.append(f"pixi run snakemake {stage.target} \\")
        o.append(f"    --configfile {STEP_CMD_CONFIG} --cores {STEP_CMD_CORES}")
    elif file_targets:
        o.append(f"pixi run snakemake --configfile {STEP_CMD_CONFIG} \\")
        o.append(f"    --cores {STEP_CMD_CORES} \\")
        o.append(f"    {_wrap_targets(file_targets, indent='    ')}")
    else:
        o.append(f"pixi run snakemake --configfile {STEP_CMD_CONFIG} \\")
        o.append(f"    --cores {STEP_CMD_CORES} \\")
        o.append(f"    --until {_wrap_targets([r.name for r in term], indent='    ')}")
    o.append("```")
    o.append("")
    if stage.target:
        o.append(f"(`{stage.target}` is a named target rather than part of "
                 "`rule all`, so it is given directly instead of via `--until`.)")
    else:
        o.append("Those are this stage's own output files. Asking for files "
                 "rather than `--until <rule>` is deliberate: `--until` cannot "
                 "cross the `assign_clusters` checkpoint and deadlocks on any "
                 "stage downstream of it. Snakemake builds only the part of the "
                 "graph these files need, so work through the stages in order "
                 "and each command picks up where the last one stopped.")
        o.append("")
        o.append("The paths come from `paths.outputs` in the config above — if "
                 "yours differs, they move with it.")
    o.append("")
    o.append("Re-run this whole stage — after editing a threshold in the config, "
             "say — and everything downstream of it:")
    o.append("")
    o.append("```bash")
    o.append(f"pixi run snakemake --configfile {STEP_CMD_CONFIG} \\")
    o.append(f"    --cores {STEP_CMD_CORES} \\")
    o.append(f"    --forcerun {_wrap_targets([r.name for r in runnable], indent='    ')}")
    o.append("```")
    o.append("")
    o.append("That re-runs the stage *and everything downstream of it*, which "
             "is usually what a changed threshold requires. To redo a single "
             "step without touching the rest, each rule below carries a "
             "`--force <file>` command that does only that job.")
    o.append("")
    if inactive:
        o.append(f"This config sets `cohort.input_type: {input_type}`, so "
                 + ", ".join(f"`{r.name}`" for r in inactive)
                 + " (the other seam path) is documented below but is not in the "
                 "DAG — the commands above leave it out, and naming it would fail "
                 "with `MissingRuleException`.")
        o.append("")
    o.append("---")
    o.append("")
    o.append("## The steps")
    o.append("")
    for r in rules:
        o.append(f"### `{r.name}`")
        o.append("")
        o.append(f"<sub>`{r.file}`</sub>")
        o.append("")
        render_rule_block(r, cfg, o)
        if r in runnable:
            own = None
            for tpl in r.output_paths:
                own = resolve_output_template(tpl, cfg)
                if own:
                    break
            if own:
                o.append("**Re-run just this step** — forces this one job and "
                         "nothing else:")
                o.append("")
                o.append("```bash")
                o.append(f"pixi run snakemake --configfile {STEP_CMD_CONFIG} "
                         f"--cores {STEP_CMD_CORES} --force {own}")
                o.append("```")
                o.append("")
                o.append("**Re-run it and everything downstream** — what you want "
                         "after changing a setting this step depends on, because "
                         "the results drawn from it are now stale too:")
            else:
                o.append("This rule's output carries a wildcard, so there is no "
                         "single file to force. Re-run it — and everything "
                         "downstream — by rule name:")
            o.append("")
            o.append("```bash")
            o.append(f"pixi run snakemake --configfile {STEP_CMD_CONFIG} "
                     f"--cores {STEP_CMD_CORES} --forcerun {r.name}")
            o.append("```")
            o.append("")
        else:
            o.append(f"*Not in the DAG for `cohort.input_type: {input_type}` — "
                     f"this is the other seam path. Switch `cohort.input_type` to "
                     f"`{rule_seam(r)}` to run it.*")
            o.append("")
        o.append("---")
        o.append("")
    if next_stage:
        o.append(f"Next: [{next_stage.title}]({next_stage.slug}.md)")
    else:
        o.append("That is the last stage. Back to [the index](README.md).")
    o.append("")
    return "\n".join(o).rstrip() + "\n"


def render_steps_index(cfg: dict, config_path: Path, commit_hash: str,
                       rules_by_file: Dict[str, List[Rule]]) -> str:
    cfgp = (config_path.relative_to(AGNOSTIC)
            if config_path.is_relative_to(AGNOSTIC) else config_path)
    o: List[str] = []
    o.append("# Running ZOOMAL-Flow one step at a time")
    o.append("")
    o.append("*Generated by `scripts/py/render_walkthrough.py` — do NOT hand-edit. "
             "Regenerate with `snakemake walkthrough`.*")
    o.append("")
    o.append(f"- Config the TUNABLES were resolved from: `{cfgp}`")
    o.append(f"- Commit: `{commit_hash}`")
    o.append("")
    o.append("One page per stage, in run order. Each page says what the stage "
             "does, what it needs from the stage before it, the exact command to "
             "run up to it, and — per rule — what that rule does, which config "
             "keys steer it, what it writes, and how to re-run just that rule.")
    o.append("")
    o.append("For the same material as one continuous read, see "
             "[`../WALKTHROUGH.md`](../WALKTHROUGH.md). Both are generated from "
             "the same rule docstrings in the same pass, so they cannot disagree.")
    o.append("")
    o.append("## Stepping through manually")
    o.append("")
    o.append("The pipeline does not have to be run in one go. Snakemake will stop "
             "where you tell it to, and pick up from there once you have looked at "
             "the output — which is the point of these pages: run a stage, check "
             "the figure or the table it produced, change a threshold or drop a "
             "bad sample, re-run that stage, carry on.")
    o.append("")
    o.append("Four things do almost all of it:")
    o.append("")
    o.append("| Invocation | What it does |")
    o.append("|---|---|")
    o.append("| `snakemake <output files>` | Build exactly those files and whatever they need, then stop. This is \"run up to step N\", and it is what each page prints. |")
    o.append("| `snakemake --force <output file>` | Re-run the one job that makes that file, and only it. This is \"redo this step\". |")
    o.append("| `snakemake --forcerun <rule>` | Re-run that rule *and everything downstream of it*. This is \"I changed a setting, so redo this and everything computed from it\". |")
    o.append("| `-n` (`--dry-run`) | Print what would run without running it. Worth doing before any `--forcerun`. |")
    o.append("")
    o.append("Three things worth knowing before you rely on this:")
    o.append("")
    o.append("- **Editing the config does not by itself trigger a re-run.** "
             "Snakemake decides what is out of date mostly from file "
             "timestamps, so after changing a threshold you have to force the "
             "rules it affects. That is what the force commands on each page "
             "are for.")
    o.append("- **`--forcerun` re-runs what comes after, too.** That is usually "
             "what you want — a changed threshold invalidates every figure "
             "drawn from it — but on a large cohort it is worth a `-n` first "
             "to see the size of what you just asked for. When you truly want "
             "one job and nothing else, use `--force <file>`.")
    o.append("- **`--until` is the obvious flag here, and it does not work "
             "past Stage 3.** Stage 3 ends in the `assign_clusters` "
             "checkpoint, which decides how many clusters exist and therefore "
             "which downstream jobs exist at all. `--until` prunes the graph "
             "before that is known and deadlocks with \"Out of jobs ready to be "
             "started, but not all files built yet\" "
             "([snakemake#823](https://github.com/snakemake/snakemake/issues/823)). "
             "Naming the output files avoids the problem entirely, which is "
             "why these pages do that.")
    o.append("")
    o.append("To edit an intermediate file by hand and keep it, `touch` it "
             "afterwards so it is newer than its inputs — otherwise the next run "
             "regenerates it and your edit is gone.")
    o.append("")
    o.append("## The stages")
    o.append("")
    o.append("| # | Stage | Purpose |")
    o.append("|---|---|---|")
    for i, st in enumerate(STAGES, start=1):
        o.append(f"| {i} | [{st.title}]({st.slug}.md) | {st.purpose} |")
    o.append("")
    o.append("## Rules per stage")
    o.append("")
    for st in STAGES:
        rules = stage_rules(st, rules_by_file)
        names = ", ".join(f"`{r.name}`" for r in rules) if rules else "*(none)*"
        o.append(f"- **[{st.title}]({st.slug}.md)** — {names}")
    o.append("")
    return "\n".join(o).rstrip() + "\n"


def render_steps(cfg: dict, config_path: Path, commit_hash: str,
                 rules_by_file: Dict[str, List[Rule]]) -> Dict[str, str]:
    """Every file under docs/steps/, as {relative filename: content}."""
    pages: Dict[str, str] = {}
    total = len(STAGES)
    for i, st in enumerate(STAGES):
        pages[f"{st.slug}.md"] = render_step_page(
            st, cfg, config_path, commit_hash, rules_by_file,
            position=i + 1, total=total,
            prev_stage=STAGES[i - 1] if i > 0 else None,
            next_stage=STAGES[i + 1] if i + 1 < total else None,
        )
    pages["README.md"] = render_steps_index(cfg, config_path, commit_hash,
                                            rules_by_file)
    return pages


def _strip_stamp(t: str) -> str:
    """
    Drop the header stamp lines before diffing. They record which config and
    which commit a render came from, so they differ legitimately between a
    committed copy and a local regeneration.
    """
    return "\n".join(
        ln for ln in t.splitlines()
        if not ln.startswith("- Config:")
        and not ln.startswith("- Commit:")
        and not ln.startswith("- Step ")
    )


def _orphan_step_files(generated: set) -> List[Path]:
    """Markdown files under docs/steps/ that this run did not generate."""
    if not STEPS_DIR.is_dir():
        return []
    return sorted(p for p in STEPS_DIR.glob("*.md") if p.name not in generated)


def _check(expected: Dict[Path, str]) -> int:
    """
    Compare every generated document against its committed copy. Reports all
    stale files, not just the first, so one regeneration fixes the lot.
    """
    stale: List[str] = []
    for path, content in expected.items():
        rel = path.relative_to(AGNOSTIC)
        if not path.exists():
            stale.append(f"{rel}: missing")
            continue
        if _strip_stamp(path.read_text()) == _strip_stamp(content):
            continue
        stale.append(f"{rel}: differs")
        diff = list(difflib.unified_diff(
            _strip_stamp(path.read_text()).splitlines(keepends=True),
            _strip_stamp(content).splitlines(keepends=True),
            fromfile=f"committed/{rel}",
            tofile=f"regenerated/{rel}",
            n=2,
        ))
        sys.stderr.write("".join(diff[:60]))
        if len(diff) > 60:
            sys.stderr.write(f"... ({len(diff) - 60} more diff lines suppressed) ...\n")

    orphans = _orphan_step_files({p.name for p in expected if p.parent == STEPS_DIR})
    for o in orphans:
        stale.append(f"{o.relative_to(AGNOSTIC)}: orphan (no stage generates it)")

    if not stale:
        print(f"Docs up to date ({len(expected)} files checked).")
        return 0
    print("", file=sys.stderr)
    for entry in stale:
        print(f"STALE  {entry}", file=sys.stderr)
    print("Regenerate with `snakemake walkthrough`.", file=sys.stderr)
    return 1


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
    p.add_argument("--write",  action="store_true", help="Regenerate WALKTHROUGH.md + docs/steps/ in place.")
    p.add_argument("--check",  action="store_true", help="Fail if any committed doc is stale.")
    args = p.parse_args(argv)

    if not args.write and not args.check:
        args.write = True  # default

    cfg = load_config(args.config)
    commit = _git_commit()

    rules_by_file: Dict[str, List[Rule]] = {}
    for smk in sorted(RULES_DIR.glob("*.smk")):
        rules_by_file[str(smk.relative_to(AGNOSTIC))] = parse_rules_from_smk(smk)

    # One pass, two products: the single-file linear read and the per-step
    # folder. Generating them together is what stops them disagreeing.
    rendered = render_walkthrough(cfg, args.config, commit, rules_by_file)
    steps    = render_steps(cfg, args.config, commit, rules_by_file)

    # {path relative to AGNOSTIC: expected content} — everything under the
    # drift check, in the order a human would want to see failures reported.
    expected: Dict[Path, str] = {WALKTHROUGH_PATH: rendered}
    for fname, content in steps.items():
        expected[STEPS_DIR / fname] = content

    if args.check:
        return _check(expected)

    WALKTHROUGH_PATH.parent.mkdir(parents=True, exist_ok=True)
    STEPS_DIR.mkdir(parents=True, exist_ok=True)
    for path, content in expected.items():
        path.write_text(content)
    # Drop step pages left behind by a stage that no longer exists, so the
    # folder is exactly what this run generated and --check cannot pass on a
    # directory holding an orphan.
    for stale in _orphan_step_files(set(steps)):
        stale.unlink()
        print(f"Removed stale {stale.relative_to(AGNOSTIC)}.")
    total = sum(len(c) for c in expected.values())
    print(f"Wrote {WALKTHROUGH_PATH.relative_to(AGNOSTIC)} + "
          f"{len(steps)} files under {STEPS_DIR.relative_to(AGNOSTIC)}/ "
          f"({total} bytes).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
