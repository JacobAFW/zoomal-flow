#!/usr/bin/env python3
"""
Unit tests for scripts/py/render_walkthrough.py — the docstring parser and
the config resolver. Runs with or without pytest.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
AGNOSTIC = HERE.parent
sys.path.insert(0, str(AGNOSTIC / "scripts" / "py"))

from render_walkthrough import (  # noqa: E402
    parse_docstring, parse_rules_from_smk, resolve_dotted, load_config,
    STAGES, RULES_DIR, stage_rules, terminal_rules, active_rules,
    stage_targets, resolve_output_template, render_steps, Rule,
)


def _rules_by_file():
    return {str(p.relative_to(AGNOSTIC)): parse_rules_from_smk(p)
            for p in sorted(RULES_DIR.glob("*.smk"))}


def test_parse_docstring_basic():
    body = """
    One-line summary.

    WHAT: bcftools view -m2 -M2
    WHY:  biallelic downstream
    TUNABLES: qc.min_mac, qc.filter_pass
    OUTPUT: {outputs}/qc/snps.qc.vcf.gz
    TRY:    bump qc.min_mac to 5
    """
    sec = parse_docstring(body)
    assert sec["summary"] == "One-line summary."
    assert sec["what"].startswith("bcftools")
    assert sec["why"].startswith("biallelic")
    assert "qc.min_mac" in sec["tunables"]
    assert "snps.qc.vcf.gz" in sec["output"]
    assert "bump qc.min_mac" in sec["try"]


def test_parse_docstring_multiline_why():
    body = """
    Summary.

    WHAT: x
    WHY:  first sentence.
          second sentence still in WHY.
    OUTPUT: out
    """
    sec = parse_docstring(body)
    assert "first sentence" in sec["why"]
    assert "second sentence" in sec["why"]


def test_parse_rules_from_actual_setup_smk():
    smk = AGNOSTIC / "workflow" / "rules" / "00_setup.smk"
    rules = parse_rules_from_smk(smk)
    names = {r.name for r in rules}
    assert "index_vcf" in names
    assert "extract_vcf_samples" in names
    assert "validate_metadata" in names
    r = next(r for r in rules if r.name == "extract_vcf_samples")
    assert "bcftools query -l" in r.what


def test_resolve_dotted():
    cfg = {"qc": {"min_mac": 2}, "cohort": {"input_type": "wgs"}}
    assert resolve_dotted(cfg, "qc.min_mac") == 2
    assert resolve_dotted(cfg, "cohort.input_type") == "wgs"
    assert resolve_dotted(cfg, "does.not.exist") is None


def test_load_config_layers():
    # config.yaml + cohort.example.yaml merge
    cfg = load_config(AGNOSTIC / "config" / "cohort.example.yaml")
    assert cfg["cohort"]["input_type"] == "wgs"
    assert cfg["structure"]["cluster_labelling"] == "reference"


# --------------------------------------------------------------------------
# Per-step pages (docs/steps/)
# --------------------------------------------------------------------------

def test_every_rule_belongs_to_exactly_one_stage():
    """
    A rule missing from STAGES would silently vanish from docs/steps/ while
    still appearing in WALKTHROUGH.md — the exact drift these pages exist to
    prevent. A rule in two stages would be documented twice.
    """
    rbf = _rules_by_file()
    everything = [r.name for rs in rbf.values() for r in rs]
    placed = []
    for st in STAGES:
        placed.extend(r.name for r in stage_rules(st, rbf))
    assert sorted(placed) == sorted(everything), (
        f"unplaced: {sorted(set(everything) - set(placed))}; "
        f"duplicated: {sorted(n for n in set(placed) if placed.count(n) > 1)}"
    )


def test_terminal_rules_drops_rules_consumed_in_stage():
    a = Rule(name="a", file="x.smk", output_paths=["{P}/one.txt"])
    b = Rule(name="b", file="x.smk", input_paths=["{P}/one.txt"],
             output_paths=["{P}/two.txt"])
    c = Rule(name="c", file="x.smk", dep_rules=["b"], output_paths=["{P}/three.txt"])
    term = [r.name for r in terminal_rules([a, b, c])]
    assert term == ["c"], term


def test_terminal_rules_handles_two_rules_declaring_one_path():
    """The WGS and microhap seams declare the same outputs; both must count as
    producers, or whichever parsed first is wrongly reported terminal."""
    wgs   = Rule(name="wgs_prep", file="p_wgs.smk", output_paths=["{P}/ld.bed"])
    micro = Rule(name="mh_prep", file="p_microhap.smk", output_paths=["{P}/ld.bed"])
    user  = Rule(name="user", file="x.smk", input_paths=["{P}/ld.bed"],
                 output_paths=["{P}/out.txt"])
    term = [r.name for r in terminal_rules([wgs, micro, user])]
    assert term == ["user"], term


def test_active_rules_follows_input_type():
    wgs   = Rule(name="w", file="workflow/rules/qc_wgs.smk")
    micro = Rule(name="m", file="workflow/rules/qc_microhap.smk")
    shared = Rule(name="s", file="workflow/rules/01_qc.smk")
    got = [r.name for r in active_rules([wgs, micro, shared],
                                        {"cohort": {"input_type": "wgs"}})]
    assert got == ["w", "s"], got


def test_resolve_output_template_resolves_and_rejects_wildcards():
    cfg = {"paths": {"outputs": "outputs"}, "selection": {"models": []}}
    assert resolve_output_template("{PATHS['outputs']}/qc/a.vcf.gz", cfg) == "outputs/qc/a.vcf.gz"
    # {sampleset} is known from config; {cluster} is checkpoint-derived and is not
    assert resolve_output_template("{PATHS['outputs']}/s/{{sampleset}}/b.png", cfg) == "outputs/s/full/b.png"
    assert resolve_output_template("{PATHS['outputs']}/ibd/{{cluster}}/c.txt", cfg) is None
    # a module-level constant from the real rule files resolves too
    assert resolve_output_template("{CLONALITY_DIR}/all_genotypes.txt", cfg) == "outputs/clonality/all_genotypes.txt"


def test_every_stage_offers_a_way_to_run_it():
    """
    Each page must name either concrete target files or an explicit named
    target. A stage with neither would print a --until command, and --until
    cannot cross the assign_clusters checkpoint.
    """
    cfg = load_config(AGNOSTIC / "config" / "cohort.example.yaml")
    rbf = _rules_by_file()
    for st in STAGES:
        runnable = active_rules(stage_rules(st, rbf), cfg)
        assert st.target or stage_targets(runnable, cfg), f"{st.slug} has no target"


def test_render_steps_emits_a_page_per_stage_plus_index():
    cfg = load_config(AGNOSTIC / "config" / "cohort.example.yaml")
    pages = render_steps(cfg, AGNOSTIC / "config" / "cohort.example.yaml", "abc123",
                         _rules_by_file())
    assert set(pages) == {f"{st.slug}.md" for st in STAGES} | {"README.md"}
    # every rule is documented on its stage's page
    rbf = _rules_by_file()
    for st in STAGES:
        page = pages[f"{st.slug}.md"]
        for r in stage_rules(st, rbf):
            assert f"### `{r.name}`" in page, f"{r.name} missing from {st.slug}.md"


def _main():
    failures = 0
    tests = [
        test_parse_docstring_basic,
        test_parse_docstring_multiline_why,
        test_parse_rules_from_actual_setup_smk,
        test_resolve_dotted,
        test_load_config_layers,
        test_every_rule_belongs_to_exactly_one_stage,
        test_terminal_rules_drops_rules_consumed_in_stage,
        test_terminal_rules_handles_two_rules_declaring_one_path,
        test_active_rules_follows_input_type,
        test_resolve_output_template_resolves_and_rejects_wildcards,
        test_every_stage_offers_a_way_to_run_it,
        test_render_steps_emits_a_page_per_stage_plus_index,
    ]
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
        except AssertionError as e:
            failures += 1
            print(f"  FAIL  {t.__name__}: {e}")
    if failures:
        print(f"\n{failures} test(s) failed", file=sys.stderr)
        sys.exit(1)
    print("\nAll tests passed.")


if __name__ == "__main__":
    _main()
