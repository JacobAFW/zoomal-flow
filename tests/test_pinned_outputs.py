"""Byte-identity guard for the stages NOT under active change.

Stage 6 is being replaced. Nothing in the test suite previously pinned any
pipeline output — and `tests/tiny_cohort/config.yaml` ships
`selection.models: []`, so Stage 6 never runs in CI either. That means a
Stage-6 change could perturb an upstream stage and nothing would catch it,
in either direction.

These tests pin the byte-stable outputs of Stages 0-5 plus the clonality pass
on the tiny cohort. The manifest is `tests/tiny_cohort/pinned_outputs.sha256`;
what is pinned and what deliberately is not is recorded in
`CONTEXT-pipeline.md`.

The guard SKIPS when the tiny cohort has not been run, so a fresh clone is not
failed for something it cannot know. It FAILS when the outputs exist and
disagree with the manifest. Point it elsewhere with

    TINY_COHORT_OUTPUTS=/path/to/outputs pytest tests/test_pinned_outputs.py

which is how a run from a scratch directory can be checked without touching
the committed outputs.
"""

import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "tiny_cohort"))

import pin_outputs  # noqa: E402


def _outputs_dir():
    return os.environ.get("TINY_COHORT_OUTPUTS", pin_outputs.DEFAULT_OUTPUTS)


@pytest.fixture(scope="module")
def pinned():
    manifest = pin_outputs.read_manifest()
    if manifest is None:
        pytest.fail(
            "tests/tiny_cohort/pinned_outputs.sha256 is missing. Regenerate "
            "it with: python3 tests/tiny_cohort/pin_outputs.py --write")
    return manifest


@pytest.fixture(scope="module")
def produced():
    out = _outputs_dir()
    if not os.path.isdir(out):
        pytest.skip(
            f"{out} does not exist — run the tiny cohort first:\n"
            f"  pixi run snakemake --configfile tests/tiny_cohort/config.yaml "
            f"--cores 4")
    found = pin_outputs.walk(out)
    if not found:
        pytest.skip(f"{out} has no pinnable files — the tiny cohort has not "
                    f"been run")
    return found


def test_manifest_is_not_empty(pinned):
    """A manifest that pins nothing is not a guard."""
    assert len(pinned) > 50, (
        f"only {len(pinned)} files pinned; the tiny cohort produces ~129 "
        f"byte-stable outputs, so the manifest looks truncated")


def test_no_unstable_files_were_pinned(pinned):
    """Pinning a log or a bgzipped VCF would make the guard flap."""
    bad = [rel for rel in pinned if not pin_outputs.is_pinnable(rel)]
    assert bad == [], (
        f"these pinned paths are known not to be byte-stable between runs: "
        f"{bad}")


def test_every_pinned_file_still_exists(pinned, produced):
    missing = sorted(set(pinned) - set(produced))
    assert missing == [], (
        f"{len(missing)} pinned output(s) are no longer produced: "
        f"{missing[:10]}. If an upstream stage was changed deliberately, "
        f"re-run the tiny cohort and regenerate the manifest, and record the "
        f"change in CONTEXT-pipeline.md.")


def test_no_pinned_file_changed(pinned, produced):
    changed = []
    for rel, digest in sorted(pinned.items()):
        if rel not in produced:
            continue                      # reported by the test above
        if pin_outputs.sha256(produced[rel]) != digest:
            changed.append(rel)
    assert changed == [], (
        f"{len(changed)} upstream output(s) changed: {changed[:10]}. These "
        f"stages were not supposed to move. Either the change is a "
        f"regression, or it is deliberate and the manifest needs "
        f"regenerating with a note in CONTEXT-pipeline.md.")


def test_no_unpinned_outputs_appeared(pinned, produced):
    """A new byte-stable output is fine, but it should be pinned on purpose.

    Requires a CLEAN output directory. Snakemake does not delete outputs that
    a later version of a rule stopped producing, so a long-lived outputs tree
    accumulates leftovers from older code — the committed
    `tests/tiny_cohort/outputs` carries six of them from before commit
    45c08fd changed the ADMIXTURE working-dir layout. That is what this test
    is most likely telling you about.
    """
    extra = sorted(set(produced) - set(pinned))
    assert extra == [], (
        f"{len(extra)} byte-stable output(s) are not in the manifest: "
        f"{extra[:10]}.\n"
        f"Most likely these are stale files from an older version of a rule — "
        f"Snakemake never deletes them. Re-run into a clean directory to "
        f"confirm. If they are genuinely new and intended, add them with "
        f"pin_outputs.py --write and note it in CONTEXT-pipeline.md.")
