#!/usr/bin/env bash
# env/postinstall.sh — the part of the environment conda cannot supply.
#
#   pixi run setup
#
# `pixi install` gets you everything in pixi.lock. Four things are not in there,
# because no conda channel publishes them for every platform we support:
#
#   moimix               GitHub only (bahlolab/moimix) — no conda build anywhere
#   rehh                 CRAN only  — no conda build anywhere
#   rnaturalearthhires   r-universe only — no conda build anywhere
#   SeqArray             bioconda linux-64 only; installed via BiocManager on macOS
#   SeqVarTools          bioconda linux-64 only; installed via BiocManager on macOS
#   ADMIXTURE, PLINK2    bioconda linux-64 only; upstream binaries on macOS
#
# Everything here is pinned and idempotent: re-running is a series of skips.
# Run it from the workspace root (pixi does this for you).

set -eo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
cd "$ROOT"

PLATFORM="$(uname -s)-$(uname -m)"
echo "==> ZOOMAL-Flow postinstall on ${PLATFORM}"

# Where pixi put the solved env. `pixi run` exports CONDA_PREFIX; fall back to
# the conventional path so the script also works when sourced by hand.
ENV_PREFIX="${CONDA_PREFIX:-$ROOT/.pixi/envs/default}"
ENV_BIN="$ENV_PREFIX/bin"
if [[ ! -x "$ENV_BIN/Rscript" ]]; then
  echo "ERROR: no Rscript at $ENV_BIN. Run 'pixi install' first, then 'pixi run setup'." >&2
  exit 1
fi
echo "    env prefix: $ENV_PREFIX"

# Pins. Bump deliberately and record the bump in env/README.md.
REHH_VERSION="3.2.3"
ADMIXTURE_VERSION="1.3.0"
PLINK2_BUILD="plink2_mac_arm64_20250707"          # matches V1's PLINK v2.0.0-a.6.20
CRAN="https://cloud.r-project.org"

#---------------------------------------------------------------------------
# R packages that no conda channel carries.
#---------------------------------------------------------------------------
echo "==> R packages not available from conda"
"$ENV_BIN/Rscript" --vanilla - "$REHH_VERSION" "$CRAN" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
rehh_version <- args[[1]]; cran <- args[[2]]

have <- function(p) requireNamespace(p, quietly = TRUE)

# Bioconductor. On linux-64 conda already provided SeqArray; SeqVarTools and
# BiocParallel are moimix's build-time dependencies and are needed everywhere.
bioc <- c("SeqArray", "SeqVarTools", "BiocParallel")
for (p in bioc) {
  if (have(p)) message("    [skip] ", p, " ", as.character(packageVersion(p)))
}
need <- bioc[!vapply(bioc, have, logical(1))]
if (length(need)) {
  message("    [add ] ", paste(need, collapse = ", "), " (BiocManager)")
  BiocManager::install(need, ask = FALSE, update = FALSE)
}

# rehh — CRAN only. Pinned: the iHS/EHH statistics are version-sensitive.
if (have("rehh") && as.character(packageVersion("rehh")) == rehh_version) {
  message("    [skip] rehh ", rehh_version)
} else {
  message("    [add ] rehh ", rehh_version, " (CRAN archive, pinned)")
  remotes::install_version("rehh", version = rehh_version,
                           repos = cran, upgrade = "never")
}

# rnaturalearthhires — r-universe only. rnaturalearth::ne_states() needs it;
# without it the province map silently falls back and then errors.
if (have("rnaturalearthhires")) {
  message("    [skip] rnaturalearthhires ", as.character(packageVersion("rnaturalearthhires")))
} else {
  message("    [add ] rnaturalearthhires (r-universe)")
  install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")
}

# moimix — GitHub only. This is what computes Fws in Stage 2.
if (have("moimix")) {
  message("    [skip] moimix ", as.character(packageVersion("moimix")))
} else {
  message("    [add ] moimix (github:bahlolab/moimix)")
  remotes::install_github("bahlolab/moimix", upgrade = "never")
}
RSCRIPT

#---------------------------------------------------------------------------
# Binaries with no osx-arm64 conda build. On Linux, pixi already installed
# both from bioconda and these blocks are pure skips.
#---------------------------------------------------------------------------
echo "==> Command-line tools not available from conda on this platform"

# ADMIXTURE — upstream ships x86-64 macOS only; it runs under Rosetta 2.
# (Same arrangement V1 uses. If Rosetta is absent: softwareupdate --install-rosetta)
if [[ -x "$ENV_BIN/admixture" ]]; then
  echo "    [skip] admixture already in $ENV_BIN"
else
  echo "    [add ] ADMIXTURE ${ADMIXTURE_VERSION} (x86-64 macOS binary, runs under Rosetta)"
  TMP="$(mktemp -d)"
  curl -fsSL "https://dalexander.github.io/admixture/binaries/admixture_macosx-${ADMIXTURE_VERSION}.tar.gz" \
    -o "$TMP/admixture.tar.gz"
  tar -xzf "$TMP/admixture.tar.gz" -C "$TMP"
  BIN="$(find "$TMP" -type f -name admixture | head -1)"
  [[ -n "$BIN" ]] || { echo "ERROR: no admixture binary in the downloaded archive." >&2; exit 1; }
  install -m 755 "$BIN" "$ENV_BIN/admixture"
  rm -rf "$TMP"
fi

# PLINK2 — bioconda has no osx-arm64 build; upstream ships a native arm64 one.
if [[ -x "$ENV_BIN/plink2" ]]; then
  echo "    [skip] plink2 already in $ENV_BIN"
else
  echo "    [add ] PLINK2 (${PLINK2_BUILD})"
  TMP="$(mktemp -d)"
  curl -fsSL "https://s3.amazonaws.com/plink2-assets/alpha6/${PLINK2_BUILD}.zip" \
    -o "$TMP/plink2.zip"
  ( cd "$TMP" && unzip -q plink2.zip )
  BIN="$(find "$TMP" -type f -name plink2 | head -1)"
  [[ -n "$BIN" ]] || { echo "ERROR: no plink2 binary in the downloaded archive." >&2; exit 1; }
  install -m 755 "$BIN" "$ENV_BIN/plink2"
  rm -rf "$TMP"
fi

#---------------------------------------------------------------------------
# Verify, and fail loudly if anything did not land.
#
# install.packages() and BiocManager::install() report a failed install as a
# WARNING, not an error: R exits 0 even when the package is not there. That is
# exactly how a clean-machine moimix failure once got past this script and
# printed "postinstall complete" — the student only found out several minutes
# into a run, as a cryptic Stage 2 error. "Complete" must mean usable, so the
# result is checked here rather than left to whether anyone runs check-env.
#---------------------------------------------------------------------------
echo
echo "==> Verifying"

FAILED=""

MISSING_R="$("$ENV_BIN/Rscript" --vanilla - "$REHH_VERSION" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
rehh_version <- args[[1]]

# Everything this script is responsible for putting on disk.
want <- c("SeqArray", "SeqVarTools", "BiocParallel",
          "rnaturalearthhires", "rehh", "moimix")

bad <- character(0)
for (p in want) {
  if (!requireNamespace(p, quietly = TRUE)) {
    bad <- c(bad, p)
  } else if (p == "rehh" && as.character(packageVersion(p)) != rehh_version) {
    bad <- c(bad, sprintf("rehh(have-%s-want-%s)", packageVersion(p), rehh_version))
  } else {
    message("    [ok  ] ", p, " ", as.character(packageVersion(p)))
  }
}
cat(paste(bad, collapse = " "))
RSCRIPT
)"

# Plain `if`, not `[[ ... ]] && ...`: under `set -e` a trailing test that
# evaluates false would exit the script on the success path.
if [[ -n "$MISSING_R" ]]; then
  FAILED="$FAILED $MISSING_R"
fi

for tool in admixture plink2; do
  if [[ -x "$ENV_BIN/$tool" ]]; then
    echo "    [ok  ] $tool"
  else
    FAILED="$FAILED $tool"
  fi
done

if [[ -n "$FAILED" ]]; then
  cat >&2 <<EOF

ERROR: postinstall did NOT complete. Missing:$FAILED

The environment is not usable yet — do not run the pipeline against it.

If moimix or SeqVarTools is in that list, the usual cause is that the
dependency chain below it (logistf -> mice -> mitml -> jomo -> lme4 -> nloptr)
fell back to building from source and failed. Those are pinned as conda
binaries in pixi.toml precisely so that cannot happen, so first check the
locked environment is actually the one in use:

  pixi install          # re-materialise the env from pixi.lock
  pixi run setup        # then re-run this script

Scroll up for the R error that caused it — it names the first package that
failed, which is the one to chase.
EOF
  exit 1
fi

echo
echo "==> postinstall complete — every dependency verified present."
echo "    Full environment check (conda stack included): pixi run check-env"
