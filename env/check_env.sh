#!/usr/bin/env bash
# env/check_env.sh — does this environment actually have everything the pipeline
# invokes? Prints a version for each tool and R package, and exits non-zero if
# any are missing.
#
#   pixi run check-env
#
# Keep the two lists below in sync with what the rules and scripts import: the
# point of this script is that a missing dependency surfaces here, in seconds,
# rather than three hours into a run.

set -o pipefail

MISSING=0

echo "=== command-line tools ==="
for tool in bcftools bgzip tabix samtools plink plink2 admixture hmmIBD Rscript snakemake quarto python; do
  if command -v "$tool" >/dev/null 2>&1; then
    # hmmIBD prints usage to stderr and has no --version flag.
    if [[ "$tool" == "hmmIBD" ]]; then
      ver="(no --version; responds to invocation)"
    else
      ver="$("$tool" --version 2>&1 | head -1)"
    fi
    printf "  %-11s %s\n" "$tool" "$ver"
  else
    printf "  %-11s MISSING\n" "$tool"
    MISSING=1
  fi
done

echo
echo "=== R packages ==="
Rscript --vanilla -e '
pkgs <- c(
  # core
  "tidyverse", "data.table", "jsonlite", "R.utils",
  # plotting
  "ggplot2", "viridis", "viridisLite", "ggnewscale", "ggbreak", "svglite",
  # Stage 2 MOI
  "SeqArray", "moimix", "quantreg", "mgcv",
  # Stage 3 structure / trees / maps
  "ape", "ggtree", "sf", "sp", "rnaturalearth", "rnaturalearthdata",
  "rnaturalearthhires",
  # Stage 4 IBD networks
  "igraph",
  # Stage 6 selection
  "rehh"
)
bad <- character(0)
for (p in pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
  if (is.na(v)) { bad <- c(bad, p); cat(sprintf("  %-20s MISSING\n", p)) }
  else cat(sprintf("  %-20s %s\n", p, v))
}
if (length(bad)) quit(status = 1)
' || MISSING=1

echo
echo "=== python modules ==="
python - <<'PY' || MISSING=1
import importlib, sys
mods = ["yaml", "jsonschema", "numpy", "pytest", "snakemake"]
bad = []
for m in mods:
    try:
        mod = importlib.import_module(m)
        print(f"  {m:<20} {getattr(mod, '__version__', 'ok')}")
    except ImportError:
        print(f"  {m:<20} MISSING"); bad.append(m)
sys.exit(1 if bad else 0)
PY

echo
if [[ "$MISSING" -ne 0 ]]; then
  echo "FAIL — something above is missing. Run 'pixi install' then 'pixi run setup'."
  exit 1
fi
echo "OK — every pipeline dependency resolves."
