# Makefile — thin entry points for the things that are NOT pipeline rules.
#
# The pipeline itself is Snakemake (`snakemake --configfile config/...`). What
# lands here is method-development work that cannot be a rule because it needs
# inputs the pipeline does not produce: a fixture with published ground truth.
#
# Targets:
#   make benchmark-malay   score every detection rule against the published
#                          Malaysian introgression result (Stage 5b)
#   make test              the R unit tests + the pytest suite

AGNOSTIC     := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
RSCRIPT      ?= Rscript
PYTEST       ?= pytest

BENCH_INPUTS := data/benchmark_malay/inputs
BENCH_OUT    ?= outputs/benchmark_malay
BENCH_LOG    := logs/benchmark/malay_benchmark.log

.PHONY: benchmark-malay test

## Score absolute / relative / distance / distance_adaptive against the Malay
## fixture. REPORTS ONLY — it changes no default.
##
## The fixture is real embargoed data under data/benchmark_malay/ and its
## outputs land under outputs/benchmark_malay/. BOTH are gitignored and must
## never be staged. Override the destination with BENCH_OUT=... if you want it
## somewhere else; keep it outside version control either way.
benchmark-malay:
	@test -f "$(BENCH_INPUTS)/hmmIBD.tsv" || { \
	  echo "ERROR: benchmark fixture missing at $(BENCH_INPUTS)/."; \
	  echo "       See data/benchmark_malay/README.md — the fixture is not in git."; \
	  exit 1; }
	@mkdir -p "$(dir $(BENCH_LOG))"
	$(RSCRIPT) scripts/R/introgression_benchmark_malay.R \
	    --inputs-dir "$(BENCH_INPUTS)" \
	    --truth-dir  data/benchmark_malay/truth \
	    --regen-dir  data/benchmark_malay/regen_density \
	    --out-dir    "$(BENCH_OUT)" \
	  2>&1 | tee "$(BENCH_LOG)"
	@echo
	@echo "scoreboard: $(BENCH_OUT)/scoreboard.tsv   (log: $(BENCH_LOG))"

test:
	$(RSCRIPT) tests/R/test_introgression_units.R
	$(PYTEST) tests -q
