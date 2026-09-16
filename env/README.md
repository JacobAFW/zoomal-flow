# The ZOOMAL-Flow environment

ZOOMAL-Flow owns its environment. Clone the repo, run two commands, and every
tool the pipeline invokes is present at a recorded version — no reference to
the V1 Indonesia pipeline's `vvg-box` env, which stays on its own separate
install and is untouched by anything here.

```bash
cd agnostic
pixi install        # the conda-solvable stack, exactly as recorded in pixi.lock
pixi run setup      # the six things conda cannot supply (below)
pixi run check-env  # prints a version for every tool + package, fails if any is missing
```

Then run the pipeline through pixi so it picks up the environment:

```bash
pixi run snakemake --configfile config/config.yaml --cores 8
```

or open a shell in it once and work normally:

```bash
pixi shell
snakemake --configfile config/config.yaml --cores 8
```

## Where the files live

| File | Tracked? | What it is |
|---|---|---|
| `../pixi.toml` | yes | the dependency declaration — what we ask for |
| `../pixi.lock` | yes | the exact solve for `osx-arm64` + `linux-64` — what you get |
| `postinstall.sh` | yes | installs the tools no conda channel publishes |
| `check_env.sh` | yes | verifies the result; run it after any env change |
| `environment.yml` | yes | conda/mamba fallback for sites without pixi (unpinned solve) |
| `../.pixi/` | **no** | the ~2.2 GB solved environment itself; gitignored |

`pixi.toml` sits in the repo root rather than in this directory so that
`cd agnostic && pixi install` works verbatim — pixi discovers its manifest by
walking up from the working directory, and a manifest inside `env/` would force
every command to carry `--manifest-path`. This directory holds everything else.

## What conda cannot give us

Six dependencies are not in `pixi.lock`, for reasons outside our control.
`postinstall.sh` handles all of them, pinned and idempotent:

| Dependency | Why not conda | Where it comes from |
|---|---|---|
| `moimix` | published on GitHub only, never packaged for conda | `remotes::install_github("bahlolab/moimix")` |
| `rehh` | CRAN only; no conda build on any platform | `remotes::install_version("rehh", "3.2.3")` — pinned, iHS statistics are version-sensitive |
| `rnaturalearthhires` | r-universe only (too large for CRAN) | `install.packages(repos = "https://ropensci.r-universe.dev")` |
| `SeqArray` | bioconda builds it for linux-64 but not osx-arm64 | conda on Linux; `BiocManager::install()` on macOS |
| `ADMIXTURE` | upstream ships x86-64 macOS only | bioconda on Linux; upstream x86-64 binary under Rosetta 2 on macOS |
| `PLINK2` | bioconda has no osx-arm64 build | bioconda on Linux; upstream `plink2_mac_arm64_20250707` on macOS |

`SeqVarTools` and `BiocParallel` come along too — they are moimix's build-time
dependencies, and moimix will not install without them.

On **linux-64** the last three rows are already in the lock file, so
`postinstall.sh` is nearly a no-op there: it installs the three R packages and
skips the binaries.

On **macOS**, ADMIXTURE runs through Rosetta 2. If it is not installed:

```bash
softwareupdate --install-rosetta
```

## What is pinned

Versions verified present after `pixi install && pixi run setup` on
Darwin-arm64 (M2 Max), 2026-09-16. The right-hand column is what the V1
pipeline ran with — the results this pipeline was developed against.

| Tool / package | ZOOMAL-Flow | V1 vvg-box | |
|---|---|---|---|
| R | 4.5.3 | 4.5.3 | match |
| bcftools / htslib / samtools | 1.23.1 | 1.23.1 | match |
| PLINK 1.9 | v1.9.0-b.8 | v1.9.0-b.8 | match |
| PLINK 2 | v2.0.0-a.6.20 M1 | v2.0.0-a.6.20 M1 | match |
| ADMIXTURE | 1.3.0 | 1.3.0 | match |
| snakemake | 9.22.0 | 9.22.0 | match |
| quarto | 1.9.38 | 1.9.38 | match |
| moimix | 0.0.2.9001 | 0.0.2.9001 | match |
| SeqArray | 1.50.1 | 1.50.1 | match |
| rehh | 3.2.3 | 3.2.3 | match |
| ggtree | 4.0.5 | 4.0.5 | match |
| tidyverse / data.table | 2.0.0 / 1.17.8 | 2.0.0 / 1.17.8 | match |
| quantreg / mgcv | 6.1 / 1.9.4 | 6.1 / 1.9.4 | match |
| ape / igraph | 5.8.1 / 2.3.3 | 5.8.1 / 2.3.2 | patch drift |
| sf / sp | 1.1.3 / 2.2.3 | 1.1.0 / 2.2.1 | patch drift |
| python | 3.12.14 | 3.12.13 | patch drift |
| hmmIBD | 2.1.3 (bioconda) | built from source, unversioned | **see below** |

### hmmIBD is the one real deviation

V1 compiles hmmIBD from a vendored copy of the source with no version string in
it. That is not something a collaborator can reproduce from a manifest, so
ZOOMAL-Flow takes bioconda's `hmmibd=2.1.3` instead. Both are hmmIBD 2.x and
expose the same interface (`-i`, `-o`, `-I`, `-f`, `-b`, `-g`, `-m`, `-n`,
`-r`), and the pipeline only uses `-i` / `-o`.

This was checked, not assumed: a full tiny-cohort run under this environment
produces IBD output (`*.hmm.txt`, `*.hmm_fract.txt`, the derived
`clonal_clusters.tsv`) byte-identical to the V1-env run. See "Verification".

## Verification

The environment is verified by reproducing a known output, not by inspection.
On 2026-09-16 the whole tiny synthetic cohort (`tests/tiny_cohort/`, Stages
0–5, both sample-set arms) was run end-to-end with **only** this environment on
`PATH` — V1's vvg-box stripped out entirely:

```bash
env -i HOME="$HOME" PATH=/usr/bin:/bin:/usr/sbin:/sbin \
  pixi run -- snakemake --configfile tests/tiny_cohort/config.yaml --cores 4
```

91 of 91 steps completed. Of the 118 text outputs comparable against the
previous V1-env run, **115 are byte-identical**. The three that are not:

- `outputs/moi/declonalization_comparison.tsv` — the stored baseline predates
  commit `e99fe67` ("Stage 2: composite country×geography Fws summary"), which
  changed the row labels from `country=CountryB` to `CountryB`. Regenerating
  this file under the V1 env produces exactly what this environment produces.
  Stale baseline, not an environment difference.
- `outputs/structure/{full,unique}/admixture/cv_error.tsv`, K=3 row only —
  ADMIXTURE is run with `-j2`, and its multi-threaded cross-validation is not
  reproducible run to run. Three consecutive K=3 runs of the *same* binary in
  the *same* environment gave 1.13692, 1.13653 and 1.13834. The K=2 row is
  stable across all runs. Thread nondeterminism, not an environment difference.
  (Worth knowing independently of the env: CV error is a model-selection
  diagnostic, so treat differences in its third decimal as noise. Run ADMIXTURE
  single-threaded if you need a byte-reproducible CV curve.)

The test suites also pass under this environment alone: 89 R unit tests and 44
pytest tests, the latter including the Stage 5 introgression positive control.

## Changing the environment

Add or bump a dependency in `../pixi.toml`, then:

```bash
pixi lock          # re-solve; commit the updated pixi.lock alongside pixi.toml
pixi install
pixi run check-env
```

Commit `pixi.toml` and `pixi.lock` together — a lock that does not match its
manifest is worse than no lock. If you add a package that `check_env.sh` should
watch, add it to the list in that script too, and record anything notable in
the version table above.

For the non-conda six, the pins live at the top of `postinstall.sh`.
