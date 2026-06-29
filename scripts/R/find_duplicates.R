#!/usr/bin/env Rscript
# find_duplicates.R  (agnostic port — config-driven dedup pattern)
# --------------------------------------------------------------------------
# Identify replicate samples and write the higher-missingness replicates
# to a PLINK --remove file (FID IID per line).
#
# V1 hardcoded `str_remove(V1, "_DK.*")` for the Indo cohort's Illumina lane
# suffix. The agnostic version takes the regex as the 4th argument:
#   - "" (empty) → no strip; still normalises underscores/hyphens.
#   - non-empty  → applied via str_remove() before normalisation.
#
# Direction (kept from V1): we `slice_max(F_MISS)` — the file lists the
# WORSE (higher-missing) replicate per dup group, which is what PLINK's
# --remove drops. The "lowest-missing replicate" wording in legacy
# comments refers to the one that's KEPT.
#
# Usage:
#   Rscript find_duplicates.R <Pk.fam> <Pk.smiss> <Pk.dups> <id_pattern>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript find_duplicates.R <Pk.fam> <Pk.smiss> <Pk.dups> <id_pattern>")
}
fam_in   <- args[1]
smiss_in <- args[2]
dups_out <- args[3]
pattern  <- args[4]      # "" → don't strip

fam <- read.table(fam_in, header = FALSE, stringsAsFactors = FALSE) %>%
  as_tibble()

dups <- fam %>%
  mutate(
    Sample = if (nzchar(pattern)) str_remove(V1, pattern) else V1,
    Sample = str_replace(Sample, "_", ""),
    Sample = str_replace(Sample, "-", "")
  ) %>%
  group_by(Sample) %>%
  filter(n() > 1) %>%
  ungroup()

smiss <- read_table(smiss_in, show_col_types = FALSE) %>%
  rename(FID = `#FID`) %>%
  select(FID, IID, F_MISS)

dups <- dups %>%
  left_join(smiss, by = c("V1" = "IID")) %>%
  group_by(Sample) %>%
  slice_max(F_MISS, with_ties = FALSE) %>%
  ungroup()

out <- dups %>% select(V1, V2)

write.table(out, dups_out, quote = FALSE, sep = " ",
            col.names = FALSE, row.names = FALSE)

message(sprintf("[find_duplicates] %d replicate(s) flagged for removal (pattern=%s)",
                nrow(out), if (nzchar(pattern)) pattern else "(none)"))
cat("\n=== Duplicates flagged (worst per dup-group) ===\n")
print(dups %>% select(V1, Sample, F_MISS), n = Inf)
