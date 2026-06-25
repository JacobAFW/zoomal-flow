#!/usr/bin/env Rscript
# validate_metadata.R
# --------------------------------------------------------------------------
# Load the user's samples.tsv, check the role mapping against actual columns,
# and report VCF<->metadata concordance. Writes a canonical, role-renamed copy
# to outputs/metadata/samples.tsv (canonical column names = role names) so
# downstream rules can read by role without consulting the config again.
#
# Hard errors on:
#   - sample_id role missing or its column absent from the table.
# Warnings (logged, never crash) on:
#   - any optional role whose configured column is absent from the table;
#   - VCF samples not in the table, or table samples not in the VCF.
#
# Usage:
#   Rscript validate_metadata.R <samples.tsv> <vcf_samples.txt> <roles.json> \
#       <out_tsv> <out_log>
# Where <roles.json> is a JSON dump of config.metadata.roles (written by the
# Snakefile to avoid teaching R to parse the full YAML config).
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(readr)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(paste(
    "Usage: Rscript validate_metadata.R",
    "<samples.tsv> <vcf_samples.txt> <roles.json> <out_tsv> <out_log>"
  ))
}
samples_in    <- args[1]
vcf_samples_in<- args[2]
roles_json    <- args[3]
out_tsv       <- args[4]
out_log       <- args[5]

dir.create(dirname(out_tsv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_log), showWarnings = FALSE, recursive = TRUE)

log_lines <- character()
log_emit <- function(...) {
  msg <- paste0(...)
  log_lines <<- c(log_lines, msg)
  message(msg)
}

roles <- fromJSON(roles_json, simplifyVector = TRUE)
# JSON nulls come back as missing list entries; normalise to NA character.
expected_roles <- c("sample_id", "group", "geography", "country",
                    "host", "date", "case_control")
for (r in expected_roles) {
  if (is.null(roles[[r]])) roles[[r]] <- NA_character_
}

if (is.na(roles$sample_id) || !nzchar(roles$sample_id)) {
  stop("config.metadata.roles.sample_id is required and must be a column name")
}

log_emit("[validate_metadata] reading ", samples_in)
meta <- read_tsv(samples_in, show_col_types = FALSE)
log_emit(sprintf("[validate_metadata] %d rows x %d cols", nrow(meta), ncol(meta)))

if (!(roles$sample_id %in% names(meta))) {
  stop(sprintf(
    "sample_id role maps to column '%s' but that column is not in %s. Available: %s",
    roles$sample_id, samples_in, paste(names(meta), collapse = ", ")
  ))
}

# Optional roles: log absences, don't crash. Build a column->role rename map
# of only the roles that resolve to a real column.
rename_map <- list()
rename_map[[roles$sample_id]] <- "sample_id"

optional <- setdiff(expected_roles, "sample_id")
for (r in optional) {
  col <- roles[[r]]
  if (is.na(col) || !nzchar(col)) {
    log_emit(sprintf("[role] %-12s : null  → dependent analyses will be skipped", r))
    next
  }
  if (!(col %in% names(meta))) {
    log_emit(sprintf(
      "[role] %-12s : column '%s' absent → dependent analyses will be skipped",
      r, col
    ))
    next
  }
  rename_map[[col]] <- r
}

# Apply rename: column 'roles$group' becomes 'group', etc. Other columns pass
# through untouched so cohort-specific extras survive into outputs/.
canonical <- meta
for (old_name in names(rename_map)) {
  new_name <- rename_map[[old_name]]
  if (old_name == new_name) next
  if (new_name %in% names(canonical) && new_name != old_name) {
    log_emit(sprintf(
      "[role] WARN: canonical name '%s' already present in table; appending '_orig' to the existing column",
      new_name
    ))
    names(canonical)[names(canonical) == new_name] <- paste0(new_name, "_orig")
  }
  names(canonical)[names(canonical) == old_name] <- new_name
}

# VCF concordance
vcf_samples <- readLines(vcf_samples_in)
meta_samples <- as.character(canonical$sample_id)
in_vcf_not_meta <- setdiff(vcf_samples, meta_samples)
in_meta_not_vcf <- setdiff(meta_samples, vcf_samples)

log_emit(sprintf("[concordance] VCF samples: %d", length(vcf_samples)))
log_emit(sprintf("[concordance] metadata samples: %d", length(meta_samples)))
log_emit(sprintf("[concordance] in VCF, not in metadata: %d", length(in_vcf_not_meta)))
if (length(in_vcf_not_meta) > 0) {
  log_emit("  - ", paste(head(in_vcf_not_meta, 20), collapse = ", "),
           if (length(in_vcf_not_meta) > 20) " ..." else "")
}
log_emit(sprintf("[concordance] in metadata, not in VCF: %d", length(in_meta_not_vcf)))
if (length(in_meta_not_vcf) > 0) {
  log_emit("  - ", paste(head(in_meta_not_vcf, 20), collapse = ", "),
           if (length(in_meta_not_vcf) > 20) " ..." else "")
}

write_tsv(canonical, out_tsv)
writeLines(log_lines, out_log)
log_emit(sprintf("[validate_metadata] wrote %s + %s", out_tsv, out_log))
