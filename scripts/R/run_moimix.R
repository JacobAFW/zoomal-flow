#!/usr/bin/env Rscript
# run_moimix.R  (agnostic port — generic, identical logic to V1)
# --------------------------------------------------------------------------
# Compute per-sample BAF matrix and Fws using moimix, starting from a GDS
# file. Writes the seam output fws_MOI.tsv (sample, Proportion) and a
# long-format BAF_dataframe.tsv for the Quarto narrative.
#
# Usage:
#   Rscript run_moimix.R <in.gds> <fws_out.tsv> <baf_out.tsv> <seed>
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(SeqArray)
  library(moimix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript run_moimix.R <in.gds> <fws_out.tsv> <baf_out.tsv> <seed>")
}
gds_in  <- args[1]
fws_out <- args[2]
baf_out <- args[3]
seed    <- as.integer(args[4])

stopifnot(file.exists(gds_in))
dir.create(dirname(fws_out), showWarnings = FALSE, recursive = TRUE)

message(sprintf("[run_moimix] opening %s", gds_in))
isolates <- seqOpen(gds_in)
on.exit(seqClose(isolates))

sample_ids <- seqGetData(isolates, "sample.id")
message(sprintf("[run_moimix] samples: %d", length(sample_ids)))

coords <- getCoordinates(isolates)
message(sprintf("[run_moimix] variants: %d", nrow(coords)))

# BAF matrix
message("[run_moimix] computing BAF matrix")
isolate_baf <- bafMatrix(isolates)

baf_df <- isolate_baf$coords %>%
  cbind(as.data.frame(isolate_baf$baf_site) %>% setNames("baf_site")) %>%
  left_join(
    isolate_baf$baf_matrix %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("variant.id") %>%
      mutate(variant.id = as.numeric(variant.id)),
    by = "variant.id"
  )

write_tsv(baf_df, baf_out, na = "NA")
message(sprintf("[run_moimix] wrote BAF to %s (%d rows)", baf_out, nrow(baf_df)))

# Fws
message(sprintf("[run_moimix] computing Fws with set.seed(%d)", seed))
set.seed(seed)
fws_all <- getFws(isolates) %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  rename("Proportion" = ".")

write_tsv(fws_all, fws_out, na = "NA")
message(sprintf("[run_moimix] wrote Fws to %s (%d rows)", fws_out, nrow(fws_all)))

cat("\n=== Fws summary ===\n")
cat(sprintf("Samples: %d\n", nrow(fws_all)))
cat(sprintf("Fws: median %.4f, mean %.4f, min %.4f\n",
            median(fws_all$Proportion, na.rm = TRUE),
            mean(fws_all$Proportion, na.rm = TRUE),
            min(fws_all$Proportion, na.rm = TRUE)))
cat(sprintf("Polyclonal (Fws < 0.95): %d (%.1f%%)\n",
            sum(fws_all$Proportion < 0.95, na.rm = TRUE),
            100 * mean(fws_all$Proportion < 0.95, na.rm = TRUE)))
