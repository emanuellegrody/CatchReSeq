#!/usr/bin/env Rscript
# stepFour_minimums_TenX.R  (Figure1 single-cell TenX libraries)
# Applies the TenX detection thresholds to the stepThree_process_TenX.sh per-UMI
# tables and writes, per sample:
#   <reads>/minimum/<sample>_GEX_filtered.csv              per-cell category table
#   <reads>/minimum/valid_UMI/<sample>_GEX_valid_cellID_UMI.csv
# The valid_UMI tables are what the analysis scripts semi_join against to filter
# the per-UMI tables down to detected molecules.
#
# The TenX arm keeps min_reads = 1: unlike the Ca*Re amplicon libraries, a TenX
# library has no targeted enrichment, so viral UMIs arrive at one or two reads
# and a read threshold removes the signal rather than the noise. min_reads_cell
# = 2 is the per-barcode QC the TenX arm has always used.
#
# Categories come from mode = "GEX", which carries the same six as CaRe so the
# two arms are directly comparable. set_minimums keeps only the listed
# categories, so any category left out of that list is dropped silently — the
# MS calls the stepThree correction makes have to be named there to survive.
#
# Run AFTER stepThree_process_TenX.sh has finished.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# bamsort_functions.R sits next to this script. Resolving it from the script's own
# location keeps the path working wherever the repo is cloned; Rscript exposes
# that location only through the command line it was invoked with.
script_path <- local({
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) sub("^--file=", "", hit[1]) else NA_character_
})
script_dir <- if (is.na(script_path)) getwd() else dirname(normalizePath(script_path))

source(file.path(script_dir, "bamsort_functions.R"))

tenx.root <- "/projects/b1042/GoyalLab/egrody/extractedData/Figure1/singleCell"

samples <- c("D0", "D195", "Invitro", "Pacute", "W0", "W2")

min_reads        <- 1
min_region_count <- 1
min_umi          <- 1
min_cells        <- 0
min_reads_cell   <- 2

reads.dir  <- file.path(tenx.root, "bamsort", "reads")
output.dir <- file.path(reads.dir, "minimum")
valid.dir  <- file.path(output.dir, "valid_UMI")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(valid.dir,  recursive = TRUE, showWarnings = FALSE)

for (sample_name in samples) {
  umi_file <- file.path(reads.dir,
                        paste0(sample_name, "_UMI_read_counts_full.csv"))
  if (!file.exists(umi_file)) {
    cat("Skipping", sample_name, "- no stepThree table at", umi_file, "\n")
    next
  }

  cat("\n===", sample_name, "===\n")
  filtered <- set_minimums(
    "GEX", umi_file,
    min_reads             = min_reads,
    min_region_count      = min_region_count,
    min_umi               = min_umi,
    min_cells             = min_cells,
    min_reads_cell        = min_reads_cell,
    save_valid_umi        = TRUE,
    valid_umi_output_file = file.path(valid.dir,
                              paste0(sample_name, "_GEX_valid_cellID_UMI.csv"))
  )

  out_file <- file.path(output.dir,
                        paste0(sample_name, "_GEX_filtered.csv"))
  write.csv(filtered, out_file, row.names = FALSE)
  cat("Wrote", nrow(filtered), "cells to", out_file, "\n")
  if ("MS" %in% colnames(filtered)) {
    cat("  MS+ cells:", sum(filtered$MS > 0), "\n")
  }

  writeLines(c(
    paste("min_reads =", min_reads),
    paste("min_region_count =", min_region_count),
    paste("min_umi =", min_umi),
    paste("min_cells =", min_cells),
    paste("min_reads_cell =", min_reads_cell),
    "mode = TenX",
    paste("processed_date =", Sys.Date())
  ), file.path(output.dir, "filters.txt"))
}

cat("\nstepFour_minimums_TenX.R complete.\n")
