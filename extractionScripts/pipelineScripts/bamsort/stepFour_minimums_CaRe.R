#!/usr/bin/env Rscript
# stepFour_minimums.R  (Figure1 Ca*Re amplicon libraries)
# Applies the Ca*Re detection thresholds to the stepThree per-UMI tables and
# writes <group>/bamsort/reads/minimum/<sample>_CaRe_filtered.csv, the
# per-cell category table that Figure1.Rmd and SuppFigs.Rmd read.
#
# stepThree deliberately stops before thresholding so thresholds can be varied
# without re-running it; this is that thresholding step.
#
# min_reads and min_umi are set per sample, not per infection level. Values come
# from the CD3+CD4+-gated threshold grid in SuppFigs.Rmd, read out as the point
# where the SIV+ call count stops rising under read subsampling (Panel G). W2 is
# held to min_reads = 9 rather than 7: at 7 its Jaccard agreement with the TenX
# call on CD4 T cells (0.200) falls below Invitro's at 7/3 (0.262), and 9 is the
# lowest min_reads at min_umi = 3 that reaches it (0.292). The low-titer samples
# sit at the permissive min_reads = 2 because a higher cutoff leaves them with no
# cells at all.
#
# Pacute and W2 are the two viremia libraries and are held at the same
# min_reads deliberately: they are compared against each other throughout
# Figure 3 and SuppFig 3, so a threshold that differed between them would make
# every viremia-to-viremia contrast a comparison of detection stringency as
# well as of biology. Pacute is far from insensitive to the choice -- at 7 it
# keeps 812 cells against 321 at 9 -- so moving one of the two means moving
# both and regenerating everything downstream (see the header of
# .runs/invivo_t9_production.R in plotScripts for what that is).
#
# This file and SuppFigs.Rmd's saturation `sample_config` both declare these
# thresholds and have to agree; the saturation panel reads the grid that these
# values are chosen from.
#
# Run AFTER stepThree_process.sh has finished for every group.

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

care.root <- "/projects/b1042/GoyalLab/egrody/extractedData/Figure1/CaRe"

care.group <- c(D0 = "uninfected", W0 = "uninfected", D195 = "ART",
                Invitro = "invitro", Pacute = "viremia", W2 = "viremia")

# Thresholds. min_reads and min_umi are per sample; the rest are shared.
thresholds <- list(
  D0      = c(min_reads = 2, min_umi = 3),
  W0      = c(min_reads = 2, min_umi = 3),
  D195    = c(min_reads = 2, min_umi = 3),
  Invitro = c(min_reads = 7, min_umi = 3),
  Pacute  = c(min_reads = 9, min_umi = 3),
  W2      = c(min_reads = 9, min_umi = 3)
)
min_region_count <- 1
min_cells        <- 0
min_reads_cell   <- 1

missing <- setdiff(names(care.group), names(thresholds))
if (length(missing) > 0)
  stop("No threshold entry for: ", paste(missing, collapse = ", "))

for (sample_name in names(care.group)) {
  group      <- care.group[[sample_name]]
  min_reads  <- thresholds[[sample_name]][["min_reads"]]
  min_umi    <- thresholds[[sample_name]][["min_umi"]]
  reads.dir  <- file.path(care.root, group, "bamsort", "reads")
  output.dir <- file.path(reads.dir, "minimum")
  valid.dir  <- file.path(output.dir, "valid_UMI")
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(valid.dir,  recursive = TRUE, showWarnings = FALSE)

  umi_file <- file.path(reads.dir,
                        paste0(sample_name, "_UMI_read_counts_full.csv"))
  if (!file.exists(umi_file)) {
    cat("Skipping", sample_name, "- no stepThree table at", umi_file, "\n")
    next
  }

  cat(sprintf("\n=== %s (%s): min_reads = %d, min_umi = %d ===\n",
              sample_name, group, min_reads, min_umi))
  filtered <- set_minimums(
    "CaRe", umi_file,
    min_reads             = min_reads,
    min_region_count      = min_region_count,
    min_umi               = min_umi,
    min_cells             = min_cells,
    min_reads_cell        = min_reads_cell,
    save_valid_umi        = TRUE,
    valid_umi_output_file = file.path(valid.dir,
                              paste0(sample_name, "_valid_cellID_UMI.csv"))
  )

  out_file <- file.path(output.dir,
                        paste0(sample_name, "_CaRe_filtered.csv"))
  write.csv(filtered, out_file, row.names = FALSE)
  cat("Wrote", nrow(filtered), "cells to", out_file, "\n")

  writeLines(c(
    paste("min_reads =", min_reads),
    paste("min_region_count =", min_region_count),
    paste("min_umi =", min_umi),
    paste("min_cells =", min_cells),
    paste("min_reads_cell =", min_reads_cell),
    "mode = CaRe",
    paste("processed_date =", Sys.Date())
  ), file.path(output.dir, paste0(sample_name, "_filters.txt")))
}

cat("\nstepFour_minimums.R complete.\n")
