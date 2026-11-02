#!/usr/bin/env Rscript
# stepFive_markers.R
# Writes <sample>_CD4T_counts.csv, the per-cell marker table the CD4 T gate
# reads, to the expressionDF tree the .Rmd figures load from.
# All configuration is read from environment variables set by the calling
# shell script, so one script serves every dataset.
#
# Required env vars:
#   FUNCTIONS_R       path to plotScripts/functions.R (must define cd4t_gate_genes,
#                     marker_counts and cd4t_call)
#   RDS_DIR           dir holding <sample>.rds Seurat objects
#   OUTPUT_DIR        dir to write <sample>_CD4T_counts.csv
#   MARKER_SAMPLES    |-separated sample names
#
# Two things differ from the older <sample>_CD4T.csv tables, which stay in place
# for the legacy EGS016/EGS024 notebooks:
#
#   raw counts, not normalized.  `cd4t_call` excludes on a count threshold
#   (cd4t_exclusion_min), and a log-normalized value of 2 corresponds to a
#   different number of molecules in every library. The old gate only asked
#   `> 0`, which is the same either way, so this did not matter until now.
#
#   the full `cd4t_gate_genes` set, not just CD3/CD4/CD8. The gate is negative
#   selection, so it needs every lineage's markers. A marker missing from the
#   table cannot exclude, and `cd4t_call` warns when most of them are.
#
# Run on a compute node -- the objects are large.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# ========== ENV VAR HELPERS ===================================================

get_env_required <- function(name) {
  v <- Sys.getenv(name, unset = NA)
  if (is.na(v) || v == "") {
    stop("Required environment variable not set: ", name)
  }
  v
}

get_env_list <- function(name, sep = "|") {
  strsplit(get_env_required(name), sep, fixed = TRUE)[[1]]
}

# ========== READ CONFIG =======================================================

functions_R <- get_env_required("FUNCTIONS_R")
rds.dir     <- get_env_required("RDS_DIR")
output.dir  <- get_env_required("OUTPUT_DIR")
samples     <- get_env_list("MARKER_SAMPLES")

# functions.R builds its default arguments from data.root, which config.R
# defines, so config.R has to be sourced first or sourcing functions.R fails
# with "object 'data.root' not found". config.R sits beside functions.R, so
# resolving it from there keeps this working wherever the repo is cloned.
source(file.path(dirname(functions_R), "config.R"))
source(functions_R)

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

for (sample_name in samples) {
  input_file <- file.path(rds.dir, paste0(sample_name, ".rds"))
  if (!file.exists(input_file)) {
    cat("Skipping", sample_name, "- no object at", input_file, "\n"); next
  }
  seurat_obj <- readRDS(input_file)

  # marker_counts() reads the counts layer across split layers, so this is the
  # same accessor the Seurat-object arm of the gate uses.
  counts <- marker_counts(seurat_obj, cd4t_gate_genes)
  sample_df <- data.frame(cellID = colnames(seurat_obj), t(counts),
                          sample = sample_name, check.names = FALSE)

  absent <- cd4t_gate_genes[rowSums(counts) == 0]
  if (length(absent) > 0) {
    cat(" ", sample_name, "- markers not detected:", paste(absent, collapse = ", "), "\n")
  }

  out_file <- file.path(output.dir, paste0(sample_name, "_CD4T_counts.csv"))
  write.csv(sample_df, out_file)
  cat("Wrote", nrow(sample_df), "cells to", out_file, "\n")

  # What the gate makes of it, so the log carries the numbers without a second
  # pass over the object.
  call <- cd4t_call(sample_df)
  cat(sprintf("  CD3+: %d | excluded as another lineage: %d | gate passes: %d of %d\n",
              sum(call$cd3), sum(call$excluded), sum(call$pass), length(call$pass)))

  rm(seurat_obj, counts); gc(verbose = FALSE)
}

cat("\nstepFive_markers.R complete.\n")
