#!/usr/bin/env Rscript
# bamsort_process_host.R
# Standalone script: ensure per-sample raw_cellIDs CSVs exist (generating from
# 10X if missing), then run process_bamsort_host on host-target alignment CSVs.
# Mirrors bamsort_process.R, the viral-target equivalent.
#
# Required env vars:
#   FUNCTIONS_R              path to bamsort_functions.R (must define process_bamsort_host)
#   RAW_CELLIDS_OUTPUT_DIR   dir to write per-sample raw_cellIDs CSVs
#   RAW_CELLIDS_DATASETS     |-separated 10X dataset names (folder suffix)
#   RAW_CELLIDS_SAMPLES      |-separated short sample names (parallel)
#   RAW_CELLIDS_FOLDERS      |-separated parent dirs of the 10X counts (parallel)
#   PROCESS_SAMPLES          |-separated sample names to process
#   PROCESS_INPUT_DIRS       |-separated dirs containing *_bamsort_alignment_*.csv
#   PROCESS_OUTPUT_DIR       dir to write *_UMI_read_counts_*.csv
#   GENOME_PREFIX            e.g. "GRCm39_EGL002" (10X folder prefix)
#
# Optional env vars:
#   PROCESS_NICKNAMES        |-separated, one per PROCESS_SAMPLES
#                            (defaults to PROCESS_SAMPLES; used by
#                            process_bamsort_host for alignment-file regex
#                            matching, e.g. nickname "KO" matches "*KO*_alignment_*.csv")
#   PROCESS_TARGETS          |-separated allowed targets, restricts the
#                            (target, alignment_part) categories accepted
#                            by extract_category_host. Unset = all targets.
#   MIN_UMI_COUNT            min UMI count to call a cell (default 100)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(Seurat)
})

# ========== ENV VAR HELPERS ===================================================

get_env_required <- function(name) {
  v <- Sys.getenv(name, unset = NA)
  if (is.na(v) || v == "") {
    stop("Required environment variable not set: ", name)
  }
  v
}

get_env_optional <- function(name, default = NULL) {
  v <- Sys.getenv(name, unset = NA)
  if (is.na(v) || v == "") return(default)
  v
}

get_env_list <- function(name, sep = "|") {
  strsplit(get_env_required(name), sep, fixed = TRUE)[[1]]
}

get_env_list_optional <- function(name, sep = "|", default = NULL) {
  v <- Sys.getenv(name, unset = NA)
  if (is.na(v) || v == "") return(default)
  strsplit(v, sep, fixed = TRUE)[[1]]
}

# ========== READ CONFIG =======================================================

functions_R              <- get_env_required("FUNCTIONS_R")
raw_cellIDs_output_dir   <- get_env_required("RAW_CELLIDS_OUTPUT_DIR")
raw_cellIDs_datasets     <- get_env_list("RAW_CELLIDS_DATASETS")
raw_cellIDs_samples_vec  <- get_env_list("RAW_CELLIDS_SAMPLES")
raw_cellIDs_folders      <- get_env_list("RAW_CELLIDS_FOLDERS")
process_samples          <- get_env_list("PROCESS_SAMPLES")
process_nicknames        <- get_env_list_optional("PROCESS_NICKNAMES",
                                                  default = process_samples)
process_input_dirs       <- get_env_list("PROCESS_INPUT_DIRS")
process_output_dir       <- get_env_required("PROCESS_OUTPUT_DIR")
genome_prefix            <- get_env_required("GENOME_PREFIX")
process_targets          <- get_env_list_optional("PROCESS_TARGETS",
                                                  default = NULL)
min_umi_count            <- as.integer(get_env_optional("MIN_UMI_COUNT", "100"))

stopifnot(
  length(raw_cellIDs_datasets) == length(raw_cellIDs_samples_vec),
  length(raw_cellIDs_datasets) == length(raw_cellIDs_folders),
  length(process_samples)      == length(process_nicknames),
  length(process_samples)      == length(raw_cellIDs_samples_vec)
)

source(functions_R)

dir.create(raw_cellIDs_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(process_output_dir,     recursive = TRUE, showWarnings = FALSE)

# ========== STEP 1: ensure per-sample raw_cellIDs CSVs exist =================

raw_cellIDs <- list()

for (i in seq_along(raw_cellIDs_samples_vec)) {
  dataset_name <- raw_cellIDs_datasets[i]
  sample_name  <- raw_cellIDs_samples_vec[i]
  sample_file  <- file.path(raw_cellIDs_output_dir,
                            paste0(sample_name, "_raw_cellIDs.csv"))

  if (file.exists(sample_file)) {
    df <- read.csv(sample_file, stringsAsFactors = FALSE)
    raw_cellIDs[[sample_name]] <- df$cellID
    cat("Loaded raw_cellIDs for", sample_name, ":", length(df$cellID), "cells\n")
    next
  }

  rawDataFolder <- paste0(raw_cellIDs_folders[i],
                          genome_prefix, "_", dataset_name,
                          "/outs/raw_feature_bc_matrix/")

  tryCatch({
    rawdata <- Read10X(rawDataFolder)
    # Feature-barcode runs (Gene Expression + Antibody Capture/hashtags) come
    # back as a list of matrices, one per library type. Cells are called on the
    # transcriptome, so take that matrix; plain TenX runs return it directly.
    if (is.list(rawdata) && !is.data.frame(rawdata)) {
      if (!"Gene Expression" %in% names(rawdata)) {
        stop("Read10X returned library types with no 'Gene Expression': ",
             paste(names(rawdata), collapse = ", "))
      }
      rawdata <- rawdata[["Gene Expression"]]
    }
    umi_counts <- colSums(rawdata)
    cell_ids <- names(umi_counts[umi_counts > min_umi_count])
    raw_cellIDs[[sample_name]] <- cell_ids
    write.csv(data.frame(cellID = cell_ids), sample_file, row.names = FALSE)
    cat("Generated raw_cellIDs for", sample_name, ":", length(cell_ids), "cells\n")

    rm(rawdata, umi_counts, cell_ids)
    gc(verbose = FALSE)
  }, error = function(e) {
    cat("Error processing", dataset_name, ":", conditionMessage(e), "\n")
  })
}

# ========== STEP 2: process_bamsort_host ====================================

# Alias raw_cellIDs entries under the process-sample key so process_bamsort_host
# can look up raw_cellIDs[[process_sample]] even when the cellID list was loaded
# under a shared biological-sample key (e.g. RAW_CELLIDS_SAMPLES="GFP_invitro" but
# PROCESS_SAMPLES="PTPRC"). When they're already equal this is a no-op.
for (i in seq_along(process_samples)) {
  raw_key  <- raw_cellIDs_samples_vec[i]
  proc_key <- process_samples[i]
  if (proc_key != raw_key && !is.null(raw_cellIDs[[raw_key]])) {
    raw_cellIDs[[proc_key]] <- raw_cellIDs[[raw_key]]
    cat("Aliased raw_cellIDs[[", proc_key, "]] <- [[", raw_key, "]]\n", sep = "")
  }
}

samples_df <- data.frame(
  samples   = process_samples,
  nicknames = process_nicknames,
  stringsAsFactors = FALSE
)

process_bamsort_host(
  samples_df  = samples_df,
  input.dir   = process_input_dirs,
  output.dir  = process_output_dir,
  raw_cellIDs = raw_cellIDs,
  targets     = process_targets
)
