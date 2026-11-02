#!/usr/bin/env Rscript
# bamsort_process.R
# Standalone script: ensure per-sample raw_cellIDs CSVs exist (generating from
# 10X if missing), then run process_bamsort on Ca*Re alignment CSVs.
# All configuration is read from environment variables set by the calling
# shell script. Designed to be safe for parallel per-sample jobs: each sample
# writes only to its own ${sample}_raw_cellIDs.csv (no shared combined file).
#
# Required env vars:
#   FUNCTIONS_R              path to bamsort_functions.R (must define process_bamsort)
#   RAW_CELLIDS_OUTPUT_DIR   dir to write per-sample raw_cellIDs CSVs
#   RAW_CELLIDS_DATASETS     |-separated 10X dataset names (folder suffix)
#   RAW_CELLIDS_SAMPLES      |-separated short sample names (parallel)
#   RAW_CELLIDS_FOLDERS      |-separated parent dirs of the 10X counts (parallel)
#   PROCESS_SAMPLES          |-separated sample names to process
#   PROCESS_INPUT_DIRS       |-separated dirs containing *_bamsort_alignment_*.csv
#   PROCESS_OUTPUT_DIR       dir to write *_UMI_read_counts_*.csv
#   GENOME_PREFIX            e.g. "Mmul_10_mac239" (10X folder prefix)
#
# Optional env vars:
#   PROCESS_CORRECTION_DIR   splice correction CSVs (empty/unset = no correction)
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

# ========== READ CONFIG =======================================================

functions_R              <- get_env_required("FUNCTIONS_R")
raw_cellIDs_output_dir   <- get_env_required("RAW_CELLIDS_OUTPUT_DIR")
raw_cellIDs_datasets     <- get_env_list("RAW_CELLIDS_DATASETS")
raw_cellIDs_samples_vec  <- get_env_list("RAW_CELLIDS_SAMPLES")
raw_cellIDs_folders      <- get_env_list("RAW_CELLIDS_FOLDERS")
process_samples          <- get_env_list("PROCESS_SAMPLES")
process_input_dirs       <- get_env_list("PROCESS_INPUT_DIRS")
process_output_dir       <- get_env_required("PROCESS_OUTPUT_DIR")
genome_prefix            <- get_env_required("GENOME_PREFIX")
process_correction_dir   <- get_env_optional("PROCESS_CORRECTION_DIR", NULL)
min_umi_count            <- as.integer(get_env_optional("MIN_UMI_COUNT", "100"))

stopifnot(
  length(raw_cellIDs_datasets) == length(raw_cellIDs_samples_vec),
  length(raw_cellIDs_datasets) == length(raw_cellIDs_folders)
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
    cat("Loading existing raw_cellIDs for", sample_name, ":", sample_file, "\n")
    df <- read.csv(sample_file, stringsAsFactors = FALSE)
    raw_cellIDs[[sample_name]] <- df$cellID
    next
  }

  rawDataFolder <- paste0(raw_cellIDs_folders[i],
                          genome_prefix, "_", dataset_name,
                          "/outs/raw_feature_bc_matrix/")

  tryCatch({
    cat("Generating raw_cellIDs for", sample_name, "from:", rawDataFolder, "\n")
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
    cat("Loaded", length(cell_ids), "raw cellIDs for sample", sample_name, "\n")

    raw_cellIDs[[sample_name]] <- cell_ids

    write.csv(data.frame(cellID = cell_ids), sample_file, row.names = FALSE)
    cat("Saved", length(cell_ids), "cell IDs to", sample_file, "\n")

    rm(rawdata, umi_counts, cell_ids)
    gc(verbose = FALSE)
  }, error = function(e) {
    cat("Error processing", dataset_name, ":", conditionMessage(e), "\n")
  })
}

# ========== STEP 2: process_bamsort ========================================

cat("\nRunning process_bamsort on", length(process_samples), "sample(s):",
    paste(process_samples, collapse = ", "), "\n")
cat("Input dirs:\n  ", paste(process_input_dirs, collapse = "\n   "), "\n")
cat("Output dir:", process_output_dir, "\n")
cat("Correction dir:",
    ifelse(is.null(process_correction_dir), "(none)", process_correction_dir), "\n\n")

process_bamsort(
  samples_list   = process_samples,
  input.dir      = process_input_dirs,
  output.dir     = process_output_dir,
  raw_cellIDs    = raw_cellIDs,
  correction_dir = process_correction_dir
)

cat("\nbamsort_process.R complete.\n")
