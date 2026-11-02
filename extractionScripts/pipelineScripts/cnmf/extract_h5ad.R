#!/usr/bin/env Rscript
# extract_h5ad.R
# Convert a Seurat .rds or a cells x genes counts .tsv to AnnData .h5ad,
# suitable as direct cNMF input. Mirrors the writeH5AD pattern already used
# in the project's Seurat pipeline (see plotScripts/functions.R) so no Python-side
# conversion is needed.
#
# Usage:
#   Rscript extract_h5ad.R <input_rds_or_tsv> <output_h5ad>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript extract_h5ad.R <input_rds_or_tsv> <output_h5ad>")
}
input       <- args[1]
output_h5ad <- args[2]

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(zellkonverter)
  library(Matrix)
})

cat("Reading:", input, "\n")

if (grepl("\\.rds$", input, ignore.case = TRUE)) {
  obj <- readRDS(input)
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj, assay = "RNA")
  sce <- SingleCellExperiment(
    assays  = list(counts = GetAssayData(obj, assay = "RNA", layer = "counts")),
    colData = obj@meta.data
  )
} else if (grepl("\\.tsv$", input, ignore.case = TRUE)) {
  df <- read.table(input, sep = "\t", header = TRUE,
                   check.names = FALSE, row.names = "cellID")
  counts_mat <- as(t(as.matrix(df)), "sparseMatrix")
  sce <- SingleCellExperiment(assays = list(counts = counts_mat))
} else {
  stop("Input must be .rds or .tsv: ", input)
}

dir.create(dirname(output_h5ad), recursive = TRUE, showWarnings = FALSE)
suppressMessages(writeH5AD(sce, file = output_h5ad))
cat("Wrote:", output_h5ad, "(", ncol(sce), "cells x", nrow(sce), "genes)\n")
