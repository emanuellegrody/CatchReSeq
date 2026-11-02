#!/usr/bin/env Rscript
# memoryseq_variability.R
# Calls over-dispersed ("potentially heritable") genes from the LD002 clone bulk
# RNA-seq TPM matrix, the upstream half of Figure4.Rmd Panel J.
#
# All configuration is read from environment variables set by the calling shell
# script. The defaults reproduce the list the published panel was built from, and
# three of them encode choices that should be reconsidered before the panel is
# quoted. Each is flagged at its own definition below.
#
# Required env vars:
#   TPM_MATRIX        star_salmon/salmon.merged.gene_tpm.tsv: a gene_id column,
#                     a gene_name column, one column per library
#   OUTPUT_DIR        dir to write the fit and hit tables to
#
# Optional env vars:
#   COUNTS_MATRIX     star_salmon/salmon.merged.gene_counts.tsv. When set, the
#                     expression floor is applied to it; unset, to TPM.
#   SAMPLE_PATTERN    regex selecting the library columns (default "^Clone_")
#   CONTROL_SAMPLES   |-separated libraries to fit as a separate "Controls" arm
#   HITS_FROM         "both" or "ExpSample": which arms reach the hit list
#                     (default "both")
#   MEAN_FLOOR        expression floor on the mean across libraries (default 1)
#   MEAN_TPM_CUTOFF   mean-TPM floor for a hit (default 1.5)
#   PERCENTILE_CUTOFF residual percentile defining a hit (default 0.98)
#   EXCLUDE_PATTERN   regex on GeneSymbol dropped from the filtered list
#                     (default "^(IG|TR)"); set empty to keep every hit
#
# ---- Method ------------------------------------------------------------------
#
# MemorySeq (Shaffer et al. 2020, Cell 182:947) asks which genes vary more
# across single-cell-derived clones than the mean-variance relationship of the
# transcriptome predicts. A clone is many generations of expansion from one
# cell, so expression that is stably inherited shows up as between-clone
# variance while transient fluctuation averages out within each clone.
#
#   1. Drop genes under the expression floor, which otherwise dominate the
#      high-CV tail on sampling noise alone.
#   2. Per gene and per arm, Mean_tpm and CV_tpm = sd(tpm) / mean(tpm).
#   3. glm(CV_tpm ~ log2(Mean_tpm), family = poisson(link = "log")) as the
#      expected mean-CV trend, then take the genes whose deviance residual sits
#      above PERCENTILE_CUTOFF and whose mean clears MEAN_TPM_CUTOFF. The
#      Poisson family is applied to a continuous response here, which is the
#      paper's own choice: it is a log-link mean-variance trend line, not a
#      likelihood for count data.
#
# ---- The receptor-locus exclusion -------------------------------------------
#
# Panel J reads potential_heritable_genes_filtered.csv, which is this script's
# hit list with EXCLUDE_PATTERN applied. No separate script for that step was
# ever found; the rule is recoverable from the file as dropping the hits whose
# GeneSymbol matches `^(IG|TR)`.
#
# A panel of single-cell-derived T-cell clones has a different TCR and immunoglobulin
# rearrangement in every clone, so V, D, J and constant-region expression is
# perfectly clone-specific by genomic rearrangement rather than by gene
# regulation. Those loci top the residual ranking in every configuration, and
# keeping them would mean calling a clone's own receptor sequence a heritable
# expression state.
#
# `^(IG|TR)` is broader than it needs to be: it also removes TROAP and TRIP13,
# which are not receptor genes. EXCLUDE_PATTERN is settable for that reason; the
# default is the pattern that reproduces the published panel.

suppressPackageStartupMessages({
  library(dplyr)
})

# ========== ENV VAR HELPERS ===================================================

get_env_required <- function(name) {
  value <- Sys.getenv(name, unset = NA)
  if (is.na(value) || value == "") {
    stop("Required environment variable not set: ", name)
  }
  value
}

get_env_numeric <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (value == "") return(default)
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed)) stop(name, " must be numeric, got: ", value)
  parsed
}

get_env_list <- function(name, sep = "|") {
  value <- Sys.getenv(name, unset = "")
  if (value == "") character(0) else strsplit(value, sep, fixed = TRUE)[[1]]
}

# ========== READ CONFIG =======================================================

tpm_matrix_path    <- get_env_required("TPM_MATRIX")
output.dir         <- get_env_required("OUTPUT_DIR")
counts_matrix_path <- Sys.getenv("COUNTS_MATRIX", unset = "")
sample_pattern     <- Sys.getenv("SAMPLE_PATTERN", unset = "^Clone_")
control_samples    <- get_env_list("CONTROL_SAMPLES")
hits_from          <- Sys.getenv("HITS_FROM", unset = "both")

mean_floor        <- get_env_numeric("MEAN_FLOOR", 1)
mean_tpm_cutoff   <- get_env_numeric("MEAN_TPM_CUTOFF", 1.5)
percentile_cutoff <- get_env_numeric("PERCENTILE_CUTOFF", 0.98)
exclude_pattern   <- Sys.getenv("EXCLUDE_PATTERN", unset = "^(IG|TR)")

if (!hits_from %in% c("both", "ExpSample")) {
  stop("HITS_FROM must be \"both\" or \"ExpSample\", got: ", hits_from)
}

if (!grepl("/$", output.dir)) output.dir <- paste0(output.dir, "/")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(tpm_matrix_path)) stop("TPM matrix not found: ", tpm_matrix_path)

cat("TPM matrix:   ", tpm_matrix_path, "\n")
cat("Output dir:   ", output.dir, "\n")

# ========== LOAD ==============================================================

tpm <- read.delim(tpm_matrix_path, check.names = FALSE, stringsAsFactors = FALSE)
for (required_column in c("gene_id", "gene_name")) {
  if (!required_column %in% colnames(tpm)) {
    stop("TPM matrix has no '", required_column, "' column. Columns found: ",
         paste(colnames(tpm), collapse = ", "))
  }
}

library_columns <- grep(sample_pattern,
                        setdiff(colnames(tpm), c("gene_id", "gene_name")),
                        value = TRUE)
if (length(library_columns) == 0) {
  stop("SAMPLE_PATTERN '", sample_pattern, "' matched no columns in ",
       tpm_matrix_path)
}

gene_symbol_of <- setNames(tpm$gene_name, tpm$gene_id)

# One arm per fit. CONTROL_SAMPLES splits the panel in two, which is how the
# published list was made.
#
# DEFECT 1, carried by the default in memoryseq_variability.sh. The four
# libraries named there are clones, not controls: no aliquot of a pooled
# population was ever sequenced for LD002. MemorySeq's controls exist to carry
# the technical and sampling variance the clones also carry, so that the hit
# list can be read against them. Splitting the clone panel instead fits one
# trend to 11 clones and another to 4, and a CV estimated from 4 libraries is
# both noisier and systematically different from one estimated from 11. Neither
# arm is a control for the other. Set CONTROL_SAMPLES empty to fit the panel
# whole, which is the only defensible reading when there are no controls.
unknown_controls <- setdiff(control_samples, library_columns)
if (length(unknown_controls) > 0) {
  stop("CONTROL_SAMPLES names libraries not in the matrix: ",
       paste(unknown_controls, collapse = ", "))
}
experimental_samples <- setdiff(library_columns, control_samples)
if (length(experimental_samples) == 0) {
  stop("CONTROL_SAMPLES claims every library; no ExpSample arm is left.")
}

arms <- list(ExpSample = experimental_samples)
if (length(control_samples) > 0) arms$Controls <- control_samples

cat("Libraries (", length(library_columns), "): ",
    paste(library_columns, collapse = ", "), "\n", sep = "")
for (arm_label in names(arms)) {
  cat("  ", arm_label, " (", length(arms[[arm_label]]), "): ",
      paste(arms[[arm_label]], collapse = ", "), "\n", sep = "")
}
if (length(control_samples) > 0) {
  cat("  NOTE: the Controls arm above is a subset of the clone panel, not a",
      "control population.\n")
}

# ========== EXPRESSION FLOOR ==================================================

# DEFECT 2. The published run applied the floor to salmon.merged.gene_counts.tsv
# -- raw counts, which the original script read into a column named "rpm" -- so
# the filter was "mean raw count > 1", not the 1 RPM the method asks for. Raw
# counts scale with library size and gene length, so the floor sat at a
# different expression level for every gene. Leave COUNTS_MATRIX unset to apply
# the floor to TPM instead, which is what MEAN_FLOOR is in units of.
if (counts_matrix_path != "") {
  if (!file.exists(counts_matrix_path)) {
    stop("COUNTS_MATRIX not found: ", counts_matrix_path)
  }
  counts <- read.delim(counts_matrix_path, check.names = FALSE,
                       stringsAsFactors = FALSE)
  if (!identical(counts$gene_id, tpm$gene_id)) {
    stop("COUNTS_MATRIX and TPM_MATRIX do not carry the same gene_id rows ",
         "in the same order.")
  }
  floor_values <- rowMeans(as.matrix(counts[, library_columns, drop = FALSE]))
  cat("Expression floor: mean raw count >", mean_floor,
      "(COUNTS_MATRIX set)\n")
} else {
  floor_values <- rowMeans(as.matrix(tpm[, library_columns, drop = FALSE]))
  cat("Expression floor: mean TPM >", mean_floor, "\n")
}

tpm_above_floor <- tpm[floor_values > mean_floor, ]
cat("Genes in matrix:", nrow(tpm), "| above floor:", nrow(tpm_above_floor), "\n")

# ========== VARIABILITY METRICS AND FIT =======================================

calc_variability_metrics <- function(columns, arm_label) {
  expression_matrix <- as.matrix(tpm_above_floor[, columns, drop = FALSE])
  mean_tpm <- rowMeans(expression_matrix)
  data.frame(gene_id    = tpm_above_floor$gene_id,
             GeneSymbol = gene_symbol_of[tpm_above_floor$gene_id],
             controls   = arm_label,
             Mean_tpm   = mean_tpm,
             CV_tpm     = apply(expression_matrix, 1, sd) / mean_tpm,
             stringsAsFactors = FALSE)
}

fit_cv_mean_trend <- function(metrics, arm_label) {
  fittable <- metrics[complete.cases(metrics), ]
  dropped <- nrow(metrics) - nrow(fittable)
  if (dropped > 0) {
    cat("[", arm_label, "] dropped ", dropped,
        " gene(s) with a non-finite mean or CV\n", sep = "")
  }

  # A Poisson family on a continuous response warns "non-integer x" once per
  # observation -- tens of thousands of lines that would bury any warning worth
  # reading. Only that message is muffled, and it is counted so the log still
  # records that the response was continuous; everything else glm says gets
  # through, including a non-convergence warning.
  non_integer_warnings <- 0L
  model <- withCallingHandlers(
    glm(CV_tpm ~ log2(Mean_tpm), family = poisson(link = "log"),
        data = fittable),
    warning = function(w) {
      if (grepl("non-integer", conditionMessage(w))) {
        non_integer_warnings <<- non_integer_warnings + 1L
        invokeRestart("muffleWarning")
      }
    })
  cat("[", arm_label, "] fit ", nrow(fittable), " genes; glm reported ",
      non_integer_warnings, " non-integer response warning(s)\n", sep = "")

  fittable$fitted_y <- fitted(model)
  fittable$resid    <- resid(model)
  fittable
}

fit_data <- do.call(rbind, lapply(names(arms), function(arm_label) {
  fit_cv_mean_trend(calc_variability_metrics(arms[[arm_label]], arm_label),
                    arm_label)
}))

# ========== HIT CALL ==========================================================

# The cutoff is set on the ExpSample arm and applied to every arm, so that a
# second arm is scored against the line the first one defined.
residual_cutoff <- quantile(
  fit_data$resid[fit_data$controls == "ExpSample" &
                 fit_data$Mean_tpm > mean_tpm_cutoff],
  percentile_cutoff
)
fit_data$cutoff <- (fit_data$Mean_tpm > mean_tpm_cutoff) &
                   (fit_data$resid > residual_cutoff)

cat("\nResidual cutoff at the", percentile_cutoff, "percentile:",
    round(residual_cutoff, 4), "\n")
cat("Genes tested (mean TPM >", mean_tpm_cutoff, "):",
    sum(fit_data$controls == "ExpSample" &
        fit_data$Mean_tpm > mean_tpm_cutoff), "\n")
for (arm_label in names(arms)) {
  cat("Hits in", arm_label, ":",
      sum(fit_data$cutoff & fit_data$controls == arm_label), "\n")
}

# DEFECT 3. The published run kept hits from every arm. The paper's own script
# keeps only the ExpSample arm (`filter(fit_data, cutoff==TRUE &
# controls=='ExpSample')`), because a hit in the control arm is by construction
# a gene the controls failed to hold steady -- evidence against heritability,
# not for it. With the LD002 split that reasoning does not even apply, since
# both arms are clones; what it means in practice is that most of the published
# list was called from the 4-library fit. Set HITS_FROM=ExpSample to keep one
# arm.
heritable_genes <- fit_data[fit_data$cutoff, ]
if (hits_from == "ExpSample") {
  heritable_genes <- heritable_genes[heritable_genes$controls == "ExpSample", ]
}
heritable_genes <- heritable_genes[order(-heritable_genes$resid), ]

cat("Hit list (HITS_FROM =", hits_from, "):", nrow(heritable_genes), "rows,",
    length(unique(heritable_genes$GeneSymbol)), "distinct symbols\n")
cat("Symbols that are bare Ensembl ids (no symbol in the annotation):",
    sum(grepl("^ENSG", unique(heritable_genes$GeneSymbol))), "\n")

# ========== RECEPTOR-LOCUS EXCLUSION ==========================================

if (exclude_pattern == "") {
  filtered_genes <- heritable_genes
  cat("\nEXCLUDE_PATTERN empty: no hits dropped\n")
} else {
  excluded <- grepl(exclude_pattern, heritable_genes$GeneSymbol)
  filtered_genes <- heritable_genes[!excluded, ]
  cat("\nEXCLUDE_PATTERN '", exclude_pattern, "' dropped ", sum(excluded),
      " of ", nrow(heritable_genes), " hits (",
      length(unique(heritable_genes$GeneSymbol[excluded])),
      " distinct symbols)\n", sep = "")
}

# Columns the .Rmd's heritable_df carries. The percentile is over the filtered
# list, so it runs from 100*(1 - 1/n) down to 0.
filtered_genes <- filtered_genes[order(-filtered_genes$resid), ]
filtered_genes$resid_rank <- seq_len(nrow(filtered_genes))
filtered_genes$resid_percentile <- 100 *
  (1 - filtered_genes$resid_rank / nrow(filtered_genes))

cat("Filtered list:", nrow(filtered_genes), "rows,",
    length(unique(filtered_genes$GeneSymbol)), "distinct symbols\n")

write.csv(fit_data[order(-fit_data$resid), ],
          paste0(output.dir, "memoryseq_allfitdata.csv"), row.names = FALSE)
write.csv(heritable_genes,
          paste0(output.dir, "potential_heritable_genes.csv"), row.names = FALSE)
write.csv(filtered_genes[, c("GeneSymbol", "gene_id", "Mean_tpm", "CV_tpm",
                             "resid", "resid_rank", "resid_percentile")],
          paste0(output.dir, "potential_heritable_genes_filtered.csv"),
          row.names = FALSE)

cat("\nWrote", paste0(output.dir, "memoryseq_allfitdata.csv"), "\n")
cat("Wrote", paste0(output.dir, "potential_heritable_genes.csv"), "\n")
cat("Wrote", paste0(output.dir, "potential_heritable_genes_filtered.csv"),
    "-- the file Panel J reads\n")

cat("\nTop 20 of the filtered list by residual:\n")
print(filtered_genes %>%
        select(GeneSymbol, controls, Mean_tpm, CV_tpm, resid) %>%
        head(20))
