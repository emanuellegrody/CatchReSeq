#!/usr/bin/env Rscript
# memoryseq_qc.R
# Library-level QC on the LD002 clone bulk RNA-seq TPM matrix, run before
# memoryseq_heritable_genes.R reads the same matrix.
#
# The heritable-gene call is a between-library variance test, so it inherits
# every between-library artefact: a clone sequenced shallower than the rest, or
# normalised onto a different scale, raises the CV of the genes it affects and
# those genes come out as hits. The four views written here are what rules that
# out -- the TPM distributions should be superimposable, the clone-to-clone
# Spearman correlations should be tight and comparable, and PC1 should not
# separate libraries by depth.
#
# Required env vars:
#   TPM_MATRIX        gene-level TPM table from nf-core/rnaseq
#                     (star_salmon/salmon.merged.gene_tpm.tsv)
#   PLOT_DIR          dir to write the QC PDFs to
#
# Optional env vars:
#   SAMPLE_PATTERN    regex selecting the library columns (default "^Clone_")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
})

get_env_required <- function(name) {
  value <- Sys.getenv(name, unset = NA)
  if (is.na(value) || value == "") {
    stop("Required environment variable not set: ", name)
  }
  value
}

tpm_matrix_path <- get_env_required("TPM_MATRIX")
plot.dir        <- get_env_required("PLOT_DIR")
sample_pattern  <- Sys.getenv("SAMPLE_PATTERN", unset = "^Clone_")

if (!grepl("/$", plot.dir)) plot.dir <- paste0(plot.dir, "/")
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(tpm_matrix_path)) {
  stop("TPM matrix not found: ", tpm_matrix_path)
}

# ========== LOAD ==============================================================

tpm_raw <- read.delim(tpm_matrix_path, check.names = FALSE,
                      stringsAsFactors = FALSE)
library_columns <- grep(sample_pattern,
                        setdiff(colnames(tpm_raw), c("gene_id", "gene_name")),
                        value = TRUE)
if (length(library_columns) == 0) {
  stop("SAMPLE_PATTERN '", sample_pattern, "' matched no columns in ",
       tpm_matrix_path)
}

# Same collapse the heritable-gene call uses, so both read one matrix.
tpm_by_symbol <- tpm_raw %>%
  select(gene_name, all_of(library_columns)) %>%
  group_by(gene_name) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  filter(gene_name != "", !is.na(gene_name)) %>%
  column_to_rownames("gene_name")
tpm_by_symbol <- tpm_by_symbol[rowSums(tpm_by_symbol) > 0, ]

cat("Libraries:", length(library_columns), "| genes with any expression:",
    nrow(tpm_by_symbol), "\n")

tpm_long <- tpm_by_symbol %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sampleID", values_to = "tpm") %>%
  mutate(log1p_tpm = log1p(tpm))

# ========== 1. TPM DISTRIBUTION PER LIBRARY ===================================

cat("\nPer-library log1p(TPM) summary:\n")
print(tpm_long %>%
        group_by(sampleID) %>%
        summarise(n      = n(),
                  mean   = mean(log1p_tpm),
                  median = median(log1p_tpm),
                  sd     = sd(log1p_tpm),
                  q25    = quantile(log1p_tpm, 0.25),
                  q75    = quantile(log1p_tpm, 0.75),
                  .groups = "drop"),
      n = Inf)

plot_distribution <- ggplot(tpm_long, aes(x = sampleID, y = log1p_tpm)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7, outlier.size = 0.2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()) +
  labs(x = "Library", y = "log1p(TPM)",
       title = "Expression distribution per clone library")
ggsave(paste0(plot.dir, "QC_tpm_distribution.pdf"), plot_distribution,
       device = "pdf", width = 0.4 * length(library_columns) + 2, height = 4)

# ========== 2. MEAN-CV TREND ==================================================

# The trend memoryseq_heritable_genes.R fits its GLM to. A loess that bends
# sharply at the low-expression end is the expression floor doing its job;
# what matters here is that no single library drives the spread.
gene_stats <- tpm_long %>%
  group_by(gene) %>%
  summarise(mean_tpm = mean(tpm),
            cv       = sd(tpm) / mean(tpm),
            .groups  = "drop") %>%
  filter(is.finite(cv), mean_tpm > 0)

plot_cv_mean <- ggplot(gene_stats, aes(x = log2(mean_tpm), y = cv)) +
  geom_point(alpha = 0.1, size = 0.3) +
  geom_smooth(method = "loess", formula = y ~ x, color = "red") +
  theme_classic() +
  labs(x = "log2(mean TPM)", y = "Coefficient of variation",
       title = "Coefficient of variation vs mean expression")
ggsave(paste0(plot.dir, "QC_cv_vs_mean.pdf"), plot_cv_mean,
       device = "pdf", width = 5, height = 4)

cat("\nCoefficient of variation across genes:\n")
print(summary(gene_stats$cv))
for (threshold in c(0.2, 0.5, 1.0)) {
  cat(sprintf("Genes with CV > %g: %.1f%%\n",
              threshold, 100 * mean(gene_stats$cv > threshold)))
}

# ========== 3. LIBRARY CORRELATIONS ===========================================

correlation_matrix <- cor(tpm_by_symbol, method = "spearman")
off_diagonal <- correlation_matrix[upper.tri(correlation_matrix)]

cat(sprintf("\nClone-to-clone Spearman correlation: %.3f to %.3f (median %.3f)\n",
            min(off_diagonal), max(off_diagonal), median(off_diagonal)))

# The colour scale is set to the observed range rather than 0-1, because every
# pair sits in a narrow band and a full-range scale renders the heatmap flat.
pdf(paste0(plot.dir, "QC_library_correlations.pdf"),
    width = 0.5 * length(library_columns) + 3,
    height = 0.5 * length(library_columns) + 2)
pheatmap(correlation_matrix,
         color = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100),
         breaks = seq(min(off_diagonal), max(off_diagonal), length.out = 101),
         display_numbers = TRUE,
         number_format = "%.3f",
         fontsize_number = 7,
         main = sprintf("Spearman correlation (scale %.3f-%.3f)",
                        min(off_diagonal), max(off_diagonal)))
dev.off()

# ========== 4. PCA ============================================================

pca_result    <- prcomp(t(log1p(tpm_by_symbol)), scale. = TRUE)
var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)
pca_dataframe <- as.data.frame(pca_result$x) %>%
  rownames_to_column("sampleID") %>%
  left_join(tpm_long %>%
              group_by(sampleID) %>%
              summarise(detected = sum(tpm > 0), .groups = "drop"),
            by = "sampleID")

plot_pca <- ggplot(pca_dataframe, aes(x = PC1, y = PC2, label = sampleID)) +
  geom_point(aes(color = detected), size = 3) +
  geom_text(vjust = -1, size = 2.5) +
  scale_color_viridis_c(name = "Genes\ndetected") +
  theme_classic() +
  labs(title = "PCA of clone libraries, coloured by genes detected",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)"))
ggsave(paste0(plot.dir, "QC_pca.pdf"), plot_pca,
       device = "pdf", width = 6, height = 5)

cat("\nVariance explained by the first PCs:", head(var_explained), "\n")
# PC1 tracking genes-detected is the depth artefact this view exists to catch.
# Libraries that all detect the same number of genes make the correlation
# undefined rather than zero, which is a pass, not a failure.
if (length(unique(pca_dataframe$detected)) > 1) {
  cat("Spearman correlation of PC1 with genes detected:",
      round(cor(pca_dataframe$PC1, pca_dataframe$detected,
                method = "spearman"), 3), "\n")
} else {
  cat("Spearman correlation of PC1 with genes detected: not defined --",
      "every library detects", unique(pca_dataframe$detected), "genes\n")
}

cat("\nWrote 4 QC PDFs to", plot.dir, "\n")
