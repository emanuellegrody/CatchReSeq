# splice+ vs splice- among infected CD4 T cells, in the viremic public datasets
# and the untreated SIV samples. Viral burden is the dominant confounder (a
# splice call needs a spliced molecule, so splice+ cells carry far more viral
# UMIs), so every test either adjusts for it or matches on it.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })
# Resolves plotScripts/ from this script's own location, so the repo can be
# cloned anywhere. Rscript exposes that location only through its command line.
script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) dirname(normalizePath(sub("^--file=", "", hit[1]))) else getwd()
})
plot_dir <- normalizePath(file.path(script_dir, "..", "..", "..", "plotScripts"),
                          mustWork = FALSE)

source(file.path(plot_dir, "packages.R"))
source(file.path(plot_dir, "config.R"))
source(file.path(plot_dir, "functions.R"))
# Caches the Supp Fig 3 Panel D chunk reads, so they live under extractedData
# rather than beside this script.
spliceviremic.dir <- paste0(data.root, "Figure3/spliceViremic/")
dir.create(spliceviremic.dir, recursive = TRUE, showWarnings = FALSE)
plot.dir    <- paste0(plot.root, "Figure3/Giron/")
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
nperm <- 10000
set.seed(1)

extracted <- readRDS(paste0(spliceviremic.dir, "splice_viremic_cells.rds"))
cells <- extracted$cells
siv_splice_genes <- extracted$siv_splice_genes

score_names <- c(names(cd4t_subset_markers), "OXIDATIVE_PHOSPHORYLATION", "GLYCOLYSIS",
                 "MTORC1_SIGNALING", "INTERFERON_GAMMA_RESPONSE", "INTERFERON_ALPHA_RESPONSE",
                 "TNFA_SIGNALING_VIA_NFKB", "INFLAMMATORY_RESPONSE", "PUBLIC_HIV_UP",
                 "SIV_SPLICE_UP")
score_names <- intersect(score_names, colnames(cells))
score_labels <- c(naive_Tcm = "naive/Tcm", Tem = "Tem", cytotoxic = "cytotoxic", Treg = "Treg",
                  cycling = "cycling", OXIDATIVE_PHOSPHORYLATION = "OxPhos",
                  GLYCOLYSIS = "Glycolysis", MTORC1_SIGNALING = "mTORC1",
                  INTERFERON_GAMMA_RESPONSE = "IFN-gamma", INTERFERON_ALPHA_RESPONSE = "IFN-alpha",
                  TNFA_SIGNALING_VIA_NFKB = "TNFA/NFkB", INFLAMMATORY_RESPONSE = "Inflammatory",
                  PUBLIC_HIV_UP = "viremic HIV+ signature", SIV_SPLICE_UP = "SIV splice+ signature")
dataset_levels <- c("Collora_Viremic", "HIVseq_Viremic", "Invitro", "Pacute", "W2")

# Burden bins used for matching and as a stratum in the permutation null
burden_breaks <- c(0, 2, 7, 15, 30, 100, Inf)
burden_labels <- c("1-2", "3-7", "8-15", "16-30", "31-100", ">100")
cells <- cells %>%
  mutate(dataset = factor(dataset, levels = dataset_levels),
         burden_bin = cut(viral_UMIs, burden_breaks, labels = burden_labels),
         stratum = paste(library, burden_bin))

cat("=== infected cells, splice+ and viral burden ===\n")
print(as.data.frame(cells %>% group_by(dataset) %>%
  summarize(infected = n(), splice_pos = sum(spliced),
            median_UMIs_pos = median(viral_UMIs[spliced]),
            median_UMIs_neg = median(viral_UMIs[!spliced]),
            p_burden = suppressWarnings(wilcox.test(viral_UMIs ~ spliced)$p.value),
            .groups = "drop") %>%
  mutate(p_burden = signif(p_burden, 2))), row.names = FALSE)

cat("\n=== strata where both splice+ and splice- occur (the matched comparison) ===\n")
usable <- cells %>% group_by(dataset, stratum) %>%
  filter(any(spliced), any(!spliced)) %>% ungroup()
print(as.data.frame(usable %>% group_by(dataset) %>%
  summarize(strata = n_distinct(stratum), splice_pos = sum(spliced),
            splice_neg = sum(!spliced), .groups = "drop")), row.names = FALSE)

# --- test 1: burden-adjusted linear model, per dataset -------------------------
# score ~ spliced + log2(viral UMIs) + log2(host UMIs) (+ library where there is
# more than one). p comes from permuting the splice labels within library x
# burden bin, so the null keeps the burden difference between the groups.
# Each score is residualised once on viral burden, host size and library, none
# of which involves the splice labels; the statistic is then the difference in
# mean residual between splice+ and splice- cells.
residualise <- function(df, score) {
  covariates <- c("log2(viral_UMIs)", "log2(nCount)")
  if (n_distinct(df$library) > 1) covariates <- c(covariates, "library")
  residuals(lm(reformulate(covariates, score), data = df))
}
permute_within <- function(spliced, stratum) {
  ave(spliced, stratum, FUN = sample)
}

adjusted_tests <- purrr::map_dfr(dataset_levels, function(ds) {
  df <- cells %>% filter(dataset == ds)
  if (sum(df$spliced) < 3 || sum(!df$spliced) < 3) return(NULL)
  purrr::map_dfr(score_names, function(score) {
    r <- residualise(df, score)
    effect_of <- function(spliced) mean(r[spliced]) - mean(r[!spliced])
    observed <- effect_of(df$spliced)
    null <- replicate(nperm, effect_of(permute_within(df$spliced, df$stratum)))
    tibble(dataset = ds, score = score, n_pos = sum(df$spliced), n_neg = sum(!df$spliced),
           effect = observed, null_sd = sd(null),
           z = (observed - mean(null)) / sd(null),
           p = min(1, 2 * min(mean(null >= observed), mean(null <= observed))))
  }) %>% mutate(p_adj = p.adjust(p, "BH"))
})
cat("\n=== splice+ vs splice-, burden- and size-adjusted (BH across scores within dataset) ===\n")
print(as.data.frame(adjusted_tests %>%
  mutate(score = score_labels[score], across(c(effect, z), ~ round(.x, 3)),
         across(c(p, p_adj), ~ signif(.x, 2)))), row.names = FALSE)

# --- test 2: do the datasets agree on direction? -------------------------------
# Spearman of the per-score z between every pair of datasets, and the sign
# agreement of the public datasets with each SIV sample.
z_matrix <- adjusted_tests %>% select(dataset, score, z) %>%
  pivot_wider(names_from = dataset, values_from = z) %>%
  tibble::column_to_rownames("score") %>% as.matrix()
rho_matrix <- cor(z_matrix, method = "spearman")
sign_matrix <- crossprod(sign(z_matrix) > 0) + crossprod(sign(z_matrix) < 0)
pair_rho <- purrr::map_dfr(seq_len(ncol(z_matrix) - 1), function(i)
  purrr::map_dfr((i + 1):ncol(z_matrix), function(j)
    tibble(pair = paste(colnames(z_matrix)[c(i, j)], collapse = " vs "),
           rho = rho_matrix[i, j],
           same_sign = sign_matrix[i, j] / nrow(z_matrix))))
cat("\n=== agreement of the per-score effects between datasets ===\n")
print(as.data.frame(pair_rho %>% mutate(across(where(is.double), ~ round(.x, 2)))), row.names = FALSE)

# --- test 3: the SIV Invitro splice+ genes in the human viremic cells ----------
# Each gene's normalized expression, splice+ vs splice-, matched on library x
# burden bin: the observed difference in means against the same permutation null.
gene_columns <- grep("^gene_", colnames(cells), value = TRUE)
gene_tests <- purrr::map_dfr(c("Collora_Viremic", "HIVseq_Viremic"), function(ds) {
  df <- cells %>% filter(dataset == ds)
  purrr::map_dfr(gene_columns, function(column) {
    if (all(df[[column]] == 0)) return(NULL)
    x <- df[[column]]
    effect_of <- function(spliced) mean(x[spliced]) - mean(x[!spliced])
    observed <- effect_of(df$spliced)
    null <- replicate(nperm, effect_of(permute_within(df$spliced, df$stratum)))
    tibble(dataset = ds, gene = sub("^gene_", "", column), difference = observed,
           z = (observed - mean(null)) / sd(null),
           p = min(1, 2 * min(mean(null >= observed), mean(null <= observed))))
  }) %>% mutate(p_adj = p.adjust(p, "BH"))
})
cat("\n=== SIV Invitro splice+ genes, tested in the human viremic cells ===\n")
print(as.data.frame(gene_tests %>% mutate(across(c(difference, z), ~ round(.x, 3)),
                                          across(c(p, p_adj), ~ signif(.x, 2))) %>%
  arrange(gene, dataset)), row.names = FALSE)
cat("\nSIV splice+ genes up in human viremic splice+ cells:",
    sum(gene_tests$z > 0), "of", nrow(gene_tests), "dataset-gene pairs\n")

saveRDS(list(adjusted = adjusted_tests, pairs = pair_rho, genes = gene_tests),
        paste0(spliceviremic.dir, "splice_viremic_results.rds"))

# --- plots ----------------------------------------------------------------------
heat <- adjusted_tests %>%
  mutate(score = factor(score_labels[score], levels = rev(score_labels[score_names])),
         dataset = factor(dataset, levels = dataset_levels),
         label = sprintf("%.2f%s", z, ifelse(p_adj < 0.05, "\n*",
                                      ifelse(p < 0.05, "\n(p<.05)", ""))))
plot_heat <- ggplot(heat, aes(dataset, score, fill = pmax(pmin(z, 4), -4))) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 2.4, lineheight = 0.85) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-4, 4),
                       name = "z, splice+ vs splice-\n(burden + host size adjusted)") +
  labs(x = NULL, y = NULL,
       subtitle = sprintf("splice+ / splice- cells: %s",
                          paste(sprintf("%s %d/%d", adjusted_tests$dataset[!duplicated(adjusted_tests$dataset)],
                                        adjusted_tests$n_pos[!duplicated(adjusted_tests$dataset)],
                                        adjusted_tests$n_neg[!duplicated(adjusted_tests$dataset)]),
                                collapse = ", "))) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))

burden_df <- cells %>% mutate(splice = ifelse(spliced, "splice+", "splice-"))
plot_burden <- ggplot(burden_df, aes(splice, viral_UMIs, color = splice)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "grey60", fill = NA) +
  geom_point(size = 0.8, alpha = 0.5,
             position = position_jitter(width = 0.2, height = 0, seed = 1)) +
  facet_wrap(~ dataset, nrow = 1) +
  scale_color_manual(values = c(`splice+` = "purple", `splice-` = "orange"), guide = "none") +
  scale_y_log10() +
  labs(x = NULL, y = "viral UMIs per cell",
       subtitle = "the confounder every test below adjusts for") +
  theme_classic(base_size = 8) + theme(strip.background = element_blank())

gene_plot_df <- gene_tests %>%
  mutate(gene = factor(gene, levels = sort(unique(gene))),
         dataset = factor(dataset, levels = dataset_levels))
plot_genes <- ggplot(gene_plot_df, aes(z, gene, color = dataset)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c(Collora_Viremic = "#6a3d9a", HIVseq_Viremic = "#33a02c"),
                     name = NULL) +
  labs(x = "z, splice+ vs splice- (matched permutation)", y = NULL,
       subtitle = "SIV Invitro splice+ genes, tested in human viremic cells") +
  theme_classic(base_size = 8) + theme(legend.position = "bottom")

ggsave(paste0(plot.dir, "Diag_spliceViremic_scores.pdf"),
       plot_burden / (plot_heat + plot_genes + plot_layout(widths = c(1.3, 1))) +
         plot_layout(heights = c(1, 2.2)), width = 11, height = 9)
