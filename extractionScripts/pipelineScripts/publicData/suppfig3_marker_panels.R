# Published marker panels and per-gene markers for each way of defining the
# infected population, against the same uninfected cells.
#
# The permutation null is computed from per-gene mean expression vectors rather
# than by subsetting the expression matrix inside each draw. The statistic is a
# mean over cells of a mean over genes, and both means are unweighted, so the
# order can be swapped: the arm's difference for any gene set is
# mean(arm_gene_means[set]) - mean(neg_gene_means[set]). Each arm therefore needs
# one rowMeans pass instead of one sparse row-subset per draw (row-subsetting a
# dgCMatrix is O(nnz)).
#
# Two things the null depends on:
#   - cut() drops the names of the expression vector, so the bin factor must be
#     named before it is indexed by gene symbol, or every draw comes back empty.
#   - the detection floor for the gene universe is 0.5%, not 5%: a 5% floor drops
#     most of the cytotoxic panel, which would leave the null testing a different
#     gene set than the module score. The matched draws come from genes of the
#     same detection rate.
#
# Each marker gene is tested on its own as well as inside a panel.

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(purrr); library(tibble)
})
# Configuration comes from the calling shell script, so one copy of this file
# serves any data tree. The defaults keep it runnable by hand.
get_env <- function(name, default) {
  v <- Sys.getenv(name, unset = NA)
  if (is.na(v) || v == "") default else v
}

# functions.R builds its default arguments from data.root, which config.R
# defines, so config.R has to be sourced first or sourcing functions.R fails
# with "object 'data.root' not found". config.R sits beside functions.R.
functions_R <- get_env("FUNCTIONS_R", "../../../plotScripts/functions.R")
source(file.path(dirname(functions_R), "config.R"))
source(functions_R)

rds.dir    <- get_env("RDS_DIR", "/projects/b1042/GoyalLab/egrody/extractedData/Figure1/singleCell/Seurat/RDS/")
table.dir  <- get_env("TABLE_DIR", "/projects/b1042/GoyalLab/egrody/extractedData/Figure1/combined/tiers/")
output.dir <- get_env("OUTPUT_DIR", "/projects/b1042/GoyalLab/egrody/plotData/Supplemental/SuppFig3/")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

samples       <- c("Invitro", "W2", "Pacute", "D195")
low_load_max  <- 8          # low load is 2-7 Ca*Re viral UMIs, the Panel C split
# A DNA+ population mixes silent, low- and high-transcribing cells, so the
# pooled Ca*Re arm is the closer analogue to it than the low-load arm alone.
arm_levels    <- c("all Ca*Re", "low (2-7)", "high (8+)", "10X-called")
min_cells     <- 10
n_random      <- 2000
expr_bins     <- 10
detection_min <- 0.005
set_seed      <- 1

reservoir_marker_sets <- list(
  Chomont_Tcm = c("CCR7", "SELL", "TCF7", "LEF1", "IL7R", "BACH2", "FOXP1",
                  "NOSIP", "NR3C2", "MYB", "KLF2"),
  Collora_cytotoxic = c("GZMB", "GZMA", "NKG7", "PRF1", "KLRD1", "CTSW", "CST7",
                        "CCL5", "HOPX", "KLRG1"),
  Weymar_latent_clone = c("HLA-DRA", "HLA-DRB1", "CD74", "CCL5", "GZMA", "GZMK",
                          "CST7", "LYAR", "DUSP2"),
  Sun_LN_survival = c("CD44", "CD28", "IL7R", "IL21R"),
  Sun_blood_checkpoint = c("PDCD1", "CTLA4", "TIGIT", "LAG3", "HAVCR2", "TOX"),
  Clark_death_receptor = c("TNFRSF10A", "TNFRSF10B", "FAS", "FASLG", "TRADD",
                           "FADD", "CASP8", "CASP10", "CASP3", "BID", "CFLAR",
                           "TNFRSF1A", "TNFRSF1B"),
  Clark_necroptosis = c("RIPK1", "RIPK3", "MLKL", "ZBP1", "CYLD", "TNFAIP3",
                        "BIRC2", "BIRC3", "XIAP", "CASP8", "CFLAR"))

marker_genes <- unique(c(unlist(reservoir_marker_sets), "SATB1", "GPHN",
                         "RASGRF2", "MKI67", "BIRC5", "TOX2", "FOXP3"))

tiers_all <- read_tier_tables(table.dir, samples)

results <- map(setNames(samples, samples), function(sample_name) {
  seurat_file <- paste0(rds.dir, sample_name, ".rds")
  if (!check_migrated(seurat_file)) { cat("missing object:", sample_name, "\n"); return(NULL) }
  tiers <- tiers_all[[sample_name]] %>% filter(tier %in% c("both", "CaRe_only", "GEX_only", "neg"))
  obj <- readRDS(seurat_file)
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj, assay = "RNA")
  tiers <- tiers %>% filter(cellID %in% colnames(obj))
  obj <- subset(obj, cells = tiers$cellID)

  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  universe <- names(which(Matrix::rowMeans(counts > 0) >= detection_min))
  data_layer <- GetAssayData(obj, assay = "RNA", layer = "data")[universe, ]
  expression <- Matrix::rowMeans(data_layer)
  expression_bin <- cut(expression,
                        unique(quantile(expression, seq(0, 1, length.out = expr_bins + 1))),
                        include.lowest = TRUE)
  # cut() returns an unnamed factor; the gene symbols are needed below
  names(expression_bin) <- universe

  arm_cellids <- list(
    `all Ca*Re`  = tiers$cellID[tiers$tier != "neg" & tiers$care_UMIs > 0],
    `low (2-7)`  = tiers$cellID[tiers$tier != "neg" & tiers$care_UMIs > 0 &
                                  tiers$care_UMIs < low_load_max],
    `high (8+)`  = tiers$cellID[tiers$tier != "neg" & tiers$care_UMIs >= low_load_max],
    `10X-called` = tiers$cellID[tiers$tier != "neg" & tiers$GEX_UMIs > 0])
  neg_cells <- tiers$cellID[tiers$tier == "neg"]
  arm_cellids <- arm_cellids[lengths(arm_cellids) >= min_cells]
  if (!length(arm_cellids) || length(neg_cells) < min_cells) {
    rm(obj, counts, data_layer); gc(); return(NULL)
  }

  # One pass per cell set: the per-gene mean expression across those cells
  neg_gene_means <- Matrix::rowMeans(data_layer[, neg_cells, drop = FALSE])
  arm_gene_means <- map(arm_cellids, function(test_cells) {
    Matrix::rowMeans(data_layer[, test_cells, drop = FALSE])
  })

  # Equivalent to the mean per-cell panel average, arm minus uninfected
  panel_difference <- function(genes, arm) {
    mean(arm_gene_means[[arm]][genes]) - mean(neg_gene_means[genes])
  }

  null_table <- map_dfr(names(arm_cellids), function(arm) {
    test_cells <- arm_cellids[[arm]]
    map_dfr(names(reservoir_marker_sets), function(set_name) {
      genes <- intersect(reservoir_marker_sets[[set_name]], universe)
      if (length(genes) < 2) {
        return(tibble(sample = sample_name, arm = arm, set = set_name,
                      n_cells = length(test_cells), n_genes = length(genes),
                      observed = NA_real_, null_mean = NA_real_, null_sd = NA_real_,
                      z = NA_real_, p_two_sided = NA_real_))
      }
      observed <- panel_difference(genes, arm)
      wanted <- table(expression_bin[genes])
      others <- setdiff(universe, genes)
      pools  <- split(others, expression_bin[others])
      set.seed(set_seed)
      draws <- replicate(n_random, {
        picked <- unlist(map2(pools[names(wanted)], as.integer(wanted), function(pool, n) {
          if (is.null(pool) || !length(pool) || n == 0) return(character(0))
          pool[sample.int(length(pool), min(length(pool), n))]
        }), use.names = FALSE)
        if (!length(picked)) return(NA_real_)
        panel_difference(picked, arm)
      })
      draws <- draws[!is.na(draws)]
      tibble(sample = sample_name, arm = arm, set = set_name,
             n_cells = length(test_cells), n_genes = length(genes),
             observed = observed, null_mean = mean(draws), null_sd = sd(draws),
             z = (observed - mean(draws)) / sd(draws),
             p_two_sided = (sum(abs(draws - mean(draws)) >= abs(observed - mean(draws))) + 1) /
               (length(draws) + 1))
    })
  })

  # Per-gene test, each arm against uninfected
  marker_present <- intersect(marker_genes, rownames(data_layer))
  marker_matrix <- as.matrix(data_layer[marker_present, , drop = FALSE])
  gene_table <- map_dfr(names(arm_cellids), function(arm) {
    test_cells <- arm_cellids[[arm]]
    map_dfr(marker_present, function(gene) {
      arm_values <- marker_matrix[gene, test_cells]
      neg_values <- marker_matrix[gene, neg_cells]
      tibble(sample = sample_name, arm = arm, gene = gene,
             n_cells = length(test_cells),
             pct_arm = 100 * mean(arm_values > 0), pct_neg = 100 * mean(neg_values > 0),
             log2FC = log2((mean(arm_values) + 0.01) / (mean(neg_values) + 0.01)),
             p_value = suppressWarnings(wilcox.test(arm_values, neg_values)$p.value))
    })
  }) %>% mutate(p_adj = p.adjust(p_value, "BH"))

  absent <- setdiff(marker_genes, rownames(data_layer))
  rm(obj, counts, data_layer, marker_matrix); gc()
  list(null = null_table, genes = gene_table, absent = absent,
       n_arms = lengths(arm_cellids), n_neg = length(neg_cells))
})
results <- results[!map_lgl(results, is.null)]

null_all <- map_dfr(results, "null")
gene_all <- map_dfr(results, "genes")
write.csv(null_all, paste0(output.dir, "SuppFig3_markerPanels_null.csv"), row.names = FALSE)
write.csv(gene_all, paste0(output.dir, "SuppFig3_markerPanels_perGene.csv"), row.names = FALSE)

cat("\n=== Cells per arm ===\n")
print(map(results, "n_arms"))
cat("uninfected per sample:\n"); print(map_int(results, "n_neg"))

cat("\n=== Marker genes absent from the reference ===\n")
print(map(results, "absent"))

cat("\n=== Expression-matched random-gene null, each arm vs uninfected ===\n")
print(as.data.frame(null_all %>%
  transmute(sample, arm, set, n_cells, n_genes, observed = round(observed, 4),
            null_mean = round(null_mean, 4), z = round(z, 2),
            p_two_sided = signif(p_two_sided, 2)) %>%
  arrange(sample, set, factor(arm, arm_levels))), row.names = FALSE)

cat("\n=== Per-gene, each arm vs uninfected (BACH2 and the survival markers first) ===\n")
print(as.data.frame(gene_all %>%
  filter(gene %in% c("BACH2", "SATB1", "CD28", "IL21R", "CD44", "IL7R", "CCR7",
                     "TCF7", "LEF1", "SELL", "GZMB", "GZMA", "CCL5", "MKI67")) %>%
  transmute(sample, gene, arm, log2FC = round(log2FC, 3), p_adj = signif(p_adj, 2)) %>%
  arrange(gene, sample, factor(arm, arm_levels))), row.names = FALSE)

cat("\n=== Per-gene, full table ===\n")
print(as.data.frame(gene_all %>%
  transmute(sample, gene, arm, pct_arm = round(pct_arm, 1), pct_neg = round(pct_neg, 1),
            log2FC = round(log2FC, 3), p_adj = signif(p_adj, 2)) %>%
  arrange(gene, sample, factor(arm, arm_levels))), row.names = FALSE)

cat("\nWrote tables to", output.dir, "\n")
