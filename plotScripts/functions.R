# functions.R  (plotScripts)
# The analysis and plotting library the Figure*.Rmd and SuppFigs.Rmd notebooks
# source. Holds the detection-threshold helpers (set_minimums), the CD4 T gate
# (apply_cd4t_gate, cd4t_gate_table), the Seurat wrapper (SeuratPipeline) and
# every plotting function the figures call.
#
# Expects plotScripts/config.R to have been sourced first: data.root appears in
# the default arguments below, and care.scripts.dir locates the pipeline library.

# =============================================================================
# CaRe and TenX
# =============================================================================

# stepThree/stepFour helpers (extract_category, process_bamsort, set_minimums and
# friends) are defined in the pipeline library, not here, so a SLURM worker can
# source them without config.R or the plotting stack. Edit them there.
source(file.path(care.scripts.dir, "..", "extractionScripts", "pipelineScripts",
                 "bamsort", "bamsort_functions.R"))

## SIV splice+ definition
siv_splice_cols <- function(x) {
  cols <- if (is.data.frame(x)) colnames(x) else x
  c(intersect(c("SS", "MS"), cols), grep("^SS\\.MS", cols, value = TRUE))
}
is_siv_splice_category <- function(categories) {
  cats <- make.names(categories)
  cats %in% c("SS", "MS") | grepl("^SS\\.MS", cats)
}

## Splicing
classify_site <- function(value, reference_values, reference_sites) {
  if (is.na(value)) {
    return("nc")
  }
  
  distances <- abs(value - reference_values)
  valid_distances <- distances[!is.na(distances)]
  
  if (length(valid_distances) == 0) {
    return("nc")
  }
  
  min_dist <- min(valid_distances)
  closest_site <- reference_sites[which.min(distances)]
  
  if (min_dist <= 3) {
    return(as.character(closest_site))
  } else if (min_dist <= 12) {
    return("mnc")
  } else {
    return("nc")
  }
}
## Combined TenX and CaRe
resolve_multimap_gs <- function(df) {
  # If no merge conflict occurred, category columns won't have .x/.y suffixes
  if (!"category.x" %in% names(df)) {
    if ("category" %in% names(df)) {
      cat("No multi-mapping detected (single source). Returning as-is.\n")
      return(df[, c("cellID", "UMI", "category"), drop = FALSE])
    } else {
      stop("No 'category' or 'category.x'/'category.y' columns found in input.")
    }
  }
  # Each category names the set of isoforms a read is consistent with, so two
  # calls on one molecule resolve by intersecting those sets: SS-MS is {SS, MS}
  # and US-SS is {US, SS}, so together they mean SS. An empty intersection is a
  # real disagreement, and there a category naming a single isoform (US/SS/MS)
  # takes precedence; if that does not settle it the call falls back to "any".
  # Banded SS-MS calls collapse first, so SS-MS4 against SS-MS7 is agreement.
  priority_cats <- c("US", "SS", "MS")
  cat_isoforms  <- list("US" = "US", "SS" = "SS", "MS" = "MS",
                        "SS-MS" = c("SS", "MS"), "US-SS" = c("US", "SS"),
                        "any" = c("US", "SS", "MS"))
  isoforms_to_cat <- function(iso) {
    key <- paste(sort(iso), collapse = "|")
    hit <- names(cat_isoforms)[vapply(cat_isoforms,
             function(x) identical(paste(sort(x), collapse = "|"), key), logical(1))]
    if (length(hit)) hit[1] else NA_character_
  }
  normalize_cat <- function(x) ifelse(is.na(x), x,
                                      ifelse(grepl("^SS[-.]MS", x), "SS-MS", x))
  df$category.x <- normalize_cat(df$category.x)
  df$category.y <- normalize_cat(df$category.y)

  # Initialize counters
  stats <- list(
    total = 0,
    resolved_same = 0,
    resolved_priority = 0,
    resolved_by_intersection = 0,
    resolved_to_any = 0,
    no_conflict = 0
  )
  
  # Process each row
  result <- vector("list", nrow(df))
  
  for (i in seq_len(nrow(df))) {
    cat_x <- df$category.x[i]
    cat_y <- df$category.y[i]
    # read_x <- ifelse(is.na(df$read_count.x[i]), 0, df$read_count.x[i])
    # read_y <- ifelse(is.na(df$read_count.y[i]), 0, df$read_count.y[i])
    
    # Case 1: Both categories are NA - skip this row
    if (is.na(cat_x) && is.na(cat_y)) {
      result[[i]] <- NULL
      next
    }
    
    # Case 2: One category is NA - use the non-NA category
    if (is.na(cat_x) && !is.na(cat_y)) {
      result[[i]] <- data.frame(
        cellID = df$cellID[i],
        UMI = df$UMI[i],
        category = cat_y,
        #read_count = read_x + read_y,
        stringsAsFactors = FALSE
      )
      stats$no_conflict <- stats$no_conflict + 1
      next
    }
    
    if (!is.na(cat_x) && is.na(cat_y)) {
      result[[i]] <- data.frame(
        cellID = df$cellID[i],
        UMI = df$UMI[i],
        category = cat_x,
        #read_count = read_x + read_y,
        stringsAsFactors = FALSE
      )
      stats$no_conflict <- stats$no_conflict + 1
      next
    }
    
    # Case 3: Both categories present - apply resolution rules
    stats$total <- stats$total + 1
    
    # Rule 0: If both categories are identical, keep that category
    if (cat_x == cat_y) {
      result[[i]] <- data.frame(
        cellID = df$cellID[i],
        UMI = df$UMI[i],
        category = cat_x,
        #read_count = read_x + read_y,
        stringsAsFactors = FALSE
      )
      stats$resolved_same <- stats$resolved_same + 1
      next
    }
    
    # Rule 1: intersect the isoform sets both calls allow; a non-empty
    # intersection is the narrower call both reads support.
    iso_x <- cat_isoforms[[cat_x]]
    iso_y <- cat_isoforms[[cat_y]]
    if (!is.null(iso_x) && !is.null(iso_y)) {
      shared <- intersect(iso_x, iso_y)
      if (length(shared) > 0) {
        narrowed <- isoforms_to_cat(shared)
        if (!is.na(narrowed)) {
          result[[i]] <- data.frame(
            cellID = df$cellID[i],
            UMI = df$UMI[i],
            category = narrowed,
            stringsAsFactors = FALSE
          )
          stats$resolved_by_intersection <- stats$resolved_by_intersection + 1
          next
        }
      }
    }

    # Rule 2: the calls contradict; a single-isoform category outranks one that
# is not.
    if (cat_x %in% priority_cats && !(cat_y %in% priority_cats)) {
      result[[i]] <- data.frame(
        cellID = df$cellID[i],
        UMI = df$UMI[i],
        category = cat_x,
        #read_count = read_x + read_y,
        stringsAsFactors = FALSE
      )
      stats$resolved_priority <- stats$resolved_priority + 1
      next
    }
    
    if (cat_y %in% priority_cats && !(cat_x %in% priority_cats)) {
      result[[i]] <- data.frame(
        cellID = df$cellID[i],
        UMI = df$UMI[i],
        category = cat_y,
        #read_count = read_x + read_y,
        stringsAsFactors = FALSE
      )
      stats$resolved_priority <- stats$resolved_priority + 1
      next
    }
    
    # Rule 3: two tier 1 calls that disagree, or two unresolved calls
    result[[i]] <- data.frame(
      cellID = df$cellID[i],
      UMI = df$UMI[i],
      category = "any",
      #read_count = read_x + read_y,
      stringsAsFactors = FALSE
    )
    stats$resolved_to_any <- stats$resolved_to_any + 1
  }
  
  # Combine results, excluding NULLs
  result <- do.call(rbind, result[!sapply(result, is.null)])
  rownames(result) <- NULL
  
  # Print statistics
  cat("Multi-mapped UMIs (both categories present):", stats$total, "\n")
  cat("UMIs with identical categories:", stats$resolved_same, 
      sprintf("(%.1f%%)", ifelse(stats$total > 0, stats$resolved_same/stats$total*100, 0)), "\n")
  cat("UMIs resolved by priority category (US/SS/MS):", stats$resolved_priority, 
      sprintf("(%.1f%%)", ifelse(stats$total > 0, stats$resolved_priority/stats$total*100, 0)), "\n")
  cat("UMIs narrowed by intersecting isoform sets:", stats$resolved_by_intersection,
      sprintf("(%.1f%%)", ifelse(stats$total > 0, stats$resolved_by_intersection/stats$total*100, 0)), "\n")
  cat("UMIs resolved to 'any':", stats$resolved_to_any, 
      sprintf("(%.1f%%)", ifelse(stats$total > 0, stats$resolved_to_any/stats$total*100, 0)), "\n")
  cat("UMIs with no conflict (one or both categories NA):", stats$no_conflict, "\n")
  
  
  return(result)
}

# =============================================================================
# polyA internal vs terminal priming
# =============================================================================
## A parameter set that leaves most reads unassigned to any priming site
## measures only the assignable minority, so internal:terminal statements are
## restricted to sets where at least half the reads are assigned. The sweep
## CSVs written by the bam_polya_*.py scripts keep every set, and the stacked
## bar panels show the unassigned fraction as its own category; this cap applies
## to the ratio statements drawn from those sweeps.
MAX_UNASSIGNED_FRAC <- 0.5

assignable_sets <- function(df,
                            max_unassigned  = MAX_UNASSIGNED_FRAC,
                            terminal_col    = "explained_terminal",
                            internal_col    = "explained_internal",
                            unexplained_col = "unexplained",
                            label           = "") {
  missing_cols <- setdiff(c(terminal_col, internal_col, unexplained_col),
                          colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("assignable_sets: missing columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }
  total <- df[[terminal_col]] + df[[internal_col]] + df[[unexplained_col]]
  unassigned_frac <- ifelse(total > 0, df[[unexplained_col]] / total, 1)
  keep <- unassigned_frac <= max_unassigned
  if (nzchar(label)) {
    cat(sprintf(paste0("[%s] unassigned <= %.0f%% cap: %d of %d parameter ",
                       "sets kept (%d dropped)\n"),
                label, 100 * max_unassigned, sum(keep), length(keep),
                sum(!keep)))
  }
  df[keep, , drop = FALSE]
}

## The panels show every parameter set tested; the ones the cap excludes from
## the ratio statements are named with a strikethrough, so a reader can see both
## the whole sweep and which part of it the statements rest on. The strike is a
## combining overlay character, which the base pdf() device cannot encode, so
## panels using these labels are written with device = cairo_pdf.
strike_text <- function(x) {
  vapply(x, function(s) paste0(strsplit(s, "")[[1]], "\u0336", collapse = ""),
         character(1), USE.NAMES = FALSE)
}

cap_struck_labels <- function(df,
                              set_col        = "parameters",
                              pct_col        = "unexplained_pct",
                              max_unassigned = MAX_UNASSIGNED_FRAC) {
  missing_cols <- setdiff(c(set_col, pct_col), colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("cap_struck_labels: missing columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }
  sets <- as.character(df[[set_col]])
  struck <- df[[pct_col]] > 100 * max_unassigned
  labels <- sets
  labels[struck] <- strike_text(sets[struck])
  stats::setNames(labels, sets)
}

# =============================================================================
# Seurat
# =============================================================================
# Number of leading dimensions of a linear reduction (PCA, or harmony for
# integrated objects) used everywhere downstream: neighbor graphs, clustering,
# UMAP, and every distance / kNN analysis built on an embedding. RunPCA still
# computes its default 50 components, but analyses take only these leading ones,
# so a distance in one figure means the same thing as a distance in another.
REDUCTION_DIMENSIONS <- 30

## standard processing of counts matrices with Seurat
# QC thresholds are arguments so a dataset can be filtered on the thresholds its
# own source published; the upper caps accept Inf to drop them entirely.
SeuratPipeline <- function(file_location, sample_name,
                           output_dir = "/projects/b1042/GoyalLab/egrody/extractedData/",
                           plots = FALSE, rds = TRUE, species = "macaque", remove_mac = TRUE,
                           min_features = 200, max_features = 3500,
                           max_counts = 20000, max_mito = 5) {
  # scCustomize and the SingleCellExperiment/zellkonverter pair are optional
  # groups in packages.R, so they are attached here rather than at startup.
  # load_packages() installs if needed and attaches quietly; bare require()
  # announced each package and its whole Depends chain into the chunk output.
  load_packages("scCustomize", "SingleCellExperiment", "zellkonverter")

  plots_dir <- file.path(output_dir, "plots/")
  rds_dir <- file.path(output_dir, "RDS/")
  sapply(c(plots_dir, rds_dir), dir.create, recursive = TRUE, showWarnings = FALSE)
  # Reading in feature matrices. file_location is either a cellranger matrix
  # directory or an already-loaded counts matrix (used when the matrix has to be
  # split by feature type first, e.g. Gene Expression vs Antibody Capture).
  data <- if (is.character(file_location)) Read10X(data.dir = file_location) else file_location
  # Initialize the Seurat object with the raw (non-normalized data)
  sample <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 3,
                               min.features = min_features)
  # QC
  sample <- suppressMessages(Add_Mito_Ribo(sample, species = species))
  if (plots) {
    vplot <- VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"),
                     ncol = 3, layer = "counts")
    ggsave(vplot, file = paste0(plots_dir, sample_name, "_preQC_violinplot.svg"),
           width = 7, height = 7)
  }
  # Written as a logical rather than subset(subset = ) so an Inf cap is a no-op.
  keep <- sample$nFeature_RNA > min_features &
          sample$nFeature_RNA < max_features &
          sample$nCount_RNA   < max_counts &
          sample$percent_mito < max_mito
  sample <- subset(sample, cells = colnames(sample)[keep])
  if (plots) {
    vplot <- VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"),
                     ncol = 3, layer = "counts")
    ggsave(vplot, file = paste0(plots_dir, sample_name, "_postQC_violinplot.svg"),
           width = 7, height = 7)
  }
  # Remove CD4/CD8 doublets
  if ("CD4" %in% rownames(sample) && ("CD8A" %in% rownames(sample) || "CD8B" %in% rownames(sample))) {
    cd4_expr <- GetAssayData(sample, layer = "counts")["CD4", ]
    cd8a_expr <- if ("CD8A" %in% rownames(sample)) GetAssayData(sample, layer = "counts")["CD8A", ] else rep(0, ncol(sample))
    cd8b_expr <- if ("CD8B" %in% rownames(sample)) GetAssayData(sample, layer = "counts")["CD8B", ] else rep(0, ncol(sample))
    
    cd4_pos <- cd4_expr > 0
    cd8_pos <- (cd8a_expr > 0) | (cd8b_expr > 0)
    doublets <- cd4_pos & cd8_pos
    
    sample <- subset(sample, cells = colnames(sample)[!doublets])
  }
  # Remove viral spike-in genes (mac239, HXB2) so the human/host transcriptome
  # alone drives QC, clustering, and integration; viral counts are joined
  # separately (from bamsort) as metadata downstream.
  if (remove_mac == TRUE) {
    viral_genes <- intersect(c("mac239", "HXB2"), rownames(sample))
    if (length(viral_genes) > 0) {
      sample <- sample[!rownames(sample) %in% viral_genes, ]
    }
  }
  # Normalizing
  sample <- NormalizeData(sample, verbose = FALSE)
  # Feature selection
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  # Scaling
  genes <- rownames(sample)
  sample <- ScaleData(sample, features = genes, verbose = FALSE)
  # Linear dimension reduction
  sample <- RunPCA(sample, features = VariableFeatures(object = sample), verbose = FALSE)
  if (plots) {
    eplot <- ElbowPlot(sample, ndims = 50)
    ggsave(eplot, file = paste0(plots_dir, sample_name, "_elbow.svg"),
           width = 7, height = 7)
  }
  # Clustering
  sample <- FindNeighbors(sample, dims = 1:REDUCTION_DIMENSIONS, verbose = FALSE)
  sample <- FindClusters(sample, resolution = 0.3, algorithm = 4, verbose = FALSE)
  # Dimension reduction
  sample <- RunUMAP(sample, dims = 1:REDUCTION_DIMENSIONS, verbose = FALSE)#, umap.method = "umap-learn", metric = "correlation")
  if (plots) {
    dplot <- DimPlot(sample, reduction = "umap")
    ggsave(dplot, file = paste0(plots_dir, sample_name, "_umap.pdf"),
           width = 7, height = 7)
  }
  if (rds) {
    saveRDS(sample, file = paste0(rds_dir, sample_name, ".rds"))
    sce <- SingleCellExperiment(
      assays = list(counts = GetAssayData(sample, layer = "counts")),
      colData = sample@meta.data
    )
    suppressMessages(
      writeH5AD(sce, file = paste0(rds_dir, sample_name, ".h5ad"))
    )
  }
  
  # Return the variable
  return(sample)
}

# =============================================================================
# Joint Analyses
# =============================================================================
## joining and printing numbers
read_sample_files <- function(samples, directory, filtered = TRUE) {
  # Create an empty list to store data frames
  data_list <- list()
  
  all_data <- NULL
  
  for (sample in samples) {
    if (filtered) {
      # Define file types to look for when filtered = TRUE
      file_types <- c("CaRe", "GEX", "KLRB1")
      
      for (file_type in file_types) {
        file_path <- file.path(directory, paste0(sample, "_", file_type, "_filtered.csv"))
        
        # Check if file exists
        if (file.exists(file_path)) {
          # Read the file
          tryCatch({
            sample_data <- read.csv(file_path)
            
            sample_data$sample <- sample
            
            if (is.null(all_data)) {
              all_data <- sample_data
            } else {
              # Get common columns
              common_cols <- intersect(colnames(all_data), colnames(sample_data))
              all_data <- rbind(all_data[, common_cols], sample_data[, common_cols])
            }
            
            cat(sprintf("Successfully read %s\n", paste0(sample, "_", file_type, "_filtered.csv")))
          }, error = function(e) {
            cat(sprintf("Error reading %s: %s\n", file_path, conditionMessage(e)))
          })
        }
      }
    } else {
      # When filtered = FALSE, only look for one file type per sample
      file_path <- file.path(directory, paste0(sample, "_UMI_read_counts_full.csv"))
      
      # Check if file exists
      if (file.exists(file_path)) {
        # Read the file
        tryCatch({
          sample_data <- read.csv(file_path)
          
          sample_data$sample <- sample
          
          if (is.null(all_data)) {
            all_data <- sample_data
          } else {
            # Get common columns
            common_cols <- intersect(colnames(all_data), colnames(sample_data))
            all_data <- rbind(all_data[, common_cols], sample_data[, common_cols])
          }
          
          cat(sprintf("Successfully read %s\n", paste0(sample, "_UMI_read_counts_full.csv")))
        }, error = function(e) {
          cat(sprintf("Error reading %s: %s\n", file_path, conditionMessage(e)))
        })
      }
    }
  }
  
  if (!is.null(all_data)) {
    cat("Total cells loaded:", nrow(all_data), "\n")
    return(all_data)
  } else {
    warning("No data was loaded for any of the provided samples")
    return(NULL)
  }
}
## correlation analysis with plots of relationship between CaRe and SingleCell data
correlation_plots <- function(care_data, single_cell_data, join = "inner") {
  if ("VISERcount" %in% colnames(care_data)) {
    care_data <- rename(care_data, CaRecount = VISERcount)
  }
  required_columns <- list(
    care_data = c("CaRecount", "cellID"),
    single_cell_data = c("SingleCellcount", "cellID")
  )
  
  lapply(names(required_columns), function(dataset) {
    missing_cols <- required_columns[[dataset]][!required_columns[[dataset]] %in% colnames(get(dataset))]
    if (length(missing_cols) > 0) {
      stop(sprintf("Missing required columns in %s: %s", dataset, paste(missing_cols, collapse = ", ")))
    }
  })
  
  join_func <- switch(join,
                      "inner" = inner_join,
                      "left" = left_join,
                      "right" = right_join,
                      "full" = full_join,
                      inner_join
  )

  jointTable <- join_func(care_data, single_cell_data, by = "cellID") %>%
    select(-contains("UMAP")) %>%
    mutate(
      log1pCaReCount = log1p(CaRecount),
      log2CaReCount = log2(CaRecount),
      log1pSingleCell = log1p(SingleCellcount),
    )
  
  cat(
    "Correlation between 10X data and log1p CaRe data:\n",
    "  Pearson: ", round(cor(jointTable$log1pSingleCell, jointTable$log1pCaReCount, method="pearson"), 3),
    "\n  Spearman: ", round(cor(jointTable$log1pSingleCell, jointTable$log1pCaReCount, method="spearman"), 3),
    
    "\n\nCorrelation between 10X data and log2 CaRe data:\n",
    "  Pearson: ", round(cor(jointTable$log1pSingleCell, jointTable$log2CaReCount, method="pearson"), 3),
    "\n  Spearman: ", round(cor(jointTable$log1pSingleCell, jointTable$log2CaReCount, method="spearman"), 3), 
    
    "\n\nCorrelation between 10X counts and log2 CaRe data:\n",
    "  Pearson: ", round(cor(jointTable$SingleCellcount, jointTable$log2CaReCount, method="pearson"), 3),
    "\n  Spearman: ", round(cor(jointTable$SingleCellcount, jointTable$log2CaReCount, method="spearman"), 3),
    
    "\n\nCorrelation between 10X counts and CaRe counts:\n",
    "  Pearson: ", round(cor(jointTable$SingleCellcount, jointTable$CaRecount, method="pearson"), 3),
    "\n  Spearman: ", round(cor(jointTable$SingleCellcount, jointTable$CaRecount, method="spearman"), 3)
  )
  
  transformations <- list(
    CaRecount = "Raw CaRe counts",
    log1pCaReCount = "Log1p CaRe",
    log2CaReCount = "Log2 CaRe"
  )
  
  plots <- lapply(names(transformations), function(trans) {
    ggplot(data = jointTable, aes_string(x = trans, y = "SingleCellcount")) +
      geom_point() +
      labs(
        x = transformations[[trans]],
        y = "SingleCell counts",
        title = paste("SingleCell vs", transformations[[trans]])
      ) +
      theme_minimal()
  })
  names(plots) <- names(transformations)
  
  return(list(
    data = jointTable,
    plots = plots
  ))
}
## compare host gene expression to viral expression
host_virus_correlation <- function(joint_care_data, joint_rds, genes_list, viral_col,
                                   output.dir = "/projects/b1042/GoyalLab/egrody/extractedData/EGS018/joint/hostVirus/", 
                                   outputs = TRUE) {
  if (outputs && !dir.exists(output.dir)) {
    dir.create(output.dir, recursive = TRUE)
  }
  
  for (sample_name in names(joint_care_data)) {
    current_df <- joint_care_data[[sample_name]]
    seurat_obj <- joint_rds[[sample_name]]
    expression_data <- GetAssayData(object = seurat_obj, assay = "RNA", layer = "data")
    expression_df <- data.frame(cellID = colnames(seurat_obj), stringsAsFactors = FALSE)
    for (gene in genes_list) {
      if (gene %in% rownames(seurat_obj)) {
        expression_df[[gene]] <- as.numeric(expression_data[gene, ])
      } else {
        expression_df[[gene]] <- rep(0, ncol(seurat_obj))
        warning(paste("Gene", gene, "not found in dataset. \n"))
      }
    }
    
    genes <- setdiff(colnames(expression_df), "cellID")
    if (!viral_col %in% colnames(current_df)) {
      cat("=== SAMPLE:", sample_name, "===\n")
      cat("Column", viral_col, "not found in data. Skipping.\n\n")
      next
    }
    care <- current_df[[viral_col]]
    if (sd(care) == 0) {
      cat("=== SAMPLE:", sample_name, "===\n")
      cat("CaRe has zero variance (all values identical). Skipping correlation analysis.\n\n")
      next
    }
    results <- data.frame()
    for (gene in genes) {
      gene_expr <- expression_df[[gene]]
      # Use cells where at least one is expressed
      cells_use <- (gene_expr > 0) | (care > 0)
      n_cells <- sum(cells_use)
      if (n_cells < 10) {
        next
      }
      # Check for sufficient variance in both variables
      gene_subset <- gene_expr[cells_use]
      care_subset <- care[cells_use]
      
      if (sd(gene_subset) == 0 || sd(care_subset) == 0) {
        next  # Skip genes with zero variance
      }
      
      # Spearman correlation
      cor_test <- cor.test(gene_expr[cells_use], 
                           care[cells_use],
                           method = "spearman",
                           exact = FALSE)
      
      results <- rbind(results, data.frame(
        gene = gene,
        correlation = as.numeric(cor_test$estimate),
        p_value = cor_test$p.value,
        n_cells = n_cells,
        gene_positive = sum(gene_expr > 0),
        care_positive = sum(care > 0),
        both_positive = sum((gene_expr > 0) & (care > 0))
      ))
    }
    # Adjust p-values for FDR
    results$p_adj <- p.adjust(results$p_value, method = "BH")
    # Sort by significance
    results <- results %>% 
      filter(p_adj < 0.05) %>%
      arrange(desc(abs(correlation)))
    
    # Print results
    cat("===", sample_name, "RESULTS (ordered by magnitude) ===\n\n")
    cat(sprintf("%-12s %10s %12s %12s %8s\n", 
                "Gene", "Corr", "P-value", "P-adj", "Sig"))
    cat(strrep("-", 65), "\n")
    
    for (i in 1:nrow(results)) {
      r <- results[i,]
      sig <- if (is.na(r$p_adj)) "NA" else if (r$p_adj < 0.001) "***" else if (r$p_adj < 0.01) "**" else if (r$p_adj < 0.05) "*" else "ns"
      cat(sprintf("%-12s %10.4f %12.2e %12.2e %8s\n",
                  r$gene, r$correlation, r$p_value, r$p_adj, sig))
    }
    cat(strrep("-", 65), "\n")
    
    # Save outputs if requested
    if (outputs) {
      # Add sample name column for context
      results$sample <- sample_name
      
      # Save CSV
      csv_file <- file.path(output.dir, paste0(sample_name, "_", viral_col, "_correlations.csv"))
      write.csv(results, csv_file, row.names = FALSE)
      
      # Create plots
      pdf_file <- file.path(output.dir, paste0(sample_name, "_", viral_col, "_correlation_barplots.pdf"))
      pdf(pdf_file, width = 10, height = 8)
      
      # Plot 1: Correlation coefficient barplot
      p1 <- ggplot(results, aes(x = reorder(gene, abs(correlation)), y = correlation, 
                                fill = correlation < 0)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
                          labels = c("Positive", "Negative"),
                          name = "Direction") +
        labs(title = paste0(sample_name, " Gene-", viral_col, " Correlations"),
             subtitle = "Spearman correlation coefficients (FDR < 0.05)",
             x = "Gene",
             y = "Correlation Coefficient") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold", size = 14),
              axis.text = element_text(size = 10))
      print(p1)
      dev.off()
      
      pdf_file <- file.path(output.dir, paste0(sample_name, "_", viral_col, "_correlation_scatterplots.pdf"))
      pdf(pdf_file, width = 10, height = 8)
      # Plot 2: Scatter plots for top 6 genes
      top_genes <- head(results$gene, 6)
      
      for (gene in top_genes) {
        gene_expr <- expression_df[[gene]]
        plot_data <- data.frame(
          gene_expr = gene_expr,
          care = care
        )
        
        # Add jitter for better visualization of discrete values
        p2 <- ggplot(plot_data, aes(x = gene_expr, y = care)) +
          geom_point(alpha = 0.3, size = 2) +
          geom_smooth(method = "lm", color = "#d73027", se = TRUE) +
          labs(title = paste(sample_name, "-", gene, "vs", viral_col),
               subtitle = sprintf("Spearman rho = %.3f, p-adj = %.2e",
                                  results$correlation[results$gene == gene],
                                  results$p_adj[results$gene == gene]),
               x = paste(gene, "expression"),
               y = viral_col) +
          theme_minimal() +
          theme(plot.title = element_text(face = "bold"))
        print(p2)
      }
      dev.off()
    }
  }
}

# =============================================================================
# Plotting
# =============================================================================
## volcano gene labels
# Every label keeps a leader line to its own point, so a label that repels far
# from its dot is still unambiguous. min.segment.length = 0 draws the segment at
# any distance (the ggrepel default of 0.5 lines leaves nearby labels floating
# unattached), and the force/padding pair pushes labels far enough off their
# points that the line has something to span.
# max_overlaps stays per-panel: a panel labelling every significant gene needs a
# cap, one labelling a top-N list does not.
volcano_gene_labels <- function(data, mapping, size = 3, max_overlaps = 25, ...) {
  geom_text_repel(
    data = data, mapping = mapping, size = size, fontface = "italic",
    box.padding = 0.5, point.padding = 0.25,
    force = 2, force_pull = 1,
    min.segment.length = 0, segment.size = 0.25, segment.alpha = 0.9,
    max.overlaps = max_overlaps, max.iter = 20000, seed = 42,
    show.legend = FALSE, ...)
}

## tier tables and per-tier DEG, read back from disk
# The tier assignment and the per-contrast DEG tables are written by the
# extraction chunks; these read them back so a plotting chunk needs nothing
# from another chunk's session state.
read_tier_tables <- function(table_dir, samples) {
  out <- lapply(setNames(samples, samples), function(sample_name) {
    f <- paste0(table_dir, sample_name, "_tiers.csv")
    if (!file.exists(f)) return(NULL)
    read.csv(f) %>% mutate(cellID = as.character(cellID),
                           cluster = as.character(cluster))
  })
  out[!vapply(out, is.null, logical(1))]
}

tier_deg_file <- function(deg_dir, sample_name, scope, contrast_name) {
  paste0(deg_dir, sample_name, "_", scope, "_", contrast_name, ".csv")
}

tier_read_deg <- function(deg_dir, sample_name, scope, contrast_name) {
  f <- tier_deg_file(deg_dir, sample_name, scope, contrast_name)
  if (file.exists(f)) read.csv(f) else NULL
}

# NA rather than an error for a test whose inputs are all ties or all missing
tier_safe_p <- function(expr) {
  tryCatch(suppressWarnings(expr), error = function(e) NA_real_)
}

## cNMF meta-program labels
# One name per program across the four per-sample cNMF runs, so a program
# that recurs in several samples carries the same name. Lives here rather
# than in a chunk so any plotting chunk can use it without running another.
# Shared annotation basis, defined on W2 (K=7, the largest run). Invitro, D195
# and Pacute programs are assigned to it by exact maximum-weight INJECTIVE
# matching on top-100 gene Jaccard, so no two programs of a sample share a name.
# Each name is given by its top gene loadings:
#   M1 resting/memory        W2_P1; D195_P1, Invitro_P1, Pacute_P1
#                            BACH2, BCL11B, BCL2, CAMK4, ARID1B
#   M2 activation/adhesion   W2_P2; D195_P2, Invitro_P3, Pacute_P2
#                            IQGAP2, SKAP1, ITGB1, ELMO1, ARHGAP15
#   M3 translation/cytoskel. W2_P3; D195_P3, Invitro_P2, Pacute_P3
#                            ACTB, TMSB4X, B2M, CD52, RPL/RPS
#   M4 cytotoxic EOMES/TOX   W2_P4; D195_P5, Pacute_P6
#                            EOMES, TOX, KALRN, LYST, CST7, FCRL3, GPR183
#   M5 platelet/stromal      W2_P5 only          PF4, PPBP, CLU, APOD
#   M6 B cell/MHC-II         W2_P6; D195_P6      CD74, CIITA
#   M7 myeloid               W2_P7 only          C1QB, CD68, CD163, AIF1
#   M8 granzyme effector     Invitro_P4; Pacute_P5  GZMB, GZMA, NKG7, KLRD1, CCL5
#                            Kept separate from M4: the two share only CST7.
#   M9 cell cycle            Invitro_P5 only
#   M10 EM adhesion          D195_P4 only        S100A4, LMNA, CRIP1, AHR, NIBAN1
#   M11 type I IFN           Pacute_P4 only      MX1, ISG15, IFI6, RSAD2, HERC5
# Pacute_P5 is near-tied between M10 and M8; it is named M8 for its GZMB/NKG7
# content and to match the joint Pacute + W2 run below.

cnmf_metaprogram_labels <- list(
  Invitro = c(Pattern_1 = "M1 resting/memory",
              Pattern_2 = "M3 translation/cytoskeleton",
              Pattern_3 = "M2 activation/adhesion",
              Pattern_4 = "M8 granzyme effector (in vitro)",
              Pattern_5 = "M9 cell cycle (in vitro)"),
  W2      = c(Pattern_1 = "M1 resting/memory",
              Pattern_2 = "M2 activation/adhesion",
              Pattern_3 = "M3 translation/cytoskeleton",
              Pattern_4 = "M4 cytotoxic EOMES/TOX",
              Pattern_5 = "M5 platelet/stromal (contaminant)",
              Pattern_6 = "M6 B cell/MHC-II (contaminant)",
              Pattern_7 = "M7 myeloid (contaminant)"),
  D195    = c(Pattern_1 = "M1 resting/memory",
              Pattern_2 = "M2 activation/adhesion",
              Pattern_3 = "M3 translation/cytoskeleton",
              Pattern_4 = "M10 EM adhesion (ART)",
              Pattern_5 = "M4 cytotoxic EOMES/TOX",
              Pattern_6 = "M6 B cell/MHC-II (contaminant)"),
  Pacute  = c(Pattern_1 = "M1 resting/memory",
              Pattern_2 = "M2 activation/adhesion",
              Pattern_3 = "M3 translation/cytoskeleton",
              Pattern_4 = "M11 type I IFN (acute)",
              Pattern_5 = "M8 granzyme effector",
              Pattern_6 = "M4 cytotoxic EOMES/TOX")
)
stopifnot(!any(vapply(cnmf_metaprogram_labels, anyDuplicated, integer(1)) > 0))

# Programs of the joint Pacute + W2 cNMF run (cleaned, K = 6), named by the
# metaprogram each matches in the per-sample runs: Pattern_1 carries BACH2,
# GPHN, RASGRF2 and SATB1, Pattern_2 INPP4B, FKBP5, SKAP1 and SIK3, Pattern_3
# the actin and ribosomal genes, Pattern_4 GZMB, NKG7 and PRF1, Pattern_5 the
# KALRN/FRMD3/DTHD1 set, and Pattern_6 the MHC-II and myeloid contaminant.
# The tier analysis reads both in vivo acute samples from this run, which puts
# them on identical programs.
cnmf_integrated_metaprogram_labels <- c(
  Pattern_1 = "M1 resting/memory",
  Pattern_2 = "M2 activation/adhesion",
  Pattern_3 = "M3 translation/cytoskeleton",
  Pattern_4 = "M8 granzyme effector",
  Pattern_5 = "M4 cytotoxic EOMES/TOX",
  Pattern_6 = "M6 MHC-II/myeloid (contaminant)")
stopifnot(!anyDuplicated(cnmf_integrated_metaprogram_labels))

## correlation plot
plot_jitter_marginal <- function(df, title = NULL, xlims = c(0, 9.1), ylims = c(0, 9.1)) {
  r2_all <- cor(df$CaRe, df$GEX)^2

  nz <- df %>% dplyr::filter(CaRe > 0, GEX > 0)
  r2_nz <- if (nrow(nz) >= 2) cor(nz$CaRe, nz$GEX)^2 else NA_real_

  if (is.na(r2_nz)) {
    r2_label <- sprintf("R² = %.2f", r2_all)
  } else {
    r2_label <- sprintf("R² = %.2f (all)   R² = %.2f (non-zero)",
                        r2_all, r2_nz)
  }

  p_main <- ggplot(data = df, aes(x = CaRe, y = GEX)) +
    geom_jitter(alpha = 0.7, size = 1.5, width = 0.1, height = 0.1,
                color = "#333333", shape = 16) +
    labs(x = "Ca*Re", y = "TenX") +
    theme_classic(base_size = 14) +
    scale_x_continuous(limits = xlims, expand = c(0, 0),
                       breaks = scales::breaks_width(1)) +
    scale_y_continuous(limits = ylims, expand = c(0, 0),
                       breaks = scales::breaks_width(1))

  top_bins <- hist(df$CaRe, breaks = seq(xlims[1], xlims[2], length.out = 31), plot = FALSE)
  top_df <- data.frame(x = head(top_bins$breaks, -1) + diff(top_bins$breaks) / 2,
                       count = top_bins$counts + 0.01,
                       width = diff(top_bins$breaks))

  p_top <- ggplot(top_df, aes(x = x, y = count)) +
    geom_col(width = top_df$width[1], fill = "gray55", color = "white") +
    scale_x_continuous(limits = xlims, expand = c(0, 0)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                       breaks = c(0, 1, 10, 100, 1000, 10000)) +
    labs(y = "Count", title = r2_label) +
    theme_classic(base_size = 14) +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title   = element_text(hjust = 0.5, size = 12))

  right_bins <- hist(df$GEX, breaks = seq(ylims[1], ylims[2], length.out = 31), plot = FALSE)
  right_df <- data.frame(x = head(right_bins$breaks, -1) + diff(right_bins$breaks) / 2,
                         count = right_bins$counts + 0.01,
                         width = diff(right_bins$breaks))

  p_right <- ggplot(right_df, aes(x = x, y = count)) +
    geom_col(width = right_df$width[1], fill = "gray55", color = "white") +
    scale_x_continuous(limits = ylims, expand = c(0, 0)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                       breaks = c(0, 1, 10, 100, 1000, 10000)) +
    labs(y = "Count") +
    coord_flip() +
    theme_classic(base_size = 14) +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank())

  combined <- p_top + plot_spacer() + p_main + p_right +
    plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

  if (!is.null(title) && nzchar(title)) {
    combined <- combined +
      patchwork::plot_annotation(title = title,
                                 theme = theme(plot.title = element_text(hjust = 0.5,
                                                                         size = 11)))
  }
  # The two R² the panel title carries, plus the cell counts behind them,
  # attached to the returned plot so a caller prints the values on the panel
  # instead of recomputing the correlation.
  attr(combined, "r2") <- list(all = r2_all, non_zero = r2_nz,
                               n = nrow(df), n_non_zero = nrow(nz))
  combined
}

# Print the R² of a plot_jitter_marginal panel: squared Pearson correlation of
# the log1p Ca*Re and TenX counts over all cells, and over cells positive in
# both assays (Ca*Re > 0 & TenX > 0), which is the "non-zero" value in the
# panel title. R² is NA where a group has fewer than two cells or no variance
# in either assay, as in the uninfected samples whose Ca*Re counts are all zero.
print_jitter_r2 <- function(p, label) {
  r2 <- attr(p, "r2")
  fmt <- function(x) if (is.null(x) || is.na(x)) "n/a" else sprintf("%.3f", x)
  cat(sprintf("  %s jitter R² (squared Pearson r, log1p Ca*Re vs log1p TenX):\n", label))
  cat(sprintf("    all cells:      R² = %s (n = %d)\n", fmt(r2$all), r2$n))
  cat(sprintf("    non-zero cells: R² = %s (n = %d)\n", fmt(r2$non_zero), r2$n_non_zero))
}

## plot a UMAP
# point_order controls draw order:
#   "foreground_last" draws zero cells as a flat, semi-transparent background
#                     layer, then non-zero cells on top sorted by `colorby`.
#                     Both layers use the same point size, so the red area on
#                     the panel stays proportional to the number of positive
#                     cells even though they render on top; zero_alpha lets a
#                     buried positive tint through the gray above it.
#   "shuffled"        draws every cell in one permuted layer, so non-zero cells
#                     intercalate with the gray rather than covering it
# The permutation is keyed on cellID when that column is present, so two panels
# of the same cells -- the TenX and Ca*Re views of one sample -- share an
# intercalation and differ only in color and point size.
# low_frac sets where the non-zero ramp starts, as a fraction of the way from
# white to `color`, so a count of 1 is a visible tint rather than background gray.
plotUMAP <- function(data, colorby, title, output_dir, saveas, comparison = FALSE,
                     color = "darkblue",
                     point_order = c("foreground_last", "shuffled"), seed = 1,
                     point_size = 1, zero_alpha = 0.35, low_frac = 0.35) {
  point_order <- match.arg(point_order)
  colorby_str <- if (is.character(substitute(colorby))) {
    colorby
  } else {
    deparse(substitute(colorby))
  }

  if (point_order == "shuffled") {
    key <- if ("cellID" %in% colnames(data)) order(data$cellID) else seq_len(nrow(data))
    set.seed(seed)
    data_ordered <- data[key[sample(length(key))], ]
    # Non-zero points render larger
    data_ordered$size_flag <- ifelse(data_ordered[[colorby_str]] != 0, "nonzero", "zero")

    umap <- ggplot(data_ordered, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(color = .data[[colorby_str]], size = size_flag), shape = 16) +
      scale_size_manual(values = c("zero" = 1, "nonzero" = 2)) +
      guides(size = "none")
  } else {
    zero_points    <- data[data[[colorby_str]] == 0, ]
    nonzero_points <- data[data[[colorby_str]] != 0, ]
    nonzero_points <- nonzero_points[order(nonzero_points[[colorby_str]]), ]

    umap <- ggplot(mapping = aes(x = UMAP1, y = UMAP2)) +
      geom_point(data = zero_points, color = "#E6E7E8", alpha = zero_alpha,
                 size = point_size, shape = 16) +
      geom_point(data = nonzero_points, aes(color = .data[[colorby_str]]),
                 size = point_size, shape = 16)
  }

  umap <- umap +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = rel(0.6)),
          legend.text = element_text(size = rel(0.6), angle = 30),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          aspect.ratio = 1) +
    labs(title = title, color = "")

  if (comparison) {
    mid <- median(data[[colorby_str]])
    umap <- umap + scale_color_gradient2(midpoint = mid, low = "blue", mid = "#E6E7E8", high = "red")
  } else if (point_order == "shuffled") {
    umap <- umap + scale_color_gradient(low = "#E6E7E8", high = color)
  } else {
    low_color <- colorRampPalette(c("white", color))(100)[max(1, round(low_frac * 100))]
    umap <- umap + scale_color_gradient(low = low_color, high = color)
  }
  ggsave(umap, file = paste0(output_dir, saveas),
         device = "pdf")
  
  return(umap)
}

## shared scale for a set of density UMAPs
# Panels that get compared have to smooth alike and share a colour scale, so the
# bandwidth, the contour breaks and the area unit are computed once over every
# column that will be drawn, then handed to each plotDensityUMAP call.
# area_frac is the patch the contour levels are quoted per, as a fraction of the
# UMAP bounding box: at 0.01 a level of 100 means a patch covering 1% of the
# panel holds about 100 positive cells. UMAP coordinates have no scale of their
# own, so levels are quoted per fraction of the panel rather than per unit area.
umapDensityScale <- function(data, cols, area_frac = 0.01, n_breaks = 7,
                             grid_n = 200) {
  lims <- c(range(data$UMAP1), range(data$UMAP2))
  h <- c(MASS::bandwidth.nrd(data$UMAP1), MASS::bandwidth.nrd(data$UMAP2))
  area_unit <- area_frac * diff(range(data$UMAP1)) * diff(range(data$UMAP2))
  peaks <- vapply(cols, function(col) {
    pos <- data[data[[col]] > 0, ]
    if (nrow(pos) < 3) return(0)
    max(MASS::kde2d(pos$UMAP1, pos$UMAP2, h = h, n = grid_n,
                    lims = lims)$z) * nrow(pos) * area_unit
  }, numeric(1))
  peak <- max(peaks)
  if (peak <= 0) {
    return(list(h = h, lims = lims, area_unit = area_unit, area_frac = area_frac,
                breaks = c(0, 1), legend_breaks = c(0, 1)))
  }
  # Columns sharing a scale can differ a lot in magnitude, and a surface whose
  # own peak sits below the first contour draws nothing at all. Space the
  # contours so the weakest column still clears its first one; the legend keeps
  # round numbers of its own, so the ticks stay readable however fine this gets.
  weak <- min(peaks[peaks > 0])
  n_contours <- min(40, max(n_breaks, ceiling(peak / (weak / 2))))
  brks <- seq(0, peak, length.out = n_contours + 1)
  legend_brks <- pretty(c(0, peak), n_breaks)
  list(h = h, lims = lims, area_unit = area_unit, area_frac = area_frac,
       breaks = brks, legend_breaks = legend_brks[legend_brks <= peak])
}

## plot the density of the positive cells on a UMAP
# Every cell is drawn gray, then filled contours of the kernel density of the
# cells with colorby > 0, scaled so a level reads as positive cells per
# area_frac of the panel. The band below the first contour is left unfilled, so
# regions with little signal show the cells rather than a wash of colour.
# A density needs points to smooth: with fewer than min_cells positives the
# contour layer is dropped and the panel is the gray cells alone, rather than
# implying a spatial pattern the counts do not support.
# legend_ends_only labels the colourbar at its two ends (0 and the shared peak)
# instead of at the pretty ticks, so the scale span is what a reader compares
# between panels.
plotDensityUMAP <- function(data, colorby, title, output_dir, saveas,
                            density_scale, color = "red", legend_label = NULL,
                            min_cells = 10, legend_ends_only = FALSE) {
  # colorby takes a bare column name, a literal string, or a variable holding
  # one: evaluate it when that yields a single column name, otherwise deparse
  colorby_expr <- substitute(colorby)
  colorby_str <- tryCatch({
    value <- eval(colorby_expr, parent.frame())
    if (is.character(value) && length(value) == 1) value else deparse(colorby_expr)
  }, error = function(e) deparse(colorby_expr))
  if (!colorby_str %in% colnames(data)) {
    stop(sprintf("plotDensityUMAP: '%s' is not a column of data", colorby_str))
  }
  pos  <- data[data[[colorby_str]] > 0, ]
  brks <- density_scale$breaks
  area_unit <- density_scale$area_unit
  if (is.null(legend_label)) {
    legend_label <- sprintf("positive cells per %g%% of UMAP area",
                            density_scale$area_frac * 100)
  }

  umap <- ggplot() +
    geom_point(data = data, aes(x = UMAP1, y = UMAP2),
               color = "#E6E7E8", size = 0.6, shape = 16)
  if (nrow(pos) >= min_cells) {
    # contour_var = "count" returns level as density * n; area_unit converts it
    # to positive cells per area_frac of the panel
    umap <- umap +
      stat_density_2d(data = pos, aes(x = UMAP1, y = UMAP2,
                                      fill = after_stat(level * area_unit)),
                      geom = "polygon", alpha = 0.6, contour_var = "count",
                      breaks = brks[brks > 0] / area_unit, h = density_scale$h)
  }
  legend_brks <- if (legend_ends_only) range(brks) else density_scale$legend_breaks
  # signif keeps the peak readable without inventing precision; as.character
  # drops the trailing zeros format() would pad on.
  legend_labs <- if (legend_ends_only) {
    as.character(signif(legend_brks, 3))
  } else {
    ggplot2::waiver()
  }
  umap <- umap +
    scale_fill_gradient(low = "#E6E7E8", high = color,
                        limits = range(brks), breaks = legend_brks,
                        labels = legend_labs, name = legend_label) +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = rel(0.6)),
          legend.text = element_text(size = rel(0.5)),
          legend.key.height = unit(0.25, "cm"),
          legend.key.width = unit(1.1, "cm"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          aspect.ratio = 1) +
    labs(title = title)

  ggsave(umap, file = paste0(output_dir, saveas), device = "pdf")

  return(umap)
}

# point_order controls draw order of the non-background identities:
#   "foreground_last" draws each in turn on top of the background layer
#   "shuffled"        draws every cell in one randomly permuted layer, so
#                     colored cells are not preferentially on top
create_identity_umap <- function(data, col1 = NULL, col2 = NULL,
                                 colors = NULL,
                                 title = "", saveas,
                                 point_order = c("foreground_last", "shuffled"),
                                 seed = 1) {
  point_order <- match.arg(point_order)
  required_coords <- c("UMAP1", "UMAP2")
  missing_coords <- required_coords[!required_coords %in% colnames(data)]
  if (length(missing_coords) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_coords, collapse = ", ")))
  }

  # Two-column mode: derive identity from col1/col2 and use default colors
  if (!is.null(col1) && !is.null(col2)) {
    missing_cols <- c(col1, col2)[!c(col1, col2) %in% colnames(data)]
    if (length(missing_cols) > 0) {
      stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
    }
    data <- data %>%
      mutate(identity = case_when(
        .data[[col1]] > 0 & .data[[col2]] > 0 ~ "both",
        .data[[col1]] > 0 ~ col1,
        .data[[col2]] > 0 ~ col2,
        .data[[col1]] == 0 & .data[[col2]] == 0 ~ "neither"
      ))
    if (is.null(colors)) {
      colors <- c("neither" = "#E5E5E5")
      colors[col1] <- "red"
      colors[col2] <- "blue"
      colors["both"] <- "purple"
    }
  }

  if (!"identity" %in% colnames(data)) {
    stop("data must contain an 'identity' column or col1/col2 must be supplied")
  }
  if (is.null(colors)) {
    stop("colors must be provided as a named vector mapping identity levels to colors")
  }

  # First level in colors is the background layer (drawn first, smaller)
  identity_levels <- names(colors)
  data$identity <- factor(data$identity, levels = identity_levels)

  background_level <- identity_levels[1]
  foreground_levels <- identity_levels[-1]

  if (point_order == "shuffled") {
    data$point_size <- ifelse(data$identity == background_level, 1, 2)
    set.seed(seed)
    shuffled <- data[sample(nrow(data)), ]

    identity_umap <- ggplot() +
      geom_point(data = shuffled,
                 aes(x = UMAP1, y = UMAP2, color = identity, size = point_size),
                 shape = 16) +
      scale_size_identity()
  } else {
    bg_points <- data %>% filter(identity == background_level)

    identity_umap <- ggplot() +
      geom_point(data = bg_points,
                 aes(x = UMAP1, y = UMAP2, color = identity),
                 size = 1, shape = 16)

    for (lvl in foreground_levels) {
      fg_points <- data %>% filter(identity == lvl)
      identity_umap <- identity_umap +
        geom_point(data = fg_points,
                   aes(x = UMAP1, y = UMAP2, color = identity),
                   size = 2, shape = 16)
    }
  }

  identity_umap <- identity_umap +
    scale_color_manual(values = colors) +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = rel(0.6)),
          legend.text = element_text(size = rel(0.6), angle = 30),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          aspect.ratio = 1) +
    labs(title = title, color = "Identity")

  ggsave(identity_umap, file = saveas,
         device = "pdf")

  return(identity_umap)
}
## plot FM barcode on UMAP
plotBarcodeUMAP <- function(data, title, min_count = 1, max_count = Inf,
                            output_dir = "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/",
                            saveas = paste0("UMAP_", title, ".pdf")) {
  # Count cells per barcode, excluding "none"
  barcode_counts <- table(data$barcode)
  target_barcodes <- names(barcode_counts[
    barcode_counts >= min_count &
      barcode_counts <= max_count &
      names(barcode_counts) != "none"
  ])
  n_targets <- length(target_barcodes)
  
  if (n_targets == 0) {
    cat("No barcodes found with count in [", min_count, ",", max_count, "]\n")
    return(NULL)
  }
  
  # Assign colors: unique hue per target barcode, gray93 for everything else
  target_colors <- setNames(scales::hue_pal()(n_targets), target_barcodes)
  target_colors["none"] <- "#E5E5E5"
  
  plot_data <- data
  plot_data$barcode <- factor(plot_data$barcode, levels = c("none", target_barcodes))
  plot_data <- arrange(plot_data, barcode)
  
  p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = barcode)) +
    geom_point(aes(size = barcode)) +
    scale_color_manual(values = target_colors, na.value = "#E5E5E5") +
    scale_size_manual(
      values = setNames(c(0.5, rep(0.8, n_targets)), c("none", target_barcodes)),
      na.value = 0.5
    ) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    theme(legend.position = "none") +
    ggtitle(title)
  
  ggsave(p, file = paste0(output_dir, saveas),
         device = "pdf")
  
  return(p)
}
# =============================================================================
# DEG Analysis
# =============================================================================
## perform comprehensive subset analysis
analyze_cell_subset <- function(seurat_obj, 
                                subset_cells_1, 
                                subset_cells_2 = NULL, 
                                min_pct = 0.1,
                                group_names = c("group1", "group2")) {
  
  # Input validation
  if (!is(seurat_obj, "Seurat")) {
    stop("First argument must be a Seurat object")
  }
  
  if (!is.vector(subset_cells_1)) {
    stop("subset_cells_1 must be a vector of cell IDs")
  }
  
  # Determine comparison mode
  two_group_mode <- !is.null(subset_cells_2)
  
  if (two_group_mode) {
    if (!is.vector(subset_cells_2)) {
      stop("subset_cells_2 must be a vector of cell IDs")
    }
    cat("Running two-group comparison mode\n")
  } else {
    cat("Running single-group comparison mode (subset vs rest)\n")
  }
  
  # Convert cell names to character vectors
  subset_cells_1 <- as.character(subset_cells_1)
  if (two_group_mode) {
    subset_cells_2 <- as.character(subset_cells_2)
    
    # Check for overlap between groups
    overlap <- intersect(subset_cells_1, subset_cells_2)
    if (length(overlap) > 0) {
      stop(sprintf("Found %d cells present in both subset_cells_1 and subset_cells_2. Groups must be mutually exclusive.",
                   length(overlap)))
    }
  }
  
  # Validate subset_cells_1
  valid_cells_1 <- subset_cells_1 %in% colnames(seurat_obj)
  if (!any(valid_cells_1)) {
    stop("None of the provided cell IDs in subset_cells_1 were found in the Seurat object")
  }
  
  if (!all(valid_cells_1)) {
    warning(sprintf("%d cells from subset_cells_1 were not found in the Seurat object",
                    sum(!valid_cells_1)))
    subset_cells_1 <- subset_cells_1[valid_cells_1]
  }
  
  # Validate subset_cells_2 if provided
  if (two_group_mode) {
    valid_cells_2 <- subset_cells_2 %in% colnames(seurat_obj)
    if (!any(valid_cells_2)) {
      stop("None of the provided cell IDs in subset_cells_2 were found in the Seurat object")
    }
    
    if (!all(valid_cells_2)) {
      warning(sprintf("%d cells from subset_cells_2 were not found in the Seurat object",
                      sum(!valid_cells_2)))
      subset_cells_2 <- subset_cells_2[valid_cells_2]
    }
  }
  
  # Print diagnostic information
  cat(sprintf("Found %d cells in %s\n", length(subset_cells_1), group_names[1]))
  if (two_group_mode) {
    cat(sprintf("Found %d cells in %s\n", length(subset_cells_2), group_names[2]))
  } else {
    n_rest <- ncol(seurat_obj) - length(subset_cells_1)
    cat(sprintf("Comparing against %d cells in rest of dataset\n", n_rest))
  }
  
  # Add metadata column for group membership
  if (two_group_mode) {
    seurat_obj$comparison_group <- ifelse(
      colnames(seurat_obj) %in% subset_cells_1, group_names[1],
      ifelse(colnames(seurat_obj) %in% subset_cells_2, group_names[2], "excluded")
    )
    # Subset to only include cells in one of the two groups
    seurat_obj <- subset(seurat_obj, cells = c(subset_cells_1, subset_cells_2))
    cat(sprintf("Analysis will include %d total cells\n", ncol(seurat_obj)))
  } else {
    seurat_obj$comparison_group <- ifelse(
      colnames(seurat_obj) %in% subset_cells_1, group_names[1], group_names[2]
    )
  }
  
  # Check for cluster information
  if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    stop("No 'seurat_clusters' column found in metadata. Please run FindClusters() first")
  }
  
  # Get unique clusters for each group
  clusters_group1 <- unique(seurat_obj$seurat_clusters[seurat_obj$comparison_group == group_names[1]])
  clusters_group2 <- unique(seurat_obj$seurat_clusters[seurat_obj$comparison_group == group_names[2]])
  
  cat(sprintf("%s cells are found in %d clusters: %s\n",
              group_names[1], length(clusters_group1),
              paste(sort(clusters_group1), collapse = ", ")))
  cat(sprintf("%s cells are found in %d clusters: %s\n",
              group_names[2], length(clusters_group2),
              paste(sort(clusters_group2), collapse = ", ")))
  
  # Get overlapping clusters
  overlap_clusters <- intersect(clusters_group1, clusters_group2)
  cat(sprintf("Overlapping clusters: %d (%s)\n",
              length(overlap_clusters),
              ifelse(length(overlap_clusters) > 0, 
                     paste(sort(overlap_clusters), collapse = ", "),
                     "none")))
  
  # Initialize results list
  results <- list()
  results$metadata <- list(
    mode = ifelse(two_group_mode, "two_group", "single_group"),
    comparison_summary = sprintf("Comparing %s (n=%d) vs %s (n=%d)", 
                                 group_names[1], length(subset_cells_1),
                                 group_names[2], 
                                 ifelse(two_group_mode, length(subset_cells_2), 
                                        ncol(seurat_obj) - length(subset_cells_1))),
    group1_name = group_names[1],
    group2_name = group_names[2],
    n_cells_group1 = length(subset_cells_1),
    n_cells_group2 = ifelse(two_group_mode, length(subset_cells_2), 
                            ncol(seurat_obj) - length(subset_cells_1)),
    group_names = group_names,
    overlap_clusters = overlap_clusters
  )
  
  # 1. Overall differential expression (across all clusters)
  cat("\n=== Overall Differential Expression Analysis ===\n")
  Idents(seurat_obj) <- factor(seurat_obj$comparison_group, levels = group_names)
  
  results$overall_de <- tryCatch({
    de_results <- FindMarkers(seurat_obj,
                              ident.1 = group_names[1],
                              ident.2 = group_names[2],
                              min.pct = min_pct,
                              test.use = "wilcox",
                              verbose = FALSE)
    
    cat(sprintf("Found %d differentially expressed genes (all p-values)\n", nrow(de_results)))
    de_results
  }, error = function(e) {
    warning(sprintf("Error in overall DE analysis: %s", e$message))
    NULL
  })
  
  # Filter for significant overall markers
  results$overall_de_significant <- tryCatch({
    if (!is.null(results$overall_de) && nrow(results$overall_de) > 0) {
      sig_results <- results$overall_de[
        abs(results$overall_de$avg_log2FC) > 0.5 & 
          results$overall_de$p_val_adj < 0.01,
      ]
      
      if (nrow(sig_results) > 0) {
        n_up <- sum(sig_results$avg_log2FC > 0)
        n_down <- nrow(sig_results) - n_up
        cat(sprintf("Significant markers: %d total (%d up in %s, %d down)\n",
                    nrow(sig_results), n_up, group_names[1], n_down))
        sig_results$gene <- rownames(sig_results)
        sig_results
      } else {
        cat("No markers passed significance thresholds\n")
        NULL
      }
    } else {
      NULL
    }
  }, error = function(e) {
    warning("Error filtering overall markers: ", e$message)
    NULL
  })
  
  # 2. Cluster-specific analysis (only for overlapping clusters)
  if (length(overlap_clusters) > 0) {
    cat("\n=== Cluster-Specific Differential Expression Analysis ===\n")
    cluster_results <- list()
    
    for(cluster in overlap_clusters) {
      # Get cells from this cluster for each group
      cluster_cells <- seurat_obj$seurat_clusters == cluster
      group1_in_cluster <- sum(seurat_obj$comparison_group[cluster_cells] == group_names[1])
      group2_in_cluster <- sum(seurat_obj$comparison_group[cluster_cells] == group_names[2])
      
      # Require at least 3 cells per group for statistical power
      if (group1_in_cluster < 3 || group2_in_cluster < 3) {
        warning(sprintf("Skipping cluster %s: insufficient cells (%s: %d, %s: %d)\n",
                        cluster, group_names[1], group1_in_cluster, 
                        group_names[2], group2_in_cluster))
        next
      }
      
      # Create a temporary Seurat object for this cluster
      cluster_obj <- subset(seurat_obj, seurat_clusters == cluster)
      Idents(cluster_obj) <- factor(cluster_obj$comparison_group, levels = group_names)
      
      # Perform DE analysis within this cluster
      tryCatch({
        de_results <- FindMarkers(cluster_obj,
                                  ident.1 = group_names[1],
                                  ident.2 = group_names[2],
                                  min.pct = min_pct,
                                  test.use = "wilcox",
                                  verbose = FALSE)
        
        cluster_results[[paste0("cluster_", cluster)]] <- de_results
        
        cat(sprintf("Cluster %s: %d DEGs (%s: %d cells, %s: %d cells)\n",
                    cluster, nrow(de_results), 
                    group_names[1], group1_in_cluster,
                    group_names[2], group2_in_cluster))
        
      }, error = function(e) {
        warning(sprintf("Error in DE analysis for cluster %s: %s", cluster, e$message))
      })
    }
    
    results$cluster_specific <- cluster_results
    
    # Filter cluster-specific results for significance
    results$cluster_specific_significant <- tryCatch({
      if (length(cluster_results) > 0) {
        significant_by_cluster <- lapply(names(cluster_results), function(cluster_name) {
          de_results <- cluster_results[[cluster_name]]
          
          if (!is.null(de_results) && nrow(de_results) > 0) {
            sig_results <- de_results[
              abs(de_results$avg_log2FC) > 0.5 & 
                de_results$p_val_adj < 0.01,
            ]
            
            if (nrow(sig_results) > 0) {
              sig_results$gene <- rownames(sig_results)
              sig_results
            } else {
              NULL
            }
          } else {
            NULL
          }
        })
        
        names(significant_by_cluster) <- names(cluster_results)
        significant_by_cluster <- significant_by_cluster[!sapply(significant_by_cluster, is.null)]
        
        if (length(significant_by_cluster) > 0) {
          total_sig <- sum(sapply(significant_by_cluster, nrow))
          cat(sprintf("\nCluster-specific significant markers: %d genes across %d clusters\n",
                      total_sig, length(significant_by_cluster)))
          
          for (cluster_name in names(significant_by_cluster)) {
            n_markers <- nrow(significant_by_cluster[[cluster_name]])
            n_up <- sum(significant_by_cluster[[cluster_name]]$avg_log2FC > 0)
            n_down <- n_markers - n_up
            cat(sprintf("  %s: %d markers (%d up in %s, %d down)\n", 
                        cluster_name, n_markers, group_names[1], n_up, n_down))
          }
          
          significant_by_cluster
        } else {
          cat("\nNo significant cluster-specific markers found\n")
          NULL
        }
      } else {
        NULL
      }
    }, error = function(e) {
      warning("Error filtering cluster-specific markers: ", e$message)
      NULL
    })
  } else {
    cat("\n=== Cluster-Specific Analysis Skipped ===\n")
    cat("No overlapping clusters between groups\n")
    results$cluster_specific <- NULL
    results$cluster_specific_significant <- NULL
  }
  
  # 3. Pseudo-bulk expression profiles
  cat("\n=== Pseudo-bulk Expression ===\n")
  results$pseudobulk_expression <- tryCatch({
    pb_expr <- AggregateExpression(seurat_obj,
                                   group.by = "comparison_group",
                                   assays = "RNA",
                                   return.seurat = FALSE)$RNA
    cat(sprintf("Calculated pseudo-bulk expression for %d genes\n", nrow(pb_expr)))
    pb_expr
  }, error = function(e) {
    warning("Error calculating pseudo-bulk expression: ", e$message)
    NULL
  })
  
  # Print final summary
  cat("\n=== Analysis Summary ===\n")
  cat(sprintf("Comparison: %s\n", results$metadata$comparison_summary))
  cat(sprintf("Mode: %s\n", results$metadata$mode))
  
  return(results)
}
## analyze_cell_subset, without using clusters
analyze_cell_subset_unclustered <- function(seurat_obj, 
                                            subset_cells_1, 
                                            subset_cells_2 = NULL,
                                            min_pct = 0.1,
                                            group_names = c("group1", "group2")) {
  # Check package dependencies
  rlang::check_installed("Seurat", reason = "to perform single-cell analysis")
  rlang::check_installed("dplyr", reason = "for data manipulation")
  
  require(Seurat)
  require(dplyr)
  
  # Input validation
  if (!is(seurat_obj, "Seurat")) {
    stop("First argument must be a Seurat object")
  }
  
  if (!is.vector(subset_cells_1)) {
    stop("subset_cells_1 must be a vector of cell IDs")
  }
  
  # Determine comparison mode
  two_group_mode <- !is.null(subset_cells_2)
  
  if (two_group_mode) {
    if (!is.vector(subset_cells_2)) {
      stop("subset_cells_2 must be a vector of cell IDs")
    }
    cat("Running two-group comparison mode (unclustered)\n")
  } else {
    cat("Running single-group comparison mode (subset vs rest, unclustered)\n")
  }
  
  # Convert cell names to character vectors
  subset_cells_1 <- as.character(subset_cells_1)
  if (two_group_mode) {
    subset_cells_2 <- as.character(subset_cells_2)
    
    # Check for overlap between groups
    overlap <- intersect(subset_cells_1, subset_cells_2)
    if (length(overlap) > 0) {
      stop(sprintf("Found %d cells present in both subset_cells_1 and subset_cells_2. Groups must be mutually exclusive.",
                   length(overlap)))
    }
  }
  
  # Validate subset_cells_1
  valid_cells_1 <- subset_cells_1 %in% colnames(seurat_obj)
  if (!any(valid_cells_1)) {
    stop("None of the provided cell IDs in subset_cells_1 were found in the Seurat object")
  }
  
  if (!all(valid_cells_1)) {
    n_invalid <- sum(!valid_cells_1)
    n_valid <- sum(valid_cells_1)
    warning(sprintf("%d cells from subset_cells_1 were not found in the Seurat object. %d valid cells remain.",
                    n_invalid, n_valid))
    
    # Show some examples of invalid cells for debugging
    if (n_invalid <= 5) {
      cat(sprintf("Invalid cell IDs: %s\n", paste(subset_cells_1[!valid_cells_1], collapse = ", ")))
    }
    
    subset_cells_1 <- subset_cells_1[valid_cells_1]
  }
  
  # Early check for insufficient cells
  if (length(subset_cells_1) < 3) {
    stop(sprintf("Insufficient valid cells in subset_cells_1: found %d, need at least 3 for statistical analysis.\nConsider checking your cell ID format or filtering criteria.",
                 length(subset_cells_1)))
  }
  
  # Validate subset_cells_2 if provided
  if (two_group_mode) {
    valid_cells_2 <- subset_cells_2 %in% colnames(seurat_obj)
    if (!any(valid_cells_2)) {
      stop("None of the provided cell IDs in subset_cells_2 were found in the Seurat object")
    }
    
    if (!all(valid_cells_2)) {
      n_invalid <- sum(!valid_cells_2)
      n_valid <- sum(valid_cells_2)
      warning(sprintf("%d cells from subset_cells_2 were not found in the Seurat object. %d valid cells remain.",
                      n_invalid, n_valid))
      
      # Show some examples of invalid cells for debugging
      if (n_invalid <= 5) {
        cat(sprintf("Invalid cell IDs: %s\n", paste(subset_cells_2[!valid_cells_2], collapse = ", ")))
      }
      
      subset_cells_2 <- subset_cells_2[valid_cells_2]
    }
    
    # Early check for insufficient cells
    if (length(subset_cells_2) < 3) {
      stop(sprintf("Insufficient valid cells in subset_cells_2: found %d, need at least 3 for statistical analysis.\nConsider checking your cell ID format or filtering criteria.",
                   length(subset_cells_2)))
    }
  }
  
  # Print diagnostic information
  cat(sprintf("\nDataset information:\n"))
  cat(sprintf("  - Total cells in dataset: %d\n", ncol(seurat_obj)))
  cat(sprintf("  - %s: %d cells (%.1f%% of total)\n", 
              group_names[1], 
              length(subset_cells_1),
              100 * length(subset_cells_1) / ncol(seurat_obj)))
  
  if (two_group_mode) {
    cat(sprintf("  - %s: %d cells (%.1f%% of total)\n",
                group_names[2],
                length(subset_cells_2),
                100 * length(subset_cells_2) / ncol(seurat_obj)))
  } else {
    n_rest <- ncol(seurat_obj) - length(subset_cells_1)
    cat(sprintf("  - %s: %d cells (%.1f%% of total)\n",
                group_names[2],
                n_rest,
                100 * n_rest / ncol(seurat_obj)))
  }
  
  # Add metadata column for group membership
  if (two_group_mode) {
    seurat_obj$comparison_group <- ifelse(
      colnames(seurat_obj) %in% subset_cells_1, group_names[1],
      ifelse(colnames(seurat_obj) %in% subset_cells_2, group_names[2], "excluded")
    )
    # Subset to only include cells in one of the two groups
    seurat_obj <- subset(seurat_obj, cells = c(subset_cells_1, subset_cells_2))
    cat(sprintf("  - Analysis will include %d total cells\n", ncol(seurat_obj)))
  } else {
    seurat_obj$comparison_group <- ifelse(
      colnames(seurat_obj) %in% subset_cells_1, group_names[1], group_names[2]
    )
  }
  
  # Count cells in each group
  n_group1 <- sum(seurat_obj$comparison_group == group_names[1])
  n_group2 <- sum(seurat_obj$comparison_group == group_names[2])
  
  cat(sprintf("\nComparing:\n"))
  cat(sprintf("  - %s: %d cells\n", group_names[1], n_group1))
  cat(sprintf("  - %s: %d cells\n", group_names[2], n_group2))
  
  # Initialize results list
  results <- list()
  results$metadata <- list(
    mode = ifelse(two_group_mode, "two_group", "single_group"),
    n_cells_group1 = n_group1,
    n_cells_group2 = n_group2,
    group_names = group_names,
    analysis_type = "unclustered"
  )
  
  # Check if we have enough cells for comparison
  if (n_group1 < 3 || n_group2 < 3) {
    stop(sprintf("Insufficient cells for comparison (%s: %d, %s: %d). Need at least 3 in each group.",
                 group_names[1], n_group1, group_names[2], n_group2))
  }
  
  # 1. Perform unclustered differential expression analysis
  cat("\n=== Differential Expression Analysis ===\n")
  
  # Set identities for comparison
  Idents(seurat_obj) <- factor(seurat_obj$comparison_group, levels = group_names)
  
  # Perform differential expression
  results$de_markers <- tryCatch({
    de_results <- FindMarkers(seurat_obj,
                              ident.1 = group_names[1],
                              ident.2 = group_names[2],
                              min.pct = min_pct,
                              test.use = "wilcox",
                              verbose = FALSE)
    
    # Add gene column
    de_results$gene <- rownames(de_results)
    
    cat(sprintf("Successfully tested %d genes\n", nrow(de_results)))
    cat(sprintf("  - Significant at p_val_adj < 0.05: %d genes\n", 
                sum(de_results$p_val_adj < 0.05)))
    cat(sprintf("  - Significant at p_val_adj < 0.01: %d genes\n",
                sum(de_results$p_val_adj < 0.01)))
    
    de_results
    
  }, error = function(e) {
    stop(sprintf("Error in differential expression analysis: %s", e$message))
  })
  
  # 2. Filter for significant markers
  cat("\nFiltering for significant markers (|avg_log2FC| > 0.5, p_val_adj < 0.01)...\n")
  
  results$de_markers_significant <- tryCatch({
    if (!is.null(results$de_markers) && nrow(results$de_markers) > 0) {
      sig_markers <- results$de_markers[
        abs(results$de_markers$avg_log2FC) > 0.5 & 
          results$de_markers$p_val_adj < 0.01,
      ]
      
      if (nrow(sig_markers) > 0) {
        # Count up and down
        n_up <- sum(sig_markers$avg_log2FC > 0)
        n_down <- sum(sig_markers$avg_log2FC < 0)
        
        cat(sprintf("Found %d significant markers\n", nrow(sig_markers)))
        cat(sprintf("  - Upregulated in %s: %d genes\n", group_names[1], n_up))
        cat(sprintf("  - Downregulated in %s: %d genes\n", group_names[1], n_down))
        
        sig_markers
      } else {
        cat("No markers passed significance thresholds\n")
        NULL
      }
    } else {
      cat("No differential expression results to filter\n")
      NULL
    }
  }, error = function(e) {
    warning("Error filtering significant markers: ", e$message)
    NULL
  })
  
  # 3. Additional filtering: moderate fold change threshold
  results$de_markers_moderate <- tryCatch({
    if (!is.null(results$de_markers) && nrow(results$de_markers) > 0) {
      moderate_markers <- results$de_markers[
        abs(results$de_markers$avg_log2FC) > 0.25 & 
          results$de_markers$p_val_adj < 0.05,
      ]
      
      if (nrow(moderate_markers) > 0) {
        cat(sprintf("\nAlternative threshold (|avg_log2FC| > 0.25, p_val_adj < 0.05):\n"))
        cat(sprintf("  - Found %d markers\n", nrow(moderate_markers)))
        moderate_markers
      } else {
        NULL
      }
    } else {
      NULL
    }
  }, error = function(e) {
    warning("Error applying moderate thresholds: ", e$message)
    NULL
  })
  
  # 4. Calculate pseudo-bulk expression profiles
  cat("\n=== Pseudo-bulk Expression ===\n")
  results$pseudobulk_expression <- tryCatch({
    pb_expr <- AggregateExpression(seurat_obj,
                                   group.by = "comparison_group",
                                   assays = "RNA",
                                   return.seurat = FALSE)$RNA
    cat(sprintf("Calculated pseudo-bulk expression for %d genes\n", nrow(pb_expr)))
    pb_expr
  }, error = function(e) {
    warning("Error calculating pseudo-bulk expression: ", e$message)
    NULL
  })
  
  # Summary of analysis
  cat("\n=== Analysis Summary ===\n")
  cat(sprintf("Mode: %s\n", results$metadata$mode))
  cat(sprintf("Analysis type: %s\n", results$metadata$analysis_type))
  cat(sprintf("Groups: %s (%d cells) vs %s (%d cells)\n",
              group_names[1], n_group1,
              group_names[2], n_group2))
  cat(sprintf("Total genes tested: %d\n", 
              ifelse(!is.null(results$de_markers), 
                     nrow(results$de_markers), 0)))
  cat(sprintf("Significant markers (strict): %s\n",
              ifelse(!is.null(results$de_markers_significant),
                     sprintf("%d genes", nrow(results$de_markers_significant)),
                     "None")))
  cat(sprintf("Significant markers (moderate): %s\n",
              ifelse(!is.null(results$de_markers_moderate),
                     sprintf("%d genes", nrow(results$de_markers_moderate)),
                     "None")))
  
  if (!is.null(results$de_markers_significant) && 
      nrow(results$de_markers_significant) > 0) {
    sig_markers <- results$de_markers_significant
    n_up <- sum(sig_markers$avg_log2FC > 0)
    n_down <- sum(sig_markers$avg_log2FC < 0)
    
    cat(sprintf("  - Upregulated in %s: %d\n", group_names[1], n_up))
    cat(sprintf("  - Downregulated in %s: %d\n", group_names[1], n_down))
    
    # Report top markers only if they exist
    if (n_up > 0) {
      top_up <- head(sig_markers[order(-sig_markers$avg_log2FC), ], 5)
      cat(sprintf("\nTop markers upregulated in %s:\n", group_names[1]))
      for (i in 1:nrow(top_up)) {
        cat(sprintf("  %d. %s (log2FC = %.2f, p_adj = %.2e)\n",
                    i, top_up$gene[i], top_up$avg_log2FC[i], top_up$p_val_adj[i]))
      }
    }
    
    if (n_down > 0) {
      top_down <- head(sig_markers[order(sig_markers$avg_log2FC), ], 5)
      cat(sprintf("\nTop markers downregulated in %s:\n", group_names[1]))
      for (i in 1:nrow(top_down)) {
        cat(sprintf("  %d. %s (log2FC = %.2f, p_adj = %.2e)\n",
                    i, top_down$gene[i], top_down$avg_log2FC[i], top_down$p_val_adj[i]))
      }
    }
  }
  
  cat("\n")
  
  return(results)
}
## visualize subset analysis results
plot_subset_analysis <- function(seurat_obj, analysis_results, subset_cells, save = TRUE, output_dir = "/projects/b1042/GoyalLab/egrody/") {
  require(ggplot2)
  
  if (!is(seurat_obj, "Seurat")) {
    stop("First argument must be a Seurat object")
  }
  
  if (!is.list(analysis_results)) {
    stop("analysis_results must be the output from analyze_cell_subset()")
  }
  
  plots <- list()
  
  # 1. UMAP highlighting subset cells
  if ("umap" %in% names(seurat_obj@reductions)) {
    seurat_obj$in_subset <- colnames(seurat_obj) %in% subset_cells
    plots$umap <- DimPlot(seurat_obj, 
                          group.by = "in_subset", 
                          label = TRUE) +
      ggtitle("Subset Cells in UMAP")
  } else {
    warning("No UMAP reduction found. Skipping UMAP plot.")
  }
  
  # 2. Expression heatmap of top markers
  if (!is.null(analysis_results$conserved_markers) && 
      nrow(analysis_results$conserved_markers) > 0) {
    top_conserved <- head(rownames(analysis_results$conserved_markers), 20)
    plots$heatmap <- DoHeatmap(seurat_obj, 
                               features = top_conserved,
                               group.by = c("in_subset", "seurat_clusters"))
  }
  
  # 3. Violin plots for top markers by cluster
  if (!is.null(analysis_results$top_genes)) {
    for(cluster in names(analysis_results$top_genes)) {
      top_genes <- analysis_results$top_genes[[cluster]]
      if (length(top_genes) > 0) {
        plots[[paste0("violin_", cluster)]] <- VlnPlot(seurat_obj,
                                                       features = head(top_genes, 5),
                                                       split.by = "in_subset",
                                                       group.by = "seurat_clusters")
      }
    }
  }
  
  if (save == TRUE) {
    # Create output directory if it doesn't exist
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Function to safely save a plot
    save_plot <- function(plot, filename, w = width, h = height) {
      if (!is.null(plot)) {
        ggsave(
          filename = file.path(output_dir, filename),
          plot = plot,
          width = w,
          height = h
        )
        return(TRUE)
      }
      return(FALSE)
    }
    
    # Initialize counters for summary
    saved_plots <- character()
    
    # 1. Save UMAP plot
    if (!is.null(plots$umap)) {
      if (save_plot(plots$umap, "umap_subset_cells.png", w = 8, h = 8)) {
        saved_plots <- c(saved_plots, "UMAP plot")
      }
    }
    
    # 2. Save heatmap
    if (!is.null(plots$heatmap)) {
      if (save_plot(plots$heatmap, "marker_heatmap.png", w = 12, h = 10)) {
        saved_plots <- c(saved_plots, "Marker heatmap")
      }
    }
    
    print("Plots saved")
    }
  
  return(plots)
}
## save analysis results to files
save_subset_analysis <- function(analysis_results, output_dir) {
  require(openxlsx)
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 0. Save overall differential expression results
  if (!is.null(analysis_results$overall_de) && nrow(analysis_results$overall_de) > 0) {
    wb <- createWorkbook()
    
    # All results
    de_results <- analysis_results$overall_de
    de_results$gene <- rownames(de_results)
    addWorksheet(wb, "all_genes")
    writeData(wb, "all_genes", de_results)
    
    # Save workbook
    saveWorkbook(wb, file.path(output_dir, "overall_DE.xlsx"), overwrite = TRUE)
  }
  
  # 0b. Save overall SIGNIFICANT markers
  if (!is.null(analysis_results$overall_de_significant) && 
      nrow(analysis_results$overall_de_significant) > 0) {
    
    sig_results <- analysis_results$overall_de_significant
    
    # Split into upregulated and downregulated
    upregulated <- sig_results[sig_results$avg_log2FC > 0, ]
    downregulated <- sig_results[sig_results$avg_log2FC < 0, ]
    
    # Sort by fold change
    if (nrow(upregulated) > 0) {
      upregulated <- upregulated[order(-upregulated$avg_log2FC), ]
    }
    if (nrow(downregulated) > 0) {
      downregulated <- downregulated[order(downregulated$avg_log2FC), ]
    }
    
    wb <- createWorkbook()
    
    # Add upregulated genes tab
    if (nrow(upregulated) > 0) {
      addWorksheet(wb, "upregulated")
      writeData(wb, "upregulated", upregulated)
    }
    
    # Add downregulated genes tab
    if (nrow(downregulated) > 0) {
      addWorksheet(wb, "downregulated")
      writeData(wb, "downregulated", downregulated)
    }
    
    # Save workbook
    saveWorkbook(wb, file.path(output_dir, "overall_significant.xlsx"), overwrite = TRUE)
    
    # Print summary
    cat(sprintf("\nOverall significant markers saved:\n"))
    cat(sprintf("  - Upregulated genes: %d\n", nrow(upregulated)))
    cat(sprintf("  - Downregulated genes: %d\n", nrow(downregulated)))
  }
  
  # 1. Save cluster-specific differential expression results
  if (length(analysis_results$cluster_specific) > 0) {
    # Create a workbook for cluster-specific results
    wb <- createWorkbook()
    
    for (cluster_name in names(analysis_results$cluster_specific)) {
      de_results <- analysis_results$cluster_specific[[cluster_name]]
      if (!is.null(de_results) && nrow(de_results) > 0) {
        # Add cluster results as a worksheet
        addWorksheet(wb, cluster_name)
        # Add data with gene names as a column
        de_results$gene <- rownames(de_results)
        writeData(wb, cluster_name, de_results)
      }
    }
    
    # Save the workbook
    saveWorkbook(wb, file.path(output_dir, "cluster_specific_DE.xlsx"), overwrite = TRUE)
  }
  
  # 1b. Save cluster-specific SIGNIFICANT markers
  if (!is.null(analysis_results$cluster_specific_significant) &&
      length(analysis_results$cluster_specific_significant) > 0) {
    
    wb <- createWorkbook()
    
    for (cluster_name in names(analysis_results$cluster_specific_significant)) {
      sig_results <- analysis_results$cluster_specific_significant[[cluster_name]]
      
      if (!is.null(sig_results) && nrow(sig_results) > 0) {
        # Split into upregulated and downregulated
        upregulated <- sig_results[sig_results$avg_log2FC > 0, ]
        downregulated <- sig_results[sig_results$avg_log2FC < 0, ]
        
        # Sort by fold change
        if (nrow(upregulated) > 0) {
          upregulated <- upregulated[order(-upregulated$avg_log2FC), ]
        }
        if (nrow(downregulated) > 0) {
          downregulated <- downregulated[order(downregulated$avg_log2FC), ]
        }
        
        # Add worksheet for upregulated genes
        if (nrow(upregulated) > 0) {
          sheet_name_up <- paste0(cluster_name, "_up")
          # Truncate sheet name if too long (Excel limit is 31 characters)
          if (nchar(sheet_name_up) > 31) {
            sheet_name_up <- substr(sheet_name_up, 1, 31)
          }
          addWorksheet(wb, sheet_name_up)
          writeData(wb, sheet_name_up, upregulated)
        }
        
        # Add worksheet for downregulated genes
        if (nrow(downregulated) > 0) {
          sheet_name_down <- paste0(cluster_name, "_down")
          # Truncate sheet name if too long
          if (nchar(sheet_name_down) > 31) {
            sheet_name_down <- substr(sheet_name_down, 1, 31)
          }
          addWorksheet(wb, sheet_name_down)
          writeData(wb, sheet_name_down, downregulated)
        }
      }
    }
    
    # Save workbook
    saveWorkbook(wb, file.path(output_dir, "cluster_specific_significant.xlsx"), overwrite = TRUE)
  }
  
  # 2. Save conserved markers
  if (!is.null(analysis_results$conserved_markers) && 
      nrow(analysis_results$conserved_markers) > 0) {
    # Add gene names as a column
    conserved_markers <- analysis_results$conserved_markers
    conserved_markers$gene <- rownames(conserved_markers)
    
    # Create a new workbook for conserved markers
    wb <- createWorkbook()
    addWorksheet(wb, "conserved_markers")
    writeData(wb, "conserved_markers", conserved_markers)
    saveWorkbook(wb, file.path(output_dir, "conserved_markers.xlsx"), overwrite = TRUE)
  }
  
  # 2b. Save significant conserved markers
  if (!is.null(analysis_results$conserved_markers_significant) && 
      nrow(analysis_results$conserved_markers_significant) > 0) {
    # Add gene names as a column
    conserved_markers_sig <- analysis_results$conserved_markers_significant
    conserved_markers_sig$gene <- rownames(conserved_markers_sig)
    
    # Determine which fold change column to use for sorting
    fc_column <- if ("mean_log2FC" %in% colnames(conserved_markers_sig)) {
      "mean_log2FC"
    } else if ("avg_log2FC" %in% colnames(conserved_markers_sig)) {
      "avg_log2FC"
    } else {
      stop("No log2FC column found in conserved_markers_significant")
    }
    
    # Split into upregulated and downregulated
    upregulated <- conserved_markers_sig[conserved_markers_sig[[fc_column]] > 0, ]
    downregulated <- conserved_markers_sig[conserved_markers_sig[[fc_column]] < 0, ]
    
    # Sort by fold change (descending for up, ascending for down)
    if (nrow(upregulated) > 0) {
      upregulated <- upregulated[order(-upregulated[[fc_column]]), ]
    }
    if (nrow(downregulated) > 0) {
      downregulated <- downregulated[order(downregulated[[fc_column]]), ]
    }
    
    # Create a new workbook for significant conserved markers
    wb <- createWorkbook()
    
    # Add upregulated genes tab
    if (nrow(upregulated) > 0) {
      addWorksheet(wb, "upregulated")
      writeData(wb, "upregulated", upregulated)
    }
    
    # Add downregulated genes tab
    if (nrow(downregulated) > 0) {
      addWorksheet(wb, "downregulated")
      writeData(wb, "downregulated", downregulated)
    }
    
    # Save workbook
    saveWorkbook(wb, file.path(output_dir, "conserved_markers_significant.xlsx"), overwrite = TRUE)
    
    # Print summary
    cat(sprintf("\nSignificant conserved markers saved:\n"))
    cat(sprintf("  - Upregulated genes: %d\n", nrow(upregulated)))
    cat(sprintf("  - Downregulated genes: %d\n", nrow(downregulated)))
  }
  
  # 3. Save average expression profiles
  if (!is.null(analysis_results$avg_expression)) {
    avg_exp <- analysis_results$avg_expression
    avg_exp_df <- as.data.frame(avg_exp)
    avg_exp_df$gene <- rownames(avg_exp_df)
    
    wb <- createWorkbook()
    addWorksheet(wb, "average_expression")
    writeData(wb, "average_expression", avg_exp_df)
    saveWorkbook(wb, file.path(output_dir, "average_expression.xlsx"), overwrite = TRUE)
  }
  
  # 4. Save pseudobulk expression profiles if present
  if (!is.null(analysis_results$pseudobulk_expression)) {
    pseudobulk_exp <- analysis_results$pseudobulk_expression
    pseudobulk_exp_df <- as.data.frame(pseudobulk_exp)
    pseudobulk_exp_df$gene <- rownames(pseudobulk_exp_df)
    
    wb <- createWorkbook()
    addWorksheet(wb, "pseudobulk_expression")
    writeData(wb, "pseudobulk_expression", pseudobulk_exp_df)
    saveWorkbook(wb, file.path(output_dir, "pseudobulk_expression.xlsx"), overwrite = TRUE)
  }
  
  # 5. Save top genes list
  if (length(analysis_results$top_genes) > 0) {
    # Convert list to data frame
    max_genes <- max(sapply(analysis_results$top_genes, length))
    top_genes_df <- data.frame(matrix(NA, nrow = max_genes, 
                                      ncol = length(analysis_results$top_genes)))
    colnames(top_genes_df) <- names(analysis_results$top_genes)
    
    for (cluster in names(analysis_results$top_genes)) {
      top_genes_df[[cluster]] <- c(analysis_results$top_genes[[cluster]],
                                   rep(NA, max_genes - length(analysis_results$top_genes[[cluster]])))
    }
    
    wb <- createWorkbook()
    addWorksheet(wb, "top_genes")
    writeData(wb, "top_genes", top_genes_df)
    saveWorkbook(wb, file.path(output_dir, "top_genes.xlsx"), overwrite = TRUE)
  }
  
  # Create a summary file
  summary_text <- c(
    "Subset Analysis Results Summary",
    "============================",
    "",
    "Comparison Details:",
    paste("  Group 1:", analysis_results$metadata$group1_name, 
          sprintf("(n=%d cells)", analysis_results$metadata$n_cells_group1)),
    paste("  Group 2:", analysis_results$metadata$group2_name,
          sprintf("(n=%d cells)", analysis_results$metadata$n_cells_group2)),
    paste("  Mode:", analysis_results$metadata$mode),
    "",
    "Results:",
    paste("  Overall DE genes:", 
          ifelse(!is.null(analysis_results$overall_de), 
                 nrow(analysis_results$overall_de), 0)),
    paste("  Significant overall markers:", 
          ifelse(!is.null(analysis_results$overall_de_significant), 
                 nrow(analysis_results$overall_de_significant), 0)),
    paste("  Clusters analyzed:", length(analysis_results$cluster_specific)),
    paste("  Clusters with significant markers:",
          ifelse(!is.null(analysis_results$cluster_specific_significant),
                 length(analysis_results$cluster_specific_significant), 0)),
    paste("  Conserved markers:", 
          ifelse(!is.null(analysis_results$conserved_markers), 
                 nrow(analysis_results$conserved_markers), 0)),
    paste("  Significant conserved markers:", 
          ifelse(!is.null(analysis_results$conserved_markers_significant), 
                 nrow(analysis_results$conserved_markers_significant), 0)),
    "",
    "Files generated:",
    "- overall_DE.xlsx: Overall differential expression results (all genes)",
    paste0("    * Positive log2FC = higher in ", analysis_results$metadata$group1_name),
    paste0("    * Negative log2FC = higher in ", analysis_results$metadata$group2_name),
    "- overall_significant.xlsx: Significant overall markers (|avg_log2FC| > 0.5, p_val_adj < 0.01)",
    paste0("    * upregulated tab: Genes higher in ", analysis_results$metadata$group1_name),
    paste0("    * downregulated tab: Genes higher in ", analysis_results$metadata$group2_name),
    "- cluster_specific_DE.xlsx: Differential expression results for each cluster (all genes)",
    "- cluster_specific_significant.xlsx: Significant markers per cluster (|avg_log2FC| > 0.5, p_val_adj < 0.01)",
    "    * One tab per cluster and direction (e.g., cluster_0_up, cluster_0_down)",
    "    * Each tab sorted by fold change magnitude",
    "- conserved_markers.xlsx: Markers conserved across clusters",
    "- conserved_markers_significant.xlsx: Significant conserved markers (|avg_log2FC| > 0.5, p_val_adj < 0.01)",
    "    * upregulated tab: Genes with positive fold change (sorted high to low)",
    "    * downregulated tab: Genes with negative fold change (sorted low to high)",
    "- pseudobulk_expression.xlsx: Pseudobulk expression profiles",
    "- top_genes.xlsx: Top differential genes by cluster"
  )
  
  writeLines(summary_text, file.path(output_dir, "analysis_summary.txt"))
  
  cat("Results saved to directory:", output_dir, "\n")
  cat("Generated files:\n")
  list.files(output_dir, pattern = "\\.xlsx$|\\.txt$") %>% 
    paste0("- ", .) %>%
    cat(sep = "\n")
}

# =============================================================================
# cNMF
# =============================================================================

run_preranked_gsea <- function(spectra_scores, gene_sets, pval_cutoff = 0.05) {
  # fgsea is in the "enrichment" optional group, not attached at startup.
  require(fgsea)
  patterns <- colnames(spectra_scores)
  all_results <- list()
  
  for (pat in patterns) {
    ranks <- spectra_scores[[pat]]
    names(ranks) <- rownames(spectra_scores)
    ranks <- ranks[is.finite(ranks)]
    ranks <- sort(ranks, decreasing = TRUE)
    
    res <- fgsea(
      pathways = gene_sets,
      stats = ranks,
      minSize = 10,
      maxSize = 500,
      nPermSimple = 10000,
      eps = 0
    )
    res$pattern <- pat
    all_results[[pat]] <- res
  }
  
  combined <- rbindlist(all_results)
  combined <- combined[combined$padj < pval_cutoff, ]
  combined <- combined[order(combined$pattern, combined$padj), ]
  return(as.data.frame(combined))
}
get_top_n_genes <- function(spectra_scores, n = 50) {
  top_list <- list()
  for (pat in colnames(spectra_scores)) {
    scores <- spectra_scores[[pat]]
    names(scores) <- rownames(spectra_scores)
    top_list[[pat]] <- names(sort(scores, decreasing = TRUE))[1:n]
  }
  return(top_list)
}
run_ora_hypergeometric <- function(query_genes, gene_sets, background_genes,
                                   min_size = 10, max_size = 500,
                                   pval_cutoff = 0.05) {
  # One-sided Fisher/hypergeometric test for overrepresentation
  bg <- unique(background_genes)
  query <- intersect(query_genes, bg)
  N <- length(bg)
  n_query <- length(query)
  
  results <- list()
  for (set_name in names(gene_sets)) {
    set_genes <- intersect(gene_sets[[set_name]], bg)
    K <- length(set_genes)
    if (K < min_size || K > max_size) next
    
    overlap <- intersect(query, set_genes)
    k <- length(overlap)
    if (k == 0) next
    
    pval <- phyper(k - 1, K, N - K, n_query, lower.tail = FALSE)
    
    results[[length(results) + 1]] <- data.frame(
      pathway = set_name,
      pval = pval,
      overlap = k,
      set_size = K,
      query_size = n_query,
      background_size = N,
      overlap_genes = paste(overlap, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(results) == 0) {
    return(data.frame(pathway = character(), pval = numeric(), padj = numeric(),
                      overlap = integer(), set_size = integer(),
                      query_size = integer(), background_size = integer(),
                      overlap_genes = character(), stringsAsFactors = FALSE))
  }
  
  res_df <- do.call(rbind, results)
  res_df$padj <- p.adjust(res_df$pval, method = "BH")
  res_df <- res_df[res_df$padj < pval_cutoff, ]
  res_df <- res_df[order(res_df$padj), ]
  return(res_df)
}
add_cnmf_to_seurat <- function(seurat_obj, usage) {
  # Align cell barcodes
  common_cells <- intersect(Cells(seurat_obj), rownames(usage))
  cat("Matching cells:", length(common_cells), "of", ncol(seurat_obj), "\n")
  
  # Add pattern weights to metadata
  colnames(usage) <- paste0("Pattern_", seq_len(ncol(usage)))
  for (pat in colnames(usage)) {
    seurat_obj[[pat]] <- NA
    seurat_obj[[pat]][common_cells, ] <- usage[common_cells, pat]
  }
  
  return(seurat_obj)
}
## comparisons
find_pattern_markers <- function(spectra_tpm,
                                 spectra_score = NULL,
                                 n_markers = 50,
                                 method = "weighted_zscore",
                                 min_score = 0,
                                 exclude_shared = FALSE,
                                 shared_threshold = 1.5) {
  
  stopifnot(is.matrix(spectra_tpm) || is.data.frame(spectra_tpm))
  spectra_tpm <- as.matrix(spectra_tpm)
  stopifnot(all(spectra_tpm >= 0))
  
  patterns <- colnames(spectra_tpm)
  genes <- rownames(spectra_tpm)
  n_patterns <- ncol(spectra_tpm)
  
  if (is.null(patterns)) patterns <- paste0("pattern_", seq_len(n_patterns))
  if (is.null(genes)) stop("spectra_tpm must have rownames (gene names).")
  
  # Compute Z-scores per gene if not provided
  if (is.null(spectra_score)) {
    row_means <- rowMeans(spectra_tpm)
    row_sds <- apply(spectra_tpm, 1, sd)
    row_sds[row_sds == 0] <- 1  # avoid division by zero for constant rows
    spectra_score <- (spectra_tpm - row_means) / row_sds
  } else {
    spectra_score <- as.matrix(spectra_score)
    stopifnot(identical(dim(spectra_score), dim(spectra_tpm)))
  }
  
  # Build per-gene, per-pattern stats (vectorized)
  # For each gene and each pattern, compute the ratio to the best competing pattern
  # and other useful quantities.
  stats_list <- vector("list", n_patterns)
  
  for (i in seq_len(n_patterns)) {
    tpm_this <- spectra_tpm[, i]
    zscore_this <- spectra_score[, i]
    
    # Max loading across all OTHER patterns for each gene
    if (n_patterns == 2) {
      max_other <- spectra_tpm[, -i]
    } else {
      max_other <- apply(spectra_tpm[, -i, drop = FALSE], 1, max)
    }
    
    ratio <- tpm_this / (max_other + 1e-10)
    
    stats_list[[i]] <- data.frame(
      gene = genes,
      pattern = patterns[i],
      tpm = tpm_this,
      zscore = zscore_this,
      max_other_tpm = max_other,
      ratio = ratio,
      stringsAsFactors = FALSE
    )
  }
  
  marker_df <- do.call(rbind, stats_list)
  rownames(marker_df) <- NULL
  
  # Apply minimum score filter
  if (min_score > 0) {
    marker_df <- marker_df[marker_df$tpm >= min_score, ]
  }
  
  # Apply shared-gene exclusion filter
  if (exclude_shared) {
    marker_df <- marker_df[marker_df$ratio >= shared_threshold, ]
  }
  
  # Compute the composite marker_score based on method
  marker_df$marker_score <- switch(method,
                                   
                                   "weighted_zscore" = {
                                     # Z-score captures specificity (how unusually high this gene loads on
                                     # this pattern relative to its own distribution). Multiplying by
                                     # log1p(TPM) upweights genes that are also highly expressed.
                                     # log1p prevents a handful of very-high-TPM housekeeping genes from
                                     # dominating; the log compression keeps magnitude informative but
                                     # bounded.
                                     marker_df$zscore * log1p(marker_df$tpm)
                                   },
                                   
                                   "softmax_weighted" = {
                                     # For each gene, compute softmax across patterns on TPM values.
                                     # softmax(x_i) = exp(x_i) / sum(exp(x_j)) gives a specificity
                                     # proportion that amplifies differences more than a linear proportion
                                     # would. Multiply by TPM to weight by magnitude.
                                     # Vectorized: subtract row max for numerical stability, then compute.
                                     row_max <- apply(spectra_tpm, 1, max)
                                     scaled <- spectra_tpm - row_max  # subtract max per row
                                     exp_scaled <- exp(scaled)
                                     softmax_mat <- exp_scaled / rowSums(exp_scaled)
                                     score_mat <- softmax_mat * spectra_tpm
                                     # Map back to the long-format marker_df
                                     as.vector(score_mat[cbind(
                                       match(marker_df$gene, genes),
                                       match(marker_df$pattern, patterns)
                                     )])
                                   },
                                   
                                   "rank_product" = {
                                     # Non-parametric approach: for each pattern, rank genes by
                                     # (1) specificity (ratio to max_other) and (2) magnitude (tpm).
                                     # Combine as reciprocal of the product of ranks. This is robust
                                     # to skewed distributions and outliers. Genes that rank well on
                                     # BOTH dimensions get the highest score.
                                     rp_scores <- numeric(nrow(marker_df))
                                     for (pat in patterns) {
                                       idx <- which(marker_df$pattern == pat)
                                       pat_data <- marker_df[idx, ]
                                       n_genes_pat <- nrow(pat_data)
                                       rank_spec <- rank(-pat_data$ratio, ties.method = "average")
                                       rank_mag <- rank(-pat_data$tpm, ties.method = "average")
                                       rp_scores[idx] <- 1 / (rank_spec * rank_mag)
                                     }
                                     rp_scores
                                   },
                                   
                                   "ratio" = {
                                     # Pure specificity baseline: fold-enrichment over best competitor.
                                     # No magnitude weighting. Useful for comparison but tends to select
                                     # low-expression genes that happen to be very pattern-specific.
                                     marker_df$ratio
                                   },
                                   
                                   stop("Unknown method: ", method,
                                        ". Choose from: weighted_zscore, softmax_weighted, rank_product, ratio")
  )
  
  # Select top markers per pattern
  markers <- list()
  for (pat in patterns) {
    pat_df <- marker_df[marker_df$pattern == pat, ]
    pat_df <- pat_df[order(pat_df$marker_score, decreasing = TRUE), ]
    markers[[pat]] <- head(pat_df$gene, n_markers)
  }
  
  return(list(
    markers = markers,
    marker_stats = marker_df,
    method = method
  ))
}
compare_pattern_usage <- function(usage, group1_cellids, group2_cellids = NULL,
                                  group_labels = c("infected", "uninfected")) {
  # If only one list given, derive the second as all remaining cells
  if (is.null(group2_cellids)) {
    group2_cellids <- setdiff(rownames(usage), group1_cellids)
  }
  
  # Subset pattern matrices
  group1_patterns <- usage[rownames(usage) %in% group1_cellids, , drop = FALSE]
  group2_patterns <- usage[rownames(usage) %in% group2_cellids, , drop = FALSE]
  
  # Calculate summary statistics
  group1_means <- colMeans(group1_patterns)
  group2_means <- colMeans(group2_patterns)
  
  comparison_df <- data.frame(
    pattern = colnames(usage),
    group1_mean = group1_means,
    group2_mean = group2_means,
    log2_fold_change = log2((group1_means + 1e-10) / (group2_means + 1e-10))
  )
  
  # Statistical testing with Wilcoxon rank-sum test
  stat_results <- lapply(colnames(usage), function(pat) {
    test <- wilcox.test(group1_patterns[, pat], group2_patterns[, pat])
    data.frame(
      pattern = pat,
      p_value = test$p.value,
      group1_median = median(group1_patterns[, pat]),
      group2_median = median(group2_patterns[, pat])
    )
  })
  
  stat_df <- do.call(rbind, stat_results)
  stat_df$p_adj <- p.adjust(stat_df$p_value, method = "BH")
  
  # Prepare data for plotting
  plot_data <- rbind(
    data.frame(
      cellID = rownames(group1_patterns),
      group = group_labels[1],
      as.data.frame(group1_patterns),
      check.names = FALSE
    ),
    data.frame(
      cellID = rownames(group2_patterns),
      group = group_labels[2],
      as.data.frame(group2_patterns),
      check.names = FALSE
    )
  )
  
  return(list(
    comparison = comparison_df,
    stats = stat_df,
    plot_data = plot_data,
    group1_patterns = group1_patterns,
    group2_patterns = group2_patterns
  ))
}
## Pattern usage comparison holding a covariate fixed
# compare_pattern_usage() runs a bare Wilcoxon, which is exposed to any covariate
# that differs between the groups. For splice+ vs splice- that covariate is viral
# burden: a spliced UMI can only be seen in a cell that yields UMIs, so burden
# predicts splice status through detection. Host depth is left out of the model.
# Use this for a covariate that confounds the group label. For one that the
# group label causes, such as host library size in SIV+ vs SIV-, matching is
# the appropriate tool: see compare_pattern_usage_matched().
# This fits, per pattern,
#     usage ~ group + covariates
# and reports the group coefficient as an adjusted difference in mean usage,
# converted to a log2 fold change against the group2 mean so the value is on the
# same scale as the unadjusted row. `covariates` is a data frame of cell-level
# terms indexed by cellID; rows absent from it are dropped.
# Returns the same list shape as compare_pattern_usage(), so it can be passed to
# plot_pattern_enrichment() as another comparison.
compare_pattern_usage_adjusted <- function(usage, group1_cellids, group2_cellids = NULL,
                                           covariates,
                                           group_labels = c("infected", "uninfected")) {
  if (is.null(group2_cellids)) {
    group2_cellids <- setdiff(rownames(usage), group1_cellids)
  }
  cov_ok <- rownames(covariates)[stats::complete.cases(covariates)]
  g1 <- intersect(intersect(group1_cellids, rownames(usage)), cov_ok)
  g2 <- intersect(intersect(group2_cellids, rownames(usage)), cov_ok)
  cells <- c(g1, g2)
  if (length(g1) < 3 || length(g2) < 3) {
    stop("compare_pattern_usage_adjusted: fewer than 3 cells in a group after ",
         "intersecting with the covariate table")
  }
  U  <- usage[cells, , drop = FALSE]
  df <- data.frame(group = factor(c(rep(group_labels[1], length(g1)),
                                    rep(group_labels[2], length(g2))),
                                  levels = c(group_labels[2], group_labels[1])),
                   covariates[cells, , drop = FALSE], check.names = FALSE)
  cov_names <- colnames(covariates)

  g1_means <- colMeans(U[seq_along(g1), , drop = FALSE])
  g2_means <- colMeans(U[length(g1) + seq_along(g2), , drop = FALSE])

  res <- lapply(colnames(U), function(pat) {
    df$usage_value <- U[, pat]
    fit <- stats::lm(stats::as.formula(
      paste("usage_value ~ group +", paste(cov_names, collapse = " + "))), data = df)
    cf  <- summary(fit)$coefficients
    row <- grep("^group", rownames(cf))
    data.frame(pattern = pat,
               adjusted_difference = as.numeric(cf[row, 1]),
               p_value = as.numeric(cf[row, 4]),
               group1_median = stats::median(U[seq_along(g1), pat]),
               group2_median = stats::median(U[length(g1) + seq_along(g2), pat]))
  })
  stat_df <- do.call(rbind, res)
  stat_df$p_adj <- stats::p.adjust(stat_df$p_value, method = "BH")

  comparison_df <- data.frame(
    pattern     = colnames(U),
    group1_mean = g1_means,
    group2_mean = g2_means,
    log2_fold_change = log2((pmax(g2_means + stat_df$adjusted_difference, 1e-10)) /
                            (g2_means + 1e-10)))

  plot_data <- rbind(
    data.frame(cellID = g1, group = group_labels[1],
               as.data.frame(U[seq_along(g1), , drop = FALSE]), check.names = FALSE),
    data.frame(cellID = g2, group = group_labels[2],
               as.data.frame(U[length(g1) + seq_along(g2), , drop = FALSE]), check.names = FALSE))

  list(comparison = comparison_df, stats = stat_df, plot_data = plot_data,
       group1_patterns = U[seq_along(g1), , drop = FALSE],
       group2_patterns = U[length(g1) + seq_along(g2), , drop = FALSE],
       n_cells = c(group1 = length(g1), group2 = length(g2)))
}
## Pattern usage comparison against a matched reference group
# compare_pattern_usage_adjusted() removes a covariate by regression, which
# biases the group effect when the covariate lies on the causal path from the
# group label to the outcome rather than confounding it. Host library size is
# such a covariate for SIV+ vs SIV-: the size gap holds in the highest-burden
# cells, where a call cannot depend on sequencing depth, so it is a property of
# infected cells and not a detection artifact. cNMF usage is compositional
# (rows sum to 1), so depth does not inflate usage mechanically either, and each
# program carries some of the size axis itself. Regressing size out therefore
# removes signal that belongs to the contrast.
# This draws reference cells to match the group1 distribution of `match_values`
# over `n_bins` quantile bins, within each level of `strata` when given, and
# compares the two at equal depth without modelling depth. Up to `ratio`
# reference cells are drawn per group1 cell in a bin.
# Returns the same list shape as compare_pattern_usage(), plus n_cells.
compare_pattern_usage_matched <- function(usage, group1_cellids, group2_cellids = NULL,
                                          match_values, strata = NULL, n_bins = 10,
                                          ratio = 5, seed = 1,
                                          group_labels = c("infected", "uninfected")) {
  if (is.null(group2_cellids)) {
    group2_cellids <- setdiff(rownames(usage), group1_cellids)
  }
  match_values <- match_values[!is.na(match_values)]
  group1 <- intersect(intersect(group1_cellids, rownames(usage)), names(match_values))
  group2 <- intersect(intersect(group2_cellids, rownames(usage)), names(match_values))
  if (length(group1) < 3 || length(group2) < 3) {
    stop("compare_pattern_usage_matched: fewer than 3 cells in a group after ",
         "intersecting with the matching values")
  }
  cells  <- c(group1, group2)
  breaks <- unique(stats::quantile(match_values[cells],
                                   probs = seq(0, 1, length.out = n_bins + 1)))
  bins   <- cut(match_values[cells], breaks, include.lowest = TRUE)
  stratum <- if (is.null(strata)) {
    setNames(as.character(bins), cells)
  } else setNames(paste(strata[cells], bins), cells)

  wanted <- table(stratum[group1])
  set.seed(seed)
  matched <- unlist(lapply(names(wanted), function(s) {
    pool <- group2[stratum[group2] == s]
    pool[sample.int(length(pool), min(length(pool), ratio * wanted[[s]]))]
  }), use.names = FALSE)
  if (length(matched) < 3) {
    stop("compare_pattern_usage_matched: fewer than 3 reference cells matched")
  }

  result <- compare_pattern_usage(usage, group1, matched, group_labels = group_labels)
  result$n_cells <- c(group1 = length(group1), group2 = length(matched))
  result
}

### plotting
plot_usage_heatmap <- function(usage, metadata = NULL, output_path = NULL) {
  # Normalize usage per cell (row-wise)
  usage_norm <- usage / rowSums(usage)
  
  # Set up annotation if metadata provided
  if (!is.null(metadata)) {
    # Ensure rownames match
    common_cells <- intersect(rownames(usage_norm), rownames(metadata))
    usage_norm <- usage_norm[common_cells, ]
    annotation_row <- metadata[common_cells, , drop = FALSE]
  } else {
    annotation_row <- NA
  }
  
  # Create heatmap
  p <- pheatmap(
    t(usage_norm),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    color = viridis(100),
    main = "cNMF Pattern Usage (normalized)",
    annotation_col = if (!is.null(metadata)) annotation_row else NA,
    silent = TRUE
  )
  
  if (!is.null(output_path)) {
    pdf(output_path, width = 12, height = 6)
    print(p)
    dev.off()
    cat("Saved:", output_path, "\n")
  }
  
  return(p)
}
plot_marker_heatmap <- function(spectra, markers, n_genes = 20, output_path = NULL) {
  # Get top n genes for each pattern
  all_marker_genes <- unique(unlist(lapply(markers, head, n_genes)))
  
  # Subset spectra to marker genes
  spectra_subset <- spectra[all_marker_genes, ]
  
  # Z-score for visualization
  spectra_scaled <- t(scale(t(spectra_subset)))
  
  # Create annotation for which pattern each gene belongs to
  gene_pattern <- data.frame(
    pattern = sapply(all_marker_genes, function(g) {
      for (pat in names(markers)) {
        if (g %in% head(markers[[pat]], n_genes)) return(pat)
      }
      return(NA)
    }),
    row.names = all_marker_genes
  )
  
  p <- pheatmap(
    spectra_scaled,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = nrow(spectra_scaled) <= 50,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = paste("Top", n_genes, "Marker Genes per Pattern"),
    annotation_row = gene_pattern,
    silent = TRUE
  )
  
  if (!is.null(output_path)) {
    pdf(output_path, width = 8, height = 12)
    print(p)
    dev.off()
    cat("Saved:", output_path, "\n")
  }
  
  return(p)
}
plot_pattern_violins <- function(comparison_result, output_path = NULL,
                                 group_names = c("infected", "uninfected"),
                                 colors = c("#E41A1C", "gray70")) {
  plot_data <- comparison_result$plot_data
  stat_df <- comparison_result$stats
  
  plot_long <- pivot_longer(
    plot_data,
    cols = -c(cellID, group),
    names_to = "pattern",
    values_to = "weight"
  )
  
  # Relabel groups using provided names
  plot_long$group <- factor(plot_long$group,
                            levels = c(group_names[1], group_names[2]))
  
  max_weights <- plot_long %>%
    group_by(pattern) %>%
    summarise(ymax = max(weight), .groups = "drop")
  
  sig_data <- merge(stat_df, max_weights, by = "pattern")
  sig_data$y_pos <- sig_data$ymax * 1.05
  sig_data$label <- case_when(
    sig_data$p_adj < 0.001 ~ "***",
    sig_data$p_adj < 0.01 ~ "**",
    sig_data$p_adj < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  
  p <- ggplot(plot_long, aes(x = pattern, y = weight, fill = group)) +
    geom_violin(position = position_dodge(0.8), scale = "width", trim = TRUE) +
    geom_boxplot(
      position = position_dodge(0.8),
      width = 0.1,
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_text(
      data = sig_data,
      aes(x = pattern, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 4
    ) +
    scale_fill_manual(values = setNames(colors, group_names)) +
    labs(
      x = "Pattern",
      y = "Pattern Weight",
      fill = "Group",
      title = paste0("cNMF Pattern Usage: ", group_names[1], " vs ", group_names[2])
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  ggsave(output_path, p, width = 10, height = 6)
  
  return(p)
}
plot_log2fc <- function(comparison_result, output_path = NULL,
                        group_names = c("infected", "uninfected"),
                        colors = c("#E41A1C", "gray70")) {
  comparison_df <- comparison_result$comparison
  stat_df <- comparison_result$stats

  plot_df <- merge(comparison_df, stat_df[, c("pattern", "p_adj")], by = "pattern")

  pattern_order <- paste0("Pattern_", sort(as.numeric(gsub("Pattern_", "", plot_df$pattern))))
  plot_df$pattern <- factor(plot_df$pattern, levels = pattern_order)
  plot_df$row <- paste0(group_names[1], " / ", group_names[2])

  plot_df$stars <- cut(plot_df$p_adj,
                       breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("***", "**", "*", ""))
  plot_df$stars <- as.character(plot_df$stars)

  lim <- max(abs(plot_df$log2_fold_change), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1

  p <- ggplot(plot_df, aes(x = pattern, y = row, fill = log2_fold_change)) +
    geom_tile(color = "black", linewidth = 0.3) +
    geom_text(aes(label = stars), size = 4, vjust = 0.75) +
    scale_fill_gradient2(low = colors[2], mid = "white", high = colors[1],
                         midpoint = 0, limits = c(-lim, lim),
                         name = "log2 FC") +
    labs(x = NULL, y = NULL) +
    coord_fixed() +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 0))

  ggsave(output_path, p, width = max(4, 0.6 * nlevels(plot_df$pattern) + 2), height = 2)

  return(p)
}
# Multi-contrast cNMF pattern enrichment in a single plot object, so every
# contrast shares one pattern axis that lines up with the pattern-weight UMAP
# row above it. Contrasts keep their own log2 fold-change scale (the SIV+ and
# splice+ comparisons have very different dynamic ranges, and one shared scale
# flattens the smaller one).
#
# comparisons: list of lists, one per row, each with
#   result: a compare_pattern_usage() return value
#   label:  row label, e.g. "SIV+ vs SIV-"
#   colors: c(high, low) fill for enriched / depleted in group 1
# pattern_labels: optional named vector or character vector of program names
#   (e.g. "Naive & quiescent"), in Pattern_1..Pattern_k order.
# styles: any of "dots", "lollipop", "tiles"; each is written to
#   <output_prefix>_<style>.pdf and returned in the named list.
plot_pattern_enrichment <- function(comparisons,
                                    output_prefix = NULL,
                                    styles = c("dots", "lollipop", "tiles"),
                                    pattern_labels = NULL,
                                    show_n = TRUE,
                                    max_neglog10_p = 10,
                                    drop_pattern_labels = "contaminant|contamination",
                                    fc_limits = NULL) {
  styles <- match.arg(styles, several.ok = TRUE)

  # One tidy row per (contrast, pattern) with the fold change, the adjusted
  # p-value and the per-contrast fold-change limit used for that row's scale.
  enrichment <- dplyr::bind_rows(lapply(comparisons, function(entry) {
    comparison_df <- entry$result$comparison
    stat_df       <- entry$result$stats
    n_group1      <- nrow(entry$result$group1_patterns)
    n_group2      <- nrow(entry$result$group2_patterns)
    row_label     <- if (show_n) {
      paste0(entry$label, "\n(", n_group1, " vs ", n_group2, " cells)")
    } else entry$label
    merge(comparison_df, stat_df[, c("pattern", "p_adj")], by = "pattern") %>%
      dplyr::mutate(contrast   = row_label,
                    color_high = entry$colors[1],
                    color_low  = entry$colors[2])
  }))

  pattern_levels <- paste0("Pattern_",
                           sort(unique(as.numeric(gsub("Pattern_", "",
                                                       enrichment$pattern)))))
  if (is.null(pattern_labels)) {
    pattern_display <- setNames(gsub("Pattern_", "P", pattern_levels), pattern_levels)
  } else if (is.null(names(pattern_labels))) {
    pattern_display <- setNames(pattern_labels[seq_along(pattern_levels)], pattern_levels)
  } else {
    pattern_display <- pattern_labels[pattern_levels]
    names(pattern_display) <- pattern_levels
    pattern_display[is.na(pattern_display)] <- gsub("Pattern_", "P",
                                                    pattern_levels[is.na(pattern_display)])
  }

  # Programs the annotation marks as contaminants are dropped before anything
  # is scaled. They are not the biology the panel is about, and one of them can
  # carry a degenerate fold change (a program with near-zero usage in one group)
  # that sets the colour limits for its whole row and flattens every real
  # difference. Their numbers are still printed by the calling chunk.
  if (!is.null(drop_pattern_labels) && nzchar(drop_pattern_labels)) {
    drop <- grepl(drop_pattern_labels, pattern_display, ignore.case = TRUE)
    if (any(drop)) {
      cat("Dropping contaminant program(s) from the panel: ",
          paste(pattern_display[drop], collapse = "; "), "\n", sep = "")
      keep_patterns  <- names(pattern_display)[!drop]
      pattern_levels <- pattern_levels[pattern_levels %in% keep_patterns]
      pattern_display <- pattern_display[pattern_levels]
      enrichment <- enrichment[enrichment$pattern %in% pattern_levels, ]
    }
  }

  # A row label minus its trailing "(n1 vs n2 cells)" line, on one line.
  short_label <- function(x) {
    gsub("\n", " ", sub("\n\\([0-9]+ vs [0-9]+ cells\\)$", "", x))
  }

  contrast_levels <- unique(enrichment$contrast)
  enrichment <- enrichment %>%
    dplyr::mutate(
      pattern     = factor(pattern, levels = pattern_levels),
      program     = factor(unname(pattern_display[as.character(pattern)]),
                           levels = unname(pattern_display[pattern_levels])),
      contrast    = factor(contrast, levels = contrast_levels),
      significant = p_adj < 0.05,
      stars       = as.character(cut(p_adj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                     labels = c("***", "**", "*", ""))),
      neglog10_p  = pmin(-log10(pmax(p_adj, 1e-300)), max_neglog10_p)) %>%
    dplyr::mutate(scale_group = paste(color_high, color_low, sep = "|")) %>%
    # Rows drawn in the same colours are the same contrast before and after a
    # control, so they share one limit: a dot that pales after adjustment has
    # actually shrunk, rather than having been rescaled. `fc_limits` overrides
    # the computed limit (named by the group's legend title, or one number for
    # every group), which is how several panels are put on one scale.
    dplyr::group_by(scale_group) %>%
    dplyr::mutate(fc_limit = {
      lim <- max(abs(log2_fold_change), na.rm = TRUE)
      if (!is.finite(lim) || lim == 0) 1 else lim
    }) %>%
    dplyr::ungroup()

  group_levels <- unique(enrichment$scale_group)
  group_title  <- vapply(group_levels, function(g)
    paste0("log2 FC, ", short_label(as.character(
      enrichment$contrast[enrichment$scale_group == g][1]))),
    character(1))
  if (!is.null(fc_limits)) {
    # Names are matched on the legend title, and failing that on the title with
    # its trailing "(... controlled)" qualifier removed, so one set of limits
    # covers a panel drawn from every row and a panel drawn from the controlled
    # rows alone, where the group's first row carries the qualifier.
    family <- function(x) trimws(sub("\\s*\\([^()]*\\)$", "", x))
    wanted <- if (is.null(names(fc_limits))) {
      rep(unname(fc_limits)[1], length(group_levels))
    } else {
      hit <- match(group_title, names(fc_limits))
      hit[is.na(hit)] <- match(family(group_title), family(names(fc_limits)))[is.na(hit)]
      unname(fc_limits[hit])
    }
    for (i in seq_along(group_levels)) {
      if (!is.na(wanted[i]) && wanted[i] > 0)
        enrichment$fc_limit[enrichment$scale_group == group_levels[i]] <- wanted[i]
    }
  }

  wrap_labels <- scales::label_wrap(14)
  base_theme <- ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.text.x        = ggplot2::element_text(size = 9, lineheight = 0.95),
      axis.text.y        = ggplot2::element_text(size = 9, lineheight = 0.95, hjust = 1),
      axis.ticks         = ggplot2::element_blank(),
      strip.text.y.left  = ggplot2::element_text(angle = 0, hjust = 1, size = 9,
                                                 lineheight = 0.95),
      legend.key.height  = grid::unit(0.3, "cm"),
      legend.key.width   = grid::unit(0.8, "cm"),
      legend.position    = "bottom",
      legend.box         = "horizontal",
      legend.title       = ggplot2::element_text(size = 8),
      legend.text        = ggplot2::element_text(size = 7),
      plot.margin        = ggplot2::margin(6, 10, 4, 6))

  # Legend titles sit above their key so they cannot collide with the tick
  # labels of the neighbouring legend. Each row carries its own fold-change
  # scale, so the title names the row in full: dropping the qualifier gives two
  # rows the same title and leaves no way to tell which colorbar belongs to
  # which line. `order` fixes the legends in row order, top row first.
  colorbar_guide <- function(i) ggplot2::guide_colorbar(title.position = "top",
                                                        title.hjust = 0, order = i)

  n_patterns <- length(pattern_levels)
  n_contrasts <- length(contrast_levels)
  plots <- list()

  # Dots: fill = log2 fold change (own scale per contrast), size = -log10 padj,
  # so non-significant programs collapse to a speck instead of a full tile.
  if ("dots" %in% styles) {
    p_dots <- ggplot2::ggplot(enrichment,
                              ggplot2::aes(x = program, y = contrast))
    for (i in seq_along(group_levels)) {
      row_data <- enrichment[enrichment$scale_group == group_levels[i], ]
      p_dots <- p_dots +
        ggplot2::geom_point(data = row_data,
                            ggplot2::aes(size = neglog10_p, fill = log2_fold_change),
                            shape = 21, color = "gray25", stroke = 0.3) +
        ggplot2::scale_fill_gradient2(
          low = row_data$color_low[1], mid = "white", high = row_data$color_high[1],
          midpoint = 0, limits = c(-row_data$fc_limit[1], row_data$fc_limit[1]),
          breaks = c(-row_data$fc_limit[1], 0, row_data$fc_limit[1]),
          labels = function(x) sprintf("%.1f", x),
          guide = colorbar_guide(i), name = group_title[i]) +
        ggnewscale::new_scale_fill()
    }
    p_dots <- p_dots +
      ggplot2::scale_size_continuous(range = c(1, 9), limits = c(0, max_neglog10_p),
                                     breaks = c(1.3, 5, 10),
                                     labels = c("0.05", "1e-5", "<=1e-10"),
                                     name = "adj. p",
                                     guide = ggplot2::guide_legend(title.position = "top",
                                                                   order = length(group_levels) + 1)) +
      ggplot2::scale_x_discrete(labels = wrap_labels) +
      ggplot2::scale_y_discrete(limits = rev(contrast_levels), labels = wrap_labels) +
      ggplot2::labs(x = NULL, y = NULL) +
      base_theme +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "gray92"))
    plots$dots <- p_dots
  }

  # Lollipop: fold change on a real axis (free per contrast), filled point =
  # adj. p < 0.05, hollow = not significant. No color legend needed.
  if ("lollipop" %in% styles) {
    lolli_data <- enrichment %>%
      dplyr::mutate(point_color = ifelse(log2_fold_change >= 0, color_high, color_low),
                    point_fill  = ifelse(significant, point_color, "white"))
    p_lolli <- ggplot2::ggplot(lolli_data,
                               ggplot2::aes(x = program, y = log2_fold_change)) +
      ggplot2::geom_hline(yintercept = 0, color = "gray60", linewidth = 0.3) +
      ggplot2::geom_segment(ggplot2::aes(xend = program, yend = 0, color = point_color),
                            linewidth = 0.6) +
      ggplot2::geom_point(ggplot2::aes(color = point_color, fill = point_fill),
                          shape = 21, size = 3, stroke = 0.8) +
      ggplot2::scale_color_identity() +
      ggplot2::scale_fill_identity() +
      ggplot2::facet_grid(contrast ~ ., scales = "free_y", switch = "y",
                          labeller = ggplot2::labeller(contrast = wrap_labels)) +
      ggplot2::scale_x_discrete(labels = wrap_labels) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.18)) +
      ggplot2::labs(x = NULL, y = "log2 fold change",
                    caption = "filled = adj. p < 0.05") +
      base_theme +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "gray92"),
                     panel.background = ggplot2::element_rect(fill = "gray98", color = NA),
                     panel.spacing.y = grid::unit(0.5, "lines"),
                     strip.placement = "outside",
                     plot.caption = ggplot2::element_text(size = 7, hjust = 0))
    plots$lollipop <- p_lolli
  }

  # Tiles: the current look, cleaned up -- one plot object, no tile borders,
  # own fill scale per contrast, fold change printed in the tile and stars only
  # where significant.
  if ("tiles" %in% styles) {
    p_tiles <- ggplot2::ggplot(enrichment,
                               ggplot2::aes(x = program, y = contrast))
    for (i in seq_along(group_levels)) {
      row_data <- enrichment[enrichment$scale_group == group_levels[i], ]
      p_tiles <- p_tiles +
        ggplot2::geom_tile(data = row_data,
                           ggplot2::aes(fill = log2_fold_change, alpha = significant),
                           width = 0.96, height = 0.96) +
        ggplot2::scale_fill_gradient2(
          low = row_data$color_low[1], mid = "white", high = row_data$color_high[1],
          midpoint = 0, limits = c(-row_data$fc_limit[1], row_data$fc_limit[1]),
          breaks = c(-row_data$fc_limit[1], 0, row_data$fc_limit[1]),
          labels = function(x) sprintf("%.1f", x),
          guide = colorbar_guide(i), name = group_title[i]) +
        ggnewscale::new_scale_fill()
    }
    p_tiles <- p_tiles +
      ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 1),
                                  guide = "none") +
      ggplot2::geom_text(ggplot2::aes(label = ifelse(significant,
                                                     paste0(sprintf("%.2f", log2_fold_change),
                                                            "\n", stars),
                                                     "")),
                         size = 2.8, lineheight = 0.85, color = "gray10") +
      ggplot2::scale_x_discrete(labels = wrap_labels) +
      ggplot2::scale_y_discrete(limits = rev(contrast_levels), labels = wrap_labels) +
      ggplot2::labs(x = NULL, y = NULL) +
      base_theme +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    plots$tiles <- p_tiles
  }

  if (!is.null(output_prefix)) {
    dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
    plot_width  <- max(5, 1.05 * n_patterns + 2.2)
    for (style in names(plots)) {
      plot_height <- switch(style,
                            lollipop = 1.9 * n_contrasts + 1.6,
                            0.7 * n_contrasts + 1.9)
      output_path <- paste0(output_prefix, "_", style, ".pdf")
      ggplot2::ggsave(output_path, plots[[style]],
                      width = plot_width, height = plot_height)
      cat("Saved:", output_path, "\n")
    }
  }

  return(plots)
}
# =============================================================================
# HCR
# =============================================================================
## Q-Q plot comparing method distributions to a reference distribution
plot_qq_comparison <- function(
    method_list,
    ref_vec,
    var_name = "Variable",
    n_boot = 1000,
    seed = 123
) {
  set.seed(seed)

  method_names <- names(method_list)

  probs <- seq(0.01, 0.99, by = 0.01)
  q_ref <- quantile(ref_vec, probs)
  q_methods <- lapply(method_list, function(vec) quantile(vec, probs))

  # Bootstrap CI for reference (null expectation)
  n_boot_size <- max(sapply(method_list, length))
  boot_q <- matrix(NA, nrow = n_boot, ncol = length(probs))

  for (i in seq_len(n_boot)) {
    boot_samp <- sample(ref_vec, n_boot_size, replace = TRUE)
    boot_q[i, ] <- quantile(boot_samp, probs)
  }

  ci_lower <- apply(boot_q, 2, quantile, probs = 0.025)
  ci_upper <- apply(boot_q, 2, quantile, probs = 0.975)

  dev_stats <- lapply(method_names, function(nm) {
    dev <- q_methods[[nm]] - q_ref
    c(mad = mean(abs(dev)), max_dev = max(abs(dev)))
  })
  names(dev_stats) <- method_names

  subtitle_parts <- sapply(method_names, function(nm) {
    sprintf("%s: MAD=%.1f, Max=%.1f", nm, dev_stats[[nm]]["mad"], dev_stats[[nm]]["max_dev"])
  })
  subtitle_text <- paste(subtitle_parts, collapse = " | ")

  df_list <- lapply(method_names, function(nm) {
    data.frame(
      ref = q_ref,
      exp = q_methods[[nm]],
      method = nm
    )
  })
  df <- do.call(rbind, df_list)
  df$method <- factor(df$method, levels = method_names)

  ci_df <- data.frame(ref = q_ref, ci_lower = ci_lower, ci_upper = ci_upper)

  # Figure colour key used across the paper: Ca*Re light blue, TenX/10X dark blue,
  # longread grey. HCR is the reference and is drawn as a dashed grey diagonal.
  method_colors <- c("Ca*Re" = "#00AEEF", "TenX" = "#2B3990", "longread" = "gray55", "HCR" = "gray30")
  plot_colors <- method_colors[method_names]

  all_exp <- unlist(q_methods)

  ggplot() +
    geom_ribbon(data = ci_df, aes(x = ref, ymin = ci_lower, ymax = ci_upper),
                fill = "grey70", alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_line(data = df, aes(x = ref, y = exp, color = method), linewidth = 1) +
    scale_color_manual(values = plot_colors) +
    coord_cartesian(
      xlim = c(0, max(q_ref) * 1.05),
      ylim = c(0, max(all_exp) * 1.05)
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      plot.subtitle = element_text(size = 9)
    ) +
    labs(
      x = "HCR quantiles",
      y = "Method quantiles",
      title = var_name,
      subtitle = subtitle_text,
      color = NULL
    )
}

## Q-Q plot grid across multiple variables
plot_qq_grid <- function(
    method_data,
    ref_data,
    var_names = c("US", "spliced", "SS", "MS"),
    method_names = c("Ca*Re", "TenX", "longread"),
    output_file = NULL,
    width = 16,
    height = 4,
    n_boot = 1000,
    seed = 123
) {
  plot_list <- list()

  for (var in var_names) {
    method_list <- setNames(
      lapply(method_names, function(m) method_data[[m]][[var]]),
      method_names
    )

    p <- plot_qq_comparison(
      method_list = method_list,
      ref_vec = ref_data[[var]],
      var_name = var,
      n_boot = n_boot,
      seed = seed
    )

    plot_list[[var]] <- p
  }

  combined <- gridExtra::arrangeGrob(grobs = plot_list, nrow = 1)

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, combined, width = width, height = height)
    message("Saved: ", output_file)
  }

  grid::grid.draw(combined)
  invisible(plot_list)
}

## Compute KS, JS, and Wasserstein distances for each method vs reference,
## with bootstrap CIs and permutation p-values (vs-ref and pairwise between methods)
compare_methods_to_reference <- function(
    method_list,
    ref_vec,
    var_name = "Variable",
    n_boot = 2000,
    n_perm = 2000,
    seed = 123
) {
  set.seed(seed)

  method_names <- names(method_list)
  n_methods <- length(method_list)

  for (nm in method_names) {
    stopifnot(all(method_list[[nm]] >= 0))
  }
  stopifnot(all(ref_vec >= 0))

  normalize <- function(x, max_val) {
    tab <- table(factor(x, levels = 0:max_val))
    (as.numeric(tab) + 1e-10) / (sum(tab) + 1e-10 * (max_val + 1))
  }

  # General form, for inputs that are not counts. Kept as the reference
  # implementation and as the fallback for non-integer input.
  compute_distances_general <- function(vec1, vec2) {
    max_val <- max(c(vec1, vec2))
    p1 <- normalize(vec1, max_val)
    p2 <- normalize(vec2, max_val)

    ks_stat <- suppressWarnings(ks.test(vec1, vec2)$statistic)
    js_val <- sqrt(philentropy::JSD(rbind(p1, p2), unit = "log2", est.prob = NULL))

    result <- c(as.numeric(ks_stat), as.numeric(js_val), transport::wasserstein1d(vec1, vec2))
    names(result) <- c("ks", "js", "wasserstein")
    result
  }

  # All three metrics are functions of the two empirical CDFs, so for
  # non-negative integer counts one pair of tabulated histograms serves all of
  # them, in closed form and with no sorting:
  #   KS = max |F1 - F2|
  #   W1 = sum |F1 - F2|   (support is the integers, so each bin has width 1)
  #   JS = sqrt of the Jensen-Shannon divergence, log2
  # This replaces table(factor(...)), ks.test, philentropy::JSD and
  # transport::wasserstein1d, and agrees with compute_distances_general to
  # floating-point tolerance. The bootstrap and permutation loops call it tens of
  # thousands of times per variable, which is why the closed form is used here.
  compute_distances <- function(vec1, vec2) {
    is_count <- function(v) {
      is.numeric(v) && length(v) > 0 && all(is.finite(v)) &&
        all(v >= 0) && all(v == trunc(v)) && max(v) < 1e7
    }
    if (!is_count(vec1) || !is_count(vec2)) {
      return(compute_distances_general(vec1, vec2))
    }

    n_bins <- max(c(vec1, vec2)) + 1L
    c1 <- tabulate(vec1 + 1L, nbins = n_bins)
    c2 <- tabulate(vec2 + 1L, nbins = n_bins)
    n1 <- length(vec1)
    n2 <- length(vec2)

    # same 1e-10 smoothing as normalize(), so JS matches the general form
    p1 <- (c1 + 1e-10) / (n1 + 1e-10 * n_bins)
    p2 <- (c2 + 1e-10) / (n2 + 1e-10 * n_bins)
    f1 <- cumsum(c1) / n1
    f2 <- cumsum(c2) / n2
    cdf_gap <- abs(f1 - f2)

    m <- (p1 + p2) / 2
    kl <- function(p, q) {
      keep <- p > 0
      sum(p[keep] * log2(p[keep] / q[keep]))
    }

    result <- c(max(cdf_gap),
                sqrt(0.5 * kl(p1, m) + 0.5 * kl(p2, m)),
                sum(cdf_gap))
    names(result) <- c("ks", "js", "wasserstein")
    result
  }

  # H0: two samples come from same distribution
  perm_test_distribution <- function(vec1, vec2, obs_dist, n_perm) {
    pooled <- c(vec1, vec2)
    n1 <- length(vec1)
    n2 <- length(vec2)

    perm_stats <- matrix(NA, nrow = n_perm, ncol = 3)

    for (i in seq_len(n_perm)) {
      perm_idx <- sample(length(pooled))
      perm1 <- pooled[perm_idx[1:n1]]
      perm2 <- pooled[perm_idx[(n1 + 1):(n1 + n2)]]
      perm_stats[i, ] <- compute_distances(perm1, perm2)
    }

    # One-sided p-value (Phipson & Smyth 2010)
    sapply(seq_len(3), function(j) {
      (sum(perm_stats[, j] >= obs_dist[j]) + 1) / (n_perm + 1)
    })
  }

  # H0: method1 and method2 are equally distant from reference
  perm_test_distance_difference <- function(vec1, vec2, ref_vec, n_perm) {
    n1 <- length(vec1)
    n2 <- length(vec2)

    obs_d1 <- compute_distances(vec1, ref_vec)
    obs_d2 <- compute_distances(vec2, ref_vec)
    obs_diff <- obs_d1 - obs_d2

    pooled <- c(vec1, vec2)
    perm_diffs <- matrix(NA, nrow = n_perm, ncol = 3)

    for (i in seq_len(n_perm)) {
      perm_idx <- sample(length(pooled))
      perm1 <- pooled[perm_idx[1:n1]]
      perm2 <- pooled[perm_idx[(n1 + 1):(n1 + n2)]]

      d1 <- compute_distances(perm1, ref_vec)
      d2 <- compute_distances(perm2, ref_vec)
      perm_diffs[i, ] <- d1 - d2
    }

    pvals <- sapply(seq_len(3), function(j) {
      (sum(abs(perm_diffs[, j]) >= abs(obs_diff[j])) + 1) / (n_perm + 1)
    })
    names(pvals) <- c("ks", "js", "wasserstein")

    list(
      obs_diff = obs_diff,
      pval = pvals,
      perm_diffs = perm_diffs
    )
  }

  boot_ci <- function(vec1, vec2, n_boot) {
    n1 <- length(vec1)
    n2 <- length(vec2)

    boot_stats <- matrix(NA, nrow = n_boot, ncol = 3)

    for (i in seq_len(n_boot)) {
      b1 <- sample(vec1, n1, replace = TRUE)
      b2 <- sample(vec2, n2, replace = TRUE)
      boot_stats[i, ] <- compute_distances(b1, b2)
    }

    list(
      ci_lower = apply(boot_stats, 2, quantile, probs = 0.025),
      ci_upper = apply(boot_stats, 2, quantile, probs = 0.975),
      se = apply(boot_stats, 2, sd)
    )
  }

  dist_to_ref <- list()
  ci_to_ref <- list()
  pval_to_ref <- list()

  for (nm in method_names) {
    dist_to_ref[[nm]] <- compute_distances(method_list[[nm]], ref_vec)
    ci_to_ref[[nm]] <- boot_ci(method_list[[nm]], ref_vec, n_boot)
    pval_to_ref[[nm]] <- perm_test_distribution(method_list[[nm]], ref_vec, dist_to_ref[[nm]], n_perm)
    names(pval_to_ref[[nm]]) <- c("ks", "js", "wasserstein")
  }

  paired_distance_tests <- list()

  for (i in seq_len(n_methods - 1)) {
    for (j in (i + 1):n_methods) {
      nm1 <- method_names[i]
      nm2 <- method_names[j]
      pair_key <- paste0(nm1, "_vs_", nm2)

      test_result <- perm_test_distance_difference(
        method_list[[nm1]],
        method_list[[nm2]],
        ref_vec,
        n_perm
      )

      closer <- ifelse(test_result$obs_diff < 0, nm1, nm2)
      names(closer) <- c("ks", "js", "wasserstein")

      paired_distance_tests[[pair_key]] <- list(
        obs_diff = test_result$obs_diff,
        pval = test_result$pval,
        closer = closer
      )
    }
  }

  metrics <- c("KS", "JS", "Wasserstein")

  distance_df <- data.frame(
    variable = var_name,
    metric = rep(metrics, each = n_methods + 1),
    method = rep(c(method_names, "HCR"), 3),
    distance = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    pval_vs_ref = NA_real_
  )

  for (m in seq_along(metrics)) {
    metric_key <- c("ks", "js", "wasserstein")[m]
    for (i in seq_along(method_names)) {
      nm <- method_names[i]
      row_idx <- (m - 1) * (n_methods + 1) + i
      distance_df$distance[row_idx] <- dist_to_ref[[nm]][metric_key]
      distance_df$ci_lower[row_idx] <- ci_to_ref[[nm]]$ci_lower[metric_key]
      distance_df$ci_upper[row_idx] <- ci_to_ref[[nm]]$ci_upper[metric_key]
      distance_df$pval_vs_ref[row_idx] <- pval_to_ref[[nm]][metric_key]
    }
    ref_row_idx <- (m - 1) * (n_methods + 1) + n_methods + 1
    distance_df$distance[ref_row_idx] <- 0
    distance_df$ci_lower[ref_row_idx] <- 0
    distance_df$ci_upper[ref_row_idx] <- 0
    distance_df$pval_vs_ref[ref_row_idx] <- NA
  }

  list(
    variable = var_name,
    method_names = method_names,
    distance_df = distance_df,
    dist_to_ref = dist_to_ref,
    pval_to_ref = pval_to_ref,
    paired_distance_tests = paired_distance_tests
  )
}

## Per-variable Wasserstein bar panels with pairwise significance brackets
plot_wasserstein_bars <- function(
    results_list,
    var_names = c("US", "spliced", "SS", "MS"),
    var_labels = c("unspliced", "spliced", "singly spliced", "multiply spliced"),
    method_names = c("Ca*Re", "TenX", "longread"),
    output_file = NULL,
    width = 14,
    height = 4
) {
  format_sig <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    return("n.s.")
  }

  # Same key as plot_qq_comparison.
  method_colors <- c("Ca*Re" = "#00AEEF", "TenX" = "#2B3990", "longread" = "gray55")

  get_paired_pval <- function(result, m1, m2) {
    pair_key <- paste0(m1, "_vs_", m2)
    alt_key <- paste0(m2, "_vs_", m1)
    if (pair_key %in% names(result$paired_distance_tests)) {
      return(result$paired_distance_tests[[pair_key]]$pval["wasserstein"])
    } else if (alt_key %in% names(result$paired_distance_tests)) {
      return(result$paired_distance_tests[[alt_key]]$pval["wasserstein"])
    }
    return(NA)
  }

  plot_list <- lapply(seq_along(var_names), function(i) {
    var <- var_names[i]
    var_label <- var_labels[i]
    result <- results_list[[i]]

    df <- result$distance_df[result$distance_df$metric == "Wasserstein" &
                              result$distance_df$method %in% method_names, ]
    df$method <- factor(df$method, levels = method_names)

    y_max <- max(df$distance, na.rm = TRUE)
    y_step <- y_max * 0.12
    y_base <- y_max * 1.08

    p <- ggplot(df, aes(x = method, y = distance, fill = method)) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_manual(values = method_colors) +
      labs(title = var_label, x = NULL, y = "Wasserstein distance") +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "none",
        axis.text.x = element_text(size = 10)
      )

    bracket_count <- 0
    # Method vs HCR
    for (m_idx in seq_along(method_names)) {
      m <- method_names[m_idx]
      pval <- result$pval_to_ref[[m]]["wasserstein"]
      sig <- format_sig(pval)
      x_pos <- m_idx
      y_pos <- df$distance[df$method == m] * 0.5
      p <- p + annotate("text", x = x_pos, y = y_pos, label = sig,
                       size = 3.5, color = "white", fontface = "bold")
    }

    # Pairwise between methods
    for (m1_idx in seq_len(length(method_names) - 1)) {
      for (m2_idx in (m1_idx + 1):length(method_names)) {
        m1 <- method_names[m1_idx]
        m2 <- method_names[m2_idx]
        pval <- get_paired_pval(result, m1, m2)
        sig <- format_sig(pval)

        y_bracket <- y_base + bracket_count * y_step
        p <- p +
          annotate("segment", x = m1_idx, xend = m1_idx,
                   y = y_bracket, yend = y_bracket + y_step * 0.2, linewidth = 0.3) +
          annotate("segment", x = m1_idx, xend = m2_idx,
                   y = y_bracket + y_step * 0.2, yend = y_bracket + y_step * 0.2, linewidth = 0.3) +
          annotate("segment", x = m2_idx, xend = m2_idx,
                   y = y_bracket, yend = y_bracket + y_step * 0.2, linewidth = 0.3) +
          annotate("text", x = (m1_idx + m2_idx) / 2, y = y_bracket + y_step * 0.4,
                   label = sig, size = 3.2)
        bracket_count <- bracket_count + 1
      }
    }

    p + coord_cartesian(ylim = c(0, y_max * 1.6))
  })

  combined <- gridExtra::arrangeGrob(grobs = plot_list, nrow = 1)

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, combined, width = width, height = height)
    message("Saved: ", output_file)
  }

  grid::grid.draw(combined)
  invisible(plot_list)
}

# =============================================================================
# PTPRC / CD45 per-cell isoform classifier
# =============================================================================

# Assign a single CD45 isoform per cell from per-UMI read counts
# (<gene>_UMI_read_counts_full.csv with cellID, UMI, category, read_count).
# Each UMI's variable exons A/B/C (exon4/5/6) are scored as included (an
# inclusion read), excluded (a skip-junction read), or undetermined; a UMI is
# labelled by its included exons (RA/RB/RC/RAB/.../RABC) and is "fully
# resolved" when all three exons are determined. RO is reserved for molecules
# shown to exclude all three cassettes, matching the UCHL1 epitope, which only
# exists on the isoform lacking A, B and C. A molecule with no inclusion but
# incomplete exclusion evidence is merely RA-negative -- it cannot name an
# isoform, and is dropped rather than counted as RO.
# Per cell:
#   Rule 1: if any UMI is fully resolved, take the most abundant fully-resolved
#           isoform (the molecule carrying the most exon information wins).
#   Rule 2: otherwise take the most abundant observed-inclusion label.
# Returns a data.frame: cellID, isoform, n_umi, resolved (TRUE = Rule 1).
classify_ptprc_isoforms <- function(umi_file, min_reads_per_umi = 2) {
  reads <- utils::read.csv(umi_file)
  reads <- reads[reads$read_count >= min_reads_per_umi, , drop = FALSE]

  per_umi <- reads %>%
    dplyr::group_by(cellID, UMI) %>%
    dplyr::summarise(
      inc_A = any(category == "PTPRC_exon4"),
      inc_B = any(category == "PTPRC_exon5"),
      inc_C = any(category == "PTPRC_exon6"),
      exc_A = any(category == "PTPRC_exon4_exclusion"),
      exc_B = any(category == "PTPRC_exon5_exclusion"),
      exc_C = any(category == "PTPRC_exon6_exclusion"),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      n_det = (inc_A | exc_A) + (inc_B | exc_B) + (inc_C | exc_C),
      label = paste0("R", ifelse(inc_A, "A", ""), ifelse(inc_B, "B", ""),
                     ifelse(inc_C, "C", "")),
      label = ifelse(label == "R",
                     ifelse(exc_A & exc_B & exc_C, "RO", NA_character_),
                     label)
    ) %>%
    dplyr::filter(n_det > 0, !is.na(label))

  # most frequent label, ties broken by table()'s ordering
  top_label <- function(labels) names(sort(table(labels), decreasing = TRUE))[1]

  per_umi %>%
    dplyr::group_by(cellID) %>%
    dplyr::summarise(
      n_umi    = dplyr::n(),
      resolved = max(n_det) == 3,
      isoform  = if (max(n_det) == 3) top_label(label[n_det == 3]) else top_label(label),
      .groups  = "drop"
    )
}


# =============================================================================
# Public human HIV datasets (Collora et al. 2022, HIV-seq/Frouard 2026)
# =============================================================================

# region -> base isoform category (regular bamsort logic)
hiv_region_category <- c(LTR_D1 = "any", D1_A1 = "US", A1_D4 = "any",
                         D4_A7 = "US.SS", A7_LTR = "any")
# per-UMI isoform priority (most specific wins). SS is absent because no region
# in hiv_region_category yields it and the splice override emits only MS/SS.MS,
# so a UMI never carries it and it needs no priority.
hiv_cat_rank   <- c(MS = 1, SS.MS = 2, US = 3, US.SS = 4, any = 5)
# The same six categories as CaRe/TenX, so the public-data arms bind against the
# SIV tables without a column mismatch. SS is unreachable here, as above, and
# comes back zero-filled.
hiv_categories <- c("US", "US.SS", "SS.MS", "SS", "MS", "any")
hiv_D4_donor   <- 6044   # HXB2 D4 donor, 0-based (bamsort_splice Donor col) -> MS
hiv_donor_tol  <- 10

# Per-cell HXB2 isoform quantification
hxb2_per_cell <- function(sample, bamsort_dir) {
  # Built from hiv_categories rather than listed out, so the no-reads path
  # returns the same columns as the populated one.
  empty <- tibble(cellID = character(), hxb2 = numeric(),
                  hxb2_umis = numeric(), hxb2_top_umi = numeric())
  empty[hiv_categories] <- lapply(hiv_categories, function(x) numeric())
  dfs <- lapply(names(hiv_region_category), function(rg) {
    f <- file.path(bamsort_dir, paste0(sample, "_bamsort_alignment_", rg, ".csv"))
    if (!file.exists(f)) return(NULL)
    d <- read.csv(f, colClasses = c(UMI = "character", cellID = "character",
                                    count = "integer"))
    if (nrow(d) == 0) return(NULL)
    # Adjacent junction windows mean one read can clear min_overlap in two of
    # them and be written to both files, so summing `count` over regions counts
    # it twice. `primary_count` charges each read to the single window it
    # overlaps most. Tables written before that column existed fall back to
    # `count`, matching `bamsort_read_column` on the pipeline side.
    if ("primary_count" %in% names(d)) {
      d$count <- d$primary_count
    } else {
      warning("bamsort table predates `primary_count`, so reads spanning two ",
              "junction windows are counted twice: ", basename(f),
              ". Re-run stepTwo (its --windows flag writes the column).",
              call. = FALSE)
    }
    d$category <- unname(hiv_region_category[rg]); d
  })
  df <- dplyr::bind_rows(dfs)
  if (nrow(df) == 0) return(empty)

  # Per-UMI read counts, which both papers' HIV+ rules are written in terms of:
  # `hxb2` is their sum, `hxb2_umis` the molecule count and `hxb2_top_umi` the
  # reads of the best-covered molecule. See hiv_pos_call.
  per_umi <- df %>% group_by(cellID, UMI) %>%
    summarize(reads = sum(count), .groups = "drop")
  total <- per_umi %>% group_by(cellID) %>%
    summarize(hxb2 = sum(reads), hxb2_umis = n(), hxb2_top_umi = max(reads),
              .groups = "drop")

  # stage 1: one region category per (cell, UMI), highest priority wins
  umi_cat <- df %>%
    distinct(cellID, UMI, category) %>%
    mutate(rank = unname(hiv_cat_rank[category])) %>%
    group_by(cellID, UMI) %>%
    slice_min(rank, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(cellID, UMI, category)

  # stage 2: splice override (bamsort_splice donor -> MS / SS.MS)
  sf <- file.path(bamsort_dir, "splice", paste0(sample, "_splicesites.csv"))
  if (file.exists(sf)) {
    sp <- read.csv(sf, colClasses = c(CB = "character", UB = "character",
                                      Donor = "integer"))
    if (nrow(sp) > 0) {
      splice_cat <- sp %>%
        transmute(cellID = CB, UMI = UB,
                  is_ms = abs(Donor - hiv_D4_donor) <= hiv_donor_tol) %>%
        group_by(cellID, UMI) %>%
        summarize(splice = if (any(is_ms)) "MS" else "SS.MS", .groups = "drop")
      umi_cat <- umi_cat %>%
        left_join(splice_cat, by = c("cellID", "UMI")) %>%
        mutate(category = ifelse(is.na(splice), category, splice)) %>%
        select(cellID, UMI, category)
    }
  }

  cats <- umi_cat %>%
    count(cellID, category, name = "n") %>%
    tidyr::pivot_wider(names_from = category, values_from = n, values_fill = 0)
  out <- left_join(total, cats, by = "cellID") %>%
    mutate(cellID = paste0(sample, "_", cellID))
  for (cc in hiv_categories) if (!cc %in% names(out)) out[[cc]] <- 0
  out %>% select(cellID, hxb2, hxb2_umis, hxb2_top_umi,
                 dplyr::all_of(hiv_categories))
}


# qc: named SeuratPipeline QC arguments (min_features, max_features,
# max_counts, max_mito) for this dataset; empty keeps the pipeline defaults.
load_human_sample <- function(sample, counts_root, rds_dir, qc = list()) {
  rds <- file.path(rds_dir, paste0(sample, ".rds"))
  if (file.exists(rds)) return(readRDS(rds))
  mtx <- file.path(counts_root, paste0("GRCh38_HXB2_", sample),
                   "outs", "filtered_feature_bc_matrix")
  obj <- do.call(SeuratPipeline, c(list(mtx, sample_name = sample, rds = FALSE,
                                        species = "Human", remove_mac = TRUE),
                                   qc))
  dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, rds)
  obj
}

# Markers of the cell-type gate for calling a cell HIV+
hiv_cd3_genes   <- c("CD3D", "CD3E", "CD3G")
hiv_cd4_genes   <- "CD4"
hiv_cd8_genes   <- "CD8A"

# Per-cell counts for a handful of marker genes, works on merged objects
marker_counts <- function(obj, genes) {
  rna <- obj[["RNA"]]
  out <- matrix(0, nrow = length(genes), ncol = ncol(obj),
                dimnames = list(genes, colnames(obj)))
  for (lyr in Layers(rna, search = "counts")) {
    mat <- LayerData(rna, layer = lyr)
    present <- intersect(genes, rownames(mat))
    if (length(present) == 0) next
    out[present, colnames(mat)] <- as.matrix(mat[present, , drop = FALSE])
  }
  out
}

# Protein CD4 T gates, one cell list per dataset. Sample names are disjoint
# between datasets, so the lists are stacked and matched on cellID alone.
#
# The two lists were gated by different rules because the panels differ. HIVseq
# carries no isotype control, so adtgate_final.R uses an Otsu cut per marker
# plus a depth guard, scored against the CD4 T counts that paper reports.
# Collora carries three isotype controls, so collora_adtgate.R applies the rule
# that paper defines: positive above the 90th percentile of its own isotype
# values. Collora is gated CD3+CD4+ with no CD8 term -- its cells are CD4+ T by
# negative selection, leaving protein CD8 at isotype background where no
# threshold can separate it.
#
# A dataset with no list here is left ungated, handled by `missing` below.
adt_cd4t_paths <- c(
  HIVseq  = paste0("/projects/b1042/GoyalLab/egrody/extractedData/",
                   "HIVseq/ADT/CD4T_cellIDs.csv"),
  Collora = paste0("/projects/b1042/GoyalLab/egrody/extractedData/",
                   "Collora/ADT/CD4T_cellIDs.csv"))
adt_cd4t_env <- new.env(parent = emptyenv())

adt_cd4t_table <- function(paths = adt_cd4t_paths) {
  key <- paste(paths, collapse = "|")
  if (is.null(adt_cd4t_env[[key]])) {
    found <- paths[file.exists(paths)]
    if (length(found) == 0) {
      stop("No ADT CD4 T list found: ", paste(paths, collapse = ", "),
           " -- run adtgate_final.R / collora_adtgate.R,",
           " or set hiv_tcell_gate to FALSE.")
    }
    if (length(found) < length(paths)) {
      message("ADT gate: no cell list for ",
              paste(names(paths)[!file.exists(paths)], collapse = ", "),
              " -- those samples fall back to the `missing` gate")
    }
    # Only the three columns the gate needs, since the two lists carry
    # different per-dataset diagnostics alongside them
    adt_cd4t_env[[key]] <- bind_rows(lapply(found, function(f)
      read.csv(f, stringsAsFactors = FALSE) %>%
        select(sample, cellID, adt_ok)))
  }
  adt_cd4t_env[[key]]
}

# Per-cell CD4 T call from the ADT gate, with a per-library fallback for the
# libraries whose ADT is too shallow to gate (PID145, three of four PID8026).
# `missing` is what happens to those libraries: "gene" runs the `cd4t_pass`
# gene gate on them, "pass" keeps every cell, "drop" removes every cell.
# Returns the call and the gate that made it, one entry per cell, so a table
# built on a mix of gates can say which library used which.
#
# Protein reaches the gate two ways. HIVseq and Collora were gated by their own
# scripts into the `adt_cd4t_paths` cell lists, keyed by library. A dataset that
# instead carries its protein panel inside the Seurat object supplies a per-cell
# `cell_type_adt` call made when the object was built (Giron, whose CITE-seq
# panel arrives in the same cellranger matrix as the transcriptome); that call
# is used directly, so those cells never fall through to `missing`.
adt_cd4t_pass <- function(obj, paths = adt_cd4t_paths,
                          missing = c("gene", "pass", "drop")) {
  missing <- match.arg(missing)
  adt_cd4t <- adt_cd4t_table(paths) %>% filter(adt_ok)
  gated_samples <- unique(adt_cd4t$sample)
  cell_sample <- as.character(obj$orig.ident)
  in_obj <- if ("cell_type_adt" %in% colnames(obj@meta.data)) {
    !is.na(obj$cell_type_adt)
  } else {
    rep(FALSE, ncol(obj))
  }
  ungated <- !in_obj & !cell_sample %in% gated_samples
  pass <- colnames(obj) %in% adt_cd4t$cellID
  # Guarded: an object with no `cell_type_adt` column must not touch it, so the
  # cell-list datasets (HIVseq, Collora) reach the `missing` branch below.
  if (any(in_obj)) {
    pass[in_obj] <- as.character(obj$cell_type_adt)[in_obj] == "CD4_T"
  }
  gate <- rep("ADT", ncol(obj))
  if (any(ungated)) {
    message("ADT gate: no usable ADT for ",
            paste(unique(cell_sample[ungated]), collapse = ", "),
            " -- those cells ",
            switch(missing, gene = "gated on CD3 + lineage-exclusion genes",
                            pass = "kept", drop = "dropped"))
    pass[ungated] <- switch(missing,
                            gene = hiv_tcell_pass(obj)[ungated],
                            pass = TRUE,
                            drop = FALSE)
    gate[ungated] <- switch(missing, gene = "gene", pass = "none",
                            drop = "none")
  }
  list(pass = pass, gate = gate)
}

# --- CD4 T cell gate ---------------------------------------------------------
# Negative selection: a CD3+ cell is a CD4 T cell unless it carries another
# lineage's markers. 
#
# Exclusion requires `cd4t_exclusion_min` counts to be robust to ambient molecules
#
# NK markers minimally overlap with CD4 Ts. Left out NKG7/GNLY/KLRD1 which are 
# also expressed by cytotoxic CD4 T
cd4t_exclusion_markers <- list(
  CD8      = c("CD8A", "CD8B"),
  B        = c("MS4A1", "CD79A", "CD79B", "IGHM"),
  myeloid  = c("CD14", "LYZ", "VCAN", "S100A8", "S100A9"),
  DC       = c("CD1C", "FCER1A", "LILRA4", "CLEC4C"),
  NK       = c("KLRF1", "NCR1"),
  platelet = c("PPBP", "PF4"))
cd4t_exclusion_min <- 2

# Every gene the gate reads. Marker tables must carry all of these.
cd4t_gate_genes <- unique(c(hiv_cd3_genes,
                            unlist(cd4t_exclusion_markers, use.names = FALSE)))

# CD3-undetected cell is a detection failure, an excluded cell is another
# lineage
cd4t_call <- function(x, min_count = cd4t_exclusion_min) {
  counts <- if (inherits(x, "Seurat")) {
    marker_counts(x, cd4t_gate_genes)
  } else {
    t(as.matrix(x[, intersect(cd4t_gate_genes, colnames(x)), drop = FALSE]))
  }
  # A marker absent from the reference cannot exclude
  absent <- setdiff(unlist(cd4t_exclusion_markers, use.names = FALSE), rownames(counts))
  if (length(absent) > length(unlist(cd4t_exclusion_markers)) / 2) {
    warning("cd4t_call: ", length(absent), " of ",
            length(unlist(cd4t_exclusion_markers, use.names = FALSE)),
            " exclusion markers missing -- rebuild the marker table with",
            " `cd4t_gate_genes`, or this is a CD3-only gate", call. = FALSE)
  }
  hit <- function(genes, n) {
    genes <- intersect(genes, rownames(counts))
    if (length(genes) == 0) return(rep(FALSE, ncol(counts)))
    colSums(counts[genes, , drop = FALSE]) >= n
  }
  cd3      <- hit(hiv_cd3_genes, 1)
  excluded <- Reduce(`|`, lapply(cd4t_exclusion_markers, hit, n = min_count))
  list(cd3 = cd3, excluded = excluded, pass = cd3 & !excluded)
}

# TRUE for cells passing the CD4 T gate.
cd4t_pass <- function(x, min_count = cd4t_exclusion_min) cd4t_call(x, min_count)$pass

hiv_tcell_pass <- function(obj) cd4t_pass(obj)

# Per-cell marker tables written by stepFive_markers.R
cd4t_expression_dir <- paste0(data.root,
                              "Figure1/singleCell/Seurat/expressionDF/")

cd4t_gate_table <- function(sample_names, expression.dir = cd4t_expression_dir,
                            prefix = FALSE) {
  one <- function(sample_name) {
    f <- file.path(expression.dir, paste0(sample_name, "_CD4T_counts.csv"))
    if (!file.exists(f)) {
      warning("cd4t_gate_table: no marker table at ", f,
              " -- run stepFive_markers.R", call. = FALSE)
      return(NULL)
    }
    df <- read.csv(f, row.names = "X")
    call <- cd4t_call(df)
    cellID <- as.character(df$cellID)
    tibble::tibble(cellID = if (prefix) paste0(sample_name, "_", cellID) else cellID,
                   cd4t = call$pass, excluded = call$excluded)
  }
  out <- bind_rows(lapply(sample_names, one))
  if (nrow(out) == 0) NULL else out
}

# zero `viral_cols` outside the gate
apply_cd4t_gate <- function(df, sample_names, viral_cols,
                            expression.dir = cd4t_expression_dir,
                            prefix = FALSE) {
  gate <- cd4t_gate_table(sample_names, expression.dir, prefix)
  if (is.null(gate)) return(df)
  df %>%
    left_join(gate, by = "cellID") %>%
    filter(!excluded %in% TRUE) %>%
    mutate(across(any_of(viral_cols), ~ ifelse(cd4t %in% TRUE, .x, 0))) %>%
    select(-cd4t, -excluded)
}

# The HIV+ call, in the terms each source paper states it in. Both rules read
# the same bamsort tables; they differ in whether a cell qualifies on reads or
# on molecules.
#   "reads"   >= min_hiv_reads HXB2 reads per barcode. This is HIV-seq's own
#             published rule (Frouard et al.), and the rule our SIV and Giron
#             TenX arms use, so it is the default.
#   "collora" the rule of Collora et al. (Immunity 2022,
#             doi 10.1016/j.immuni.2022.03.004) STAR Methods: at least 2 HIV-1
#             UMIs, or 4 reads of a single UMI. Strictly more stringent than
#             "reads" at min_hiv_reads = 2, since 2 UMIs carry >= 2 reads and 4
#             reads of one UMI are 4 reads.
# The read clause needs `primary_count` in the bamsort tables: without it a read
# spanning two adjacent junction windows is charged to both, so a 2-read
# molecule reads as 4 and clears the clause on its own. hxb2_per_cell prefers
# that column and warns when a table predates it.
#   "frouard" the rule of Frouard et al. (Nat Commun 2026,
#             doi 10.1038/s41467-026-68797-3): one HIV transcript, i.e. a single
#             HXB2 UMI, which they set from HIV-uninfected control donors in whom
#             "not a single HIV transcript was identified". Looser than "reads"
#             at min_hiv_reads = 2, since a lone UMI of one read clears it.
hiv_collora_min_umis      <- 2   # HIV-1 UMIs per barcode
hiv_collora_min_umi_reads <- 4   # reads of a single UMI
hiv_frouard_min_umis      <- 1   # HIV transcripts (UMIs) per barcode

hiv_pos_call <- function(hxb2, hxb2_umis, hxb2_top_umi, rule = "reads") {
  switch(rule,
         reads   = hxb2 >= min_hiv_reads,
         collora = hxb2_umis    >= hiv_collora_min_umis |
                   hxb2_top_umi >= hiv_collora_min_umi_reads,
         frouard = hxb2_umis    >= hiv_frouard_min_umis,
         stop("unknown HIV+ rule: ", rule, call. = FALSE))
}

# One-line statement of a rule, for the prints and the SI
hiv_pos_rule_text <- function(rule = "reads") {
  switch(rule,
         reads   = paste0(">=", min_hiv_reads, " HXB2 reads per barcode"),
         collora = paste0(">=", hiv_collora_min_umis, " HIV-1 UMIs, or >=",
                          hiv_collora_min_umi_reads, " reads of a single UMI"),
         frouard = paste0(">=", hiv_frouard_min_umis,
                          " HIV transcript (UMI) per barcode"),
         stop("unknown HIV+ rule: ", rule, call. = FALSE))
}

# The rule a public dataset is called under: its own paper's, from
# `public_paths$<dataset>$hiv_rule`, defaulting to "reads".
# `[[` and not `$`: `$` partial-matches on a list, so `$hiv_rule` would return
# the value of `hiv_rule_published` for a dataset that carries only the latter,
# silently changing the primary call.
public_hiv_rule <- function(dataset) {
  rule <- public_paths[[dataset]][["hiv_rule"]]
  if (is.null(rule)) "reads" else rule
}

# The rule that dataset's own paper published, for the secondary comparison in
# SuppFigs. It equals `hiv_rule` wherever we already call a dataset under its
# own rule, and differs only where we deliberately use a stricter one.
public_hiv_rule_published <- function(dataset) {
  rule <- public_paths[[dataset]][["hiv_rule_published"]]
  if (is.null(rule)) public_hiv_rule(dataset) else rule
}

# Attach per-cell HXB2 reads, molecules and isoform categories (bamsort) and set
# HIV_status under `rule` (hiv_pos_call), in a gated CD3+CD4+ cell
attach_hxb2 <- function(obj, samples, bamsort_dir, rule = "reads") {
  hxb2 <- purrr::map_dfr(samples, hxb2_per_cell, bamsort_dir = bamsort_dir)
  m <- match(colnames(obj), hxb2$cellID)
  pick <- function(col) ifelse(is.na(m), 0, hxb2[[col]][m])
  obj$hxb2         <- pick("hxb2")
  obj$hxb2_umis    <- pick("hxb2_umis")
  obj$hxb2_top_umi <- pick("hxb2_top_umi")
  obj$HXB2_US    <- pick("US")
  obj$HXB2_US.SS <- pick("US.SS")
  obj$HXB2_any   <- pick("any")
  obj$HXB2_SS.MS <- pick("SS.MS")
  obj$HXB2_SS    <- pick("SS")
  obj$HXB2_MS    <- pick("MS")

  
  # CD4 T called by ADT protein when available, else by `cd4t_pass` gene gate
  cd4t <- if (hiv_tcell_gate) adt_cd4t_pass(obj, missing = "gene")
          else                list(pass = TRUE, gate = "none")
  obj$CD4_Tcell      <- cd4t$pass
  obj$CD4_Tcell_gate <- cd4t$gate
  hxb2_cols <- c("hxb2", "hxb2_umis", "hxb2_top_umi",
                 paste0("HXB2_", hiv_categories))
  for (cc in hxb2_cols) {
    obj@meta.data[[cc]] <- ifelse(obj$CD4_Tcell, obj@meta.data[[cc]], 0)
  }
  obj$HIV_rule   <- rule
  obj$HIV_status <- ifelse(hiv_pos_call(obj$hxb2, obj$hxb2_umis,
                                        obj$hxb2_top_umi, rule),
                           "HIV+", "HIV-")
  obj
}

# Integrate participants separately
public_group_id <- function(condition, participant) {
  paste0(condition, "_", participant)
}
public_group_parts <- function(group) {
  list(condition   = sub("_[^_]*$", "", group),
       participant = sub(".*_", "", group))
}

# HIV+ out of the gated CD4 T cells, splice+ out of the HIV+ cells
public_add_pcts <- function(df) {
  df %>%
    mutate(pct_HIV_pos    = round(100 * HIV_pos / n_CD4_Tcells, 3),
           pct_splice_pos = ifelse(HIV_pos > 0,
                                   round(100 * splice_pos / HIV_pos, 1),
                                   NA_real_))
}

# Summary of HIV+ and splice+ cells counts and percentages per sample
public_sample_table <- function(obj, group) {
  parts <- public_group_parts(group)
  tibble(condition   = parts$condition,
         participant = parts$participant,
         sample      = as.character(obj$orig.ident),
         tcell       = obj$CD4_Tcell,
         gate        = obj$CD4_Tcell_gate,
         hiv         = obj$HIV_status == "HIV+",
         # Same splice+ definition as the SIV arm: SS or MS or SS-MS.
         spliced     = obj$HXB2_SS.MS + obj$HXB2_SS + obj$HXB2_MS) %>%
    group_by(condition, participant, sample) %>%
    summarize(gate         = paste(sort(unique(gate)), collapse = "+"),
              n_cells      = n(),
              n_CD4_Tcells = sum(tcell),
              HIV_pos      = sum(hiv),
              splice_pos   = sum(hiv & spliced >  0),
              splice_neg   = sum(hiv & spliced == 0),
              .groups = "drop") %>%
    public_add_pcts()
}

# Totals over the per-sample table, summed across the grouping columns given
public_totals <- function(by_sample, ...) {
  by_sample %>%
    group_by(...) %>%
    summarize(n_samples = n(),
              gate      = paste(sort(unique(unlist(strsplit(gate, "\\+")))),
                                collapse = "+"),
              across(c(n_cells, n_CD4_Tcells, HIV_pos, splice_pos, splice_neg),
                     sum),
              .groups = "drop") %>%
    public_add_pcts()
}

# Integrate all samples of one condition then attach the HXB2 isoform join.
integrate_condition <- function(samples, counts_root, bamsort_dir,
                                sample_rds_dir, integrated_rds,
                                loader = load_human_sample, qc = list(),
                                rule = "reads") {
  if (file.exists(integrated_rds)) {
    message("Public data: loading cached integration -> ", basename(integrated_rds))
    obj <- readRDS(integrated_rds)
  } else {
    message("Public data: integrating ", length(samples), " sample(s) -> ",
            basename(integrated_rds), " (compute-heavy)")
    individual <- lapply(samples, function(s) {
      o <- loader(s, counts_root, sample_rds_dir, qc = qc)
      o$orig.ident <- s
      RenameCells(o, add.cell.id = s)
    })
    names(individual) <- samples
    obj <- Reduce(function(x, y) merge(x, y), individual)
    obj <- SCTransform(obj, vst.flavor = "v2", verbose = FALSE)
    obj <- RunPCA(obj, verbose = FALSE)
    if (length(samples) > 1) {
      obj <- IntegrateLayers(obj, method = HarmonyIntegration,
                             orig.reduction = "pca", new.reduction = "harmony",
                             normalization.method = "SCT", verbose = FALSE)
      reduction <- "harmony"
    } else {
      reduction <- "pca"
    }
    obj <- FindNeighbors(obj, reduction = reduction, dims = 1:REDUCTION_DIMENSIONS, verbose = FALSE)
    obj <- FindClusters(obj, resolution = 0.7, algorithm = 4, random.seed = 1,
                        verbose = FALSE)
    obj <- RunUMAP(obj, dims = 1:REDUCTION_DIMENSIONS, reduction = reduction, verbose = FALSE)
    dir.create(dirname(integrated_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(obj, integrated_rds)
  }
  attach_hxb2(obj, samples, bamsort_dir, rule = rule)
}

# The results object every public-data chunk reads: the per-group Seurat objects
# (one participant each), the per-sample table and its two rollups
public_results <- function(objs) {
  by_sample <- bind_rows(purrr::imap(objs, public_sample_table))
  list(objects        = objs,
       by_sample      = by_sample,
       by_participant = public_totals(by_sample, condition, participant),
       summary        = public_totals(by_sample, condition))
}

# Build and cache the per-sample objects of a dataset
public_build_rds <- function(dataset, samples, force = FALSE) {
  paths <- public_paths[[dataset]]
  qc <- public_qc(dataset)
  purrr::map_dfr(samples, function(sample) {
    rds <- file.path(paths$rds, paste0(sample, ".rds"))
    cached <- file.exists(rds)
    if (cached && force) file.remove(rds)
    message("Public data build: ", sample,
            if (cached && !force) " (cached)" else " (building)")
    obj <- load_human_sample(sample, paths$counts, paths$rds, qc = qc)
    out <- tibble(sample = sample, cells = ncol(obj), features = nrow(obj),
                  status = if (cached && !force) "cached" else "built")
    rm(obj); gc()
    out
  })
}

# The dataset's SeuratPipeline QC arguments
public_qc <- function(dataset) public_paths[[dataset]][["qc"]] %||% list()

# Integrate each participant of a dataset (or load its cached integration)
public_extract <- function(dataset, conditions) {
  paths <- public_paths[[dataset]]
  map <- public_maps[[dataset]] %>%
    filter(condition %in% conditions) %>%
    mutate(group = public_group_id(condition, participant))
  if (nrow(map) == 0) {
    warning("No ", dataset, " samples map to this pass; nothing to integrate.",
            call. = FALSE)
    return(NULL)
  }
  dir.create(paths$integrated, recursive = TRUE, showWarnings = FALSE)
  groups <- unique(map$group)
  public_results(setNames(lapply(groups, function(group) {
    integrate_condition(map$sample[map$group == group], paths$counts,
                        paths$bamsort, paths$rds,
                        file.path(paths$integrated, paste0(group, ".rds")),
                        qc = public_qc(dataset),
                        rule = public_hiv_rule(dataset))
  }), groups))
}

# Reload the cached per-participant RDS WITHOUT re-integrating, re-applying the
# HXB2 isoform join plus the CD4 T-cell gate
public_reload <- function(dataset, conditions) {
  paths <- public_paths[[dataset]]
  groups <- sub("\\.rds$", "", basename(
    Sys.glob(file.path(paths$integrated, "*.rds"))))

  # Caches written before the per-participant split are named for the condition
  # alone and pool every participant, so they are never loaded
  pooled <- groups[!grepl("_", groups)]
  if (length(pooled) > 0) {
    message("Public data: ignoring pre-participant cache in ", paths$integrated,
            ": ", paste(paste0(pooled, ".rds"), collapse = ", "))
  }
  groups <- groups[grepl("_", groups)]
  groups <- groups[vapply(groups, function(group)
    public_group_parts(group)$condition, character(1)) %in% conditions]

  if (length(groups) == 0) {
    warning("No cached per-participant ", dataset, " integration for ",
            paste(conditions, collapse = ", "), " in ", paths$integrated,
            " -- run the Integrate chunk.", call. = FALSE)
    return(NULL)
  }
  public_results(setNames(lapply(groups, function(group) {
    message("Public data reload: ", group)
    obj <- readRDS(file.path(paths$integrated, paste0(group, ".rds")))
    attach_hxb2(obj, unique(as.character(obj$orig.ident)), paths$bamsort,
                rule = public_hiv_rule(dataset))
  }), groups))
}

# Load every dataset
public_load_all <- function(datasets, conditions, loader = public_reload) {
  for (dataset in datasets) {
    res <- loader(dataset, conditions)
    if (is.null(res)) next
    assign(public_datasets[[dataset]]$results, res, envir = globalenv())
    cat("Public data:", dataset, "-", nrow(res$by_participant), "participants,",
        nrow(res$by_sample), "libraries,", sum(res$by_sample$n_cells), "cells,",
        sum(res$by_sample$HIV_pos), "HIV+ -",
        paste(names(res$objects), collapse = ", "), "\n")
  }
  invisible(NULL)
}

# Condition totals, then the participants they pool, then the libraries of each
# participant with the share of their participant's total that each contributes
public_report <- function(res, dataset, conditions) {
  ordered <- function(df) df %>%
    mutate(condition = factor(condition, levels = conditions)) %>%
    arrange(condition, desc(HIV_pos))

  cat("\n=== Public data:", dataset, "- joint GRCh38+HXB2, HIV+ on",
      hiv_pos_rule_text(public_hiv_rule(dataset)), "===\n")
  print(as.data.frame(ordered(res$summary)))

  # The condition totals average over participants whose HIV+ rates differ by
  # more than an order of magnitude, so these rates are the ones to compare.
  # Each row here is one integrated object.
  cat("\n--- per participant ---\n")
  print(as.data.frame(ordered(res$by_participant)))

  cat("\n--- per library, with its share of the participant ---\n")
  print(as.data.frame(ordered(res$by_sample) %>%
    group_by(condition, participant) %>%
    mutate(across(c(HIV_pos, splice_pos),
                  ~ if (sum(.x) > 0) round(100 * .x / sum(.x), 1)
                    else rep(NA_real_, length(.x)),
                  .names = "pct_of_{.col}")) %>%
    ungroup() %>%
    select(condition, participant, sample, gate, n_cells, n_CD4_Tcells,
           HIV_pos, pct_HIV_pos, pct_of_HIV_pos,
           splice_pos, pct_splice_pos, pct_of_splice_pos)))
  invisible(res)
}


# Minimum cells per group for a FindMarkers contrast
public_min_cells <- 3
# Genes dropped before DEG
public_exclude_genes <- c("CD4", "CD8A", "CD8B")

# splice+/splice- among the HIV+ cells; every other cell gets `other`
public_splice_call <- function(obj, other = "other") {
  spliced <- obj$HXB2_SS.MS + obj$HXB2_MS
  dplyr::case_when(obj$HIV_status == "HIV+" & spliced >  0 ~ "splice+",
                   obj$HIV_status == "HIV+" & spliced == 0 ~ "splice-",
                   TRUE                                    ~ other)
}

# Per-cell table of UMAP coordinates, HXB2 counts and the two identity calls.
public_cell_table <- function(obj) {
  coords <- Embeddings(obj, reduction = "umap")
  tibble(cellID     = colnames(obj),
         UMAP1      = coords[, 1],
         UMAP2      = coords[, 2],
         hxb2       = obj$hxb2,
         spliced    = obj$HXB2_SS.MS + obj$HXB2_MS,
         HIV_status = obj$HIV_status,
         infect_identity = ifelse(HIV_status == "HIV+",
                                  "infected", "uninfected"),
         splice_identity = public_splice_call(obj))
}

public_contrasts <- list(
  infected = list(
    suffix = "_DEG_infected.csv",
    ident.1 = "infected", ident.2 = "uninfected",
    colors = c(infected = "red", uninfected = "gray50"),
    groups = function(obj) ifelse(obj$HIV_status == "HIV+",
                                  "infected", "uninfected")),
  splice = list(
    suffix = "_DEG_splice.csv",
    ident.1 = "splice+", ident.2 = "splice-",
    colors = c(`splice+` = "purple", `splice-` = "orange"),
    groups = function(obj) public_splice_call(obj, other = NA_character_))
)

# Accessors for the global results objects the Load chunk assigns
public_res <- function(dataset) {
  get(public_datasets[[dataset]]$results, envir = globalenv())
}
public_loaded <- function() {
  loaded <- names(public_datasets)[vapply(public_datasets,
    function(d) exists(d$results, envir = globalenv()), logical(1))]
  if (length(loaded) == 0) {
    warning("No public-data results loaded; run the public_load chunk first.")
  }
  loaded
}
# Every integrated object of a dataset, in condition order then participant.
# The handler is passed the object and its "<condition>_<participant>" label.
public_walk_groups <- function(dataset, handler) {
  results <- public_res(dataset)
  groups <- names(results$objects)
  conditions <- vapply(groups, function(group)
    public_group_parts(group)$condition, character(1))
  groups <- groups[order(match(conditions,
                               public_datasets[[dataset]]$conditions), groups)]
  for (group in groups) handler(results$objects[[group]], group)
  invisible(NULL)
}

# The per-sample tables of every loaded dataset
public_all_by_sample <- function() {
  purrr::map_dfr(public_loaded(), function(dataset) {
    public_res(dataset)$by_sample %>% mutate(dataset = dataset, .before = 1)
  })
}

# Sequencing depth per library, read from the metrics_summary.csv
public_cellranger_metrics <- function(dataset, samples) {
  counts_root <- public_paths[[dataset]][["counts"]]
  purrr::map_dfr(samples, function(sample) {
    f <- file.path(counts_root, paste0("GRCh38_HXB2_", sample), "outs",
                   "metrics_summary.csv")
    if (!file.exists(f)) {
      warning("No metrics_summary.csv for ", dataset, " ", sample)
      return(NULL)
    }
    m <- read.csv(f, check.names = FALSE, colClasses = "character")
    num <- function(col) as.numeric(gsub("[,%]", "", m[[col]]))
    tibble(sample               = sample,
           cells                = num("Estimated Number of Cells"),
           reads                = num("Number of Reads"),
           reads_per_cell       = num("Mean Reads per Cell"),
           median_umis_per_cell = num("Median UMI Counts per Cell"))
  })
}

# Categorical colors for the two public datasets
public_dataset_colors <- c(Collora = "#2a78d6", HIVseq = "#eb6834")

# --- CD4 T subsets, public HIV signature, set-level tests ----------------------
# CD4 T subset marker sets. A cell's subset call is the highest of their
# AddModuleScore scores.
cd4t_subset_markers <- list(
  naive_Tcm = c("CCR7", "SELL", "TCF7", "LEF1", "IL7R", "BACH2", "NOSIP", "NR3C2"),
  Tem       = c("S100A4", "S100A6", "GZMK", "KLRB1", "ANXA1", "CXCR3", "AHR"),
  cytotoxic = c("GZMB", "GZMA", "NKG7", "PRF1", "KLRD1", "CTSW", "CST7", "CCL5", "HOPX"),
  Treg      = c("FOXP3", "IKZF2", "IL2RA", "CTLA4", "TNFRSF18"),
  cycling   = c("MKI67", "TOP2A", "CDK1", "BIRC5", "STMN1", "NUSAP1"))

# Genes called up (or down) in HIV+ vs HIV- CD4 T cells: p_val_adj < pval and
# |avg_log2FC| > fc in at least min_tables of the given per-participant
# FindMarkers tables, and never significant the other way.
public_hiv_signature <- function(files, fc, pval, min_tables) {
  calls <- purrr::map_dfr(files, function(f) {
    read.csv(f) %>%
      transmute(gene,
                direction = case_when(p_val_adj < pval & avg_log2FC >  fc ~ "up",
                                      p_val_adj < pval & avg_log2FC < -fc ~ "down"))
  }) %>%
    filter(!is.na(direction)) %>%
    count(gene, direction) %>%
    tidyr::pivot_wider(names_from = direction, values_from = n, values_fill = 0)
  for (direction in c("up", "down")) {
    if (!direction %in% colnames(calls)) calls[[direction]] <- 0
  }
  list(PUBLIC_HIV_UP   = calls$gene[calls$up   >= min_tables & calls$down == 0],
       PUBLIC_HIV_DOWN = calls$gene[calls$down >= min_tables & calls$up   == 0])
}

# Human MSigDB Hallmark gene symbols, named without the HALLMARK_ prefix. Reads
# msigdbr; if the package will not load (its curl dependency needs
# libssl.so.10, absent from these nodes), reads msigdbr's own local cache of the
# same MSigDB release.
hallmark_gene_sets <- function(cache = "~/.cache/R/msigdbr/msigdb.2026.1.Hs.H.rds") {
  hallmark <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H") %>%
      transmute(gene_symbol, gs_name),
    error = function(e) {
      message("msigdbr unavailable (", conditionMessage(e), ") -- reading ", cache)
      readRDS(cache) %>% transmute(gene_symbol = db_gene_symbol, gs_name)
    })
  split(hallmark$gene_symbol, sub("^HALLMARK_", "", hallmark$gs_name))
}

# Percentile of each value among the values of its group: (rank - 0.5) / n,
# ties averaged. Uniform on (0, 1) within every group.
within_group_percentile <- function(x, group) {
  ave(x, group, FUN = function(v) (rank(v, ties.method = "average") - 0.5) / length(v))
}

# Null means of `values` for sets drawn like the positives: each draw takes, per
# stratum, as many non-positive cells as there are positives in that stratum.
matched_mean_null <- function(values, strata, positive, nperm) {
  wanted <- table(strata[positive])
  pools <- split(values[!positive], strata[!positive])[names(wanted)]
  replicate(nperm, mean(unlist(mapply(function(pool, n) pool[sample.int(length(pool), n)],
                                      pools, wanted, SIMPLIFY = FALSE))))
}

# Bar of the per-dataset mean with one dot per library on top
public_bar_dots <- function(df, title, group_col = "dataset",
                            colors = public_dataset_colors, saveas = NULL,
                            width = NULL, height = NULL) {
  df <- df %>% filter(!is.na(value)) %>% mutate(.group = .data[[group_col]])
  if (is.null(width))  width  <- 3.4 + 2.0 * n_distinct(df$condition)
  if (is.null(height)) height <- 1.6 + 1.4 * n_distinct(df$metric)

  # Unpaired Wilcoxon per facet, skipped where only one group is present
  tests <- df %>%
    group_by(condition, metric) %>%
    group_modify(~ {
      sets <- sort(unique(.x$.group))
      vals <- split(.x$value, factor(.x$.group, levels = sets))
      paired <- length(sets) == 2
      tibble(group1 = sets[1],
             group2 = if (paired) sets[2] else NA_character_,
             n1 = length(vals[[1]]),
             n2 = if (paired) length(vals[[2]]) else 0L,
             p  = if (paired) suppressWarnings(
                    wilcox.test(vals[[1]], vals[[2]])$p.value) else NA_real_,
             y  = max(.x$value) * 1.13)
    }) %>%
    ungroup() %>%
    mutate(label = ifelse(is.na(p), "", sprintf("p = %.3g", p)))

  plot <- ggplot(df, aes(x = .group, y = value, fill = .group)) +
    stat_summary(fun = mean, geom = "col", width = 0.6, alpha = 0.85) +
    geom_point(aes(color = .group), size = 2.4, shape = 21, fill = "white",
               stroke = 0.9,
               position = position_jitter(width = 0.13, height = 0, seed = 1)) +
    geom_text(data = tests, aes(x = 1.5, y = y, label = label),
              inherit.aes = FALSE, size = 3, color = "#3d3d38") +
    facet_grid(metric ~ condition, scales = "free_y", switch = "y") +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0, 0.22))) +
    labs(title = title,
         subtitle = "Bars are group means; each dot is one library",
         x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.subtitle = element_text(color = "#5c5c57", size = 9),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0))

  if (!is.null(saveas)) {
    dir.create(dirname(saveas), recursive = TRUE, showWarnings = FALSE)
    ggsave(saveas, plot, width = width, height = height)
  }
  print(as.data.frame(tests %>%
    select(condition, metric, group1, n1, group2, n2, p)))
  plot
}

# Colors for the HIV-seq library chemistries
public_library_colors <- c(NOTHING = "#2a78d6", Primer = "#eb6834")

# Per-library HIV+/splice+ counts WITHOUT integration. `libraries` needs a
# sample, condition and participant column.
public_library_table <- function(dataset, libraries) {
  paths <- public_paths[[dataset]]
  purrr::pmap_dfr(list(libraries$sample, libraries$condition,
                       libraries$participant),
                  function(sample, condition, participant) {
    obj <- load_human_sample(sample, paths$counts, paths$rds,
                             qc = public_qc(dataset))
    obj$orig.ident <- sample
    obj <- RenameCells(obj, add.cell.id = sample)
    obj <- attach_hxb2(obj, sample, paths$bamsort,
                       rule = public_hiv_rule(dataset))
    out <- public_sample_table(obj, public_group_id(condition, participant))
    rm(obj); gc()
    out
  })
}

## TRUE when every input is readable, naming anything that is not. Every input
## now lives under its figure's own tree, so a missing file means the step that
## writes it has not been run -- there is no second location to check.
check_migrated <- function(...) {
  paths <- c(...)
  missing <- paths[!file.exists(paths)]
  for (path in missing) cat("Data missing:", path, "\n")
  length(missing) == 0
}

## read.csv/readRDS wrapper that returns NULL with the message above instead of
## erroring when the file is absent.
read_migrated <- function(path, reader = read.csv, ...) {
  if (!check_migrated(path)) return(NULL)
  reader(path, ...)
}

# The Figure1 SIV samples were counted in three separate cellranger campaigns, so
# their outs/ directories sit under three different experiment trees. `runs` is
# the single source of truth for where each sample's cellranger output came from;
# collect_cellranger_outs() flattens it into one directory so downstream code
# only needs a sample name.
figure1_cellranger_runs <- data.frame(
  samples = c("D0", "D195", "D83", "Invitro", "Pacute", "W0", "W2"),
  runs = c(
    "/projects/b1042/GoyalLab/egrody/extractedData/Figure1/singleCell/counts/Mmul_10_mac239/Mmul_10_mac239_LP29_D0",
    "/projects/b1042/GoyalLab/egrody/extractedData/Figure1/singleCell/counts/Mmul_10_mac239/Mmul_10_mac239_LP29_D195",
    "/projects/b1042/GoyalLab/egrody/extractedData/Figure1/singleCell/counts/Mmul_10_mac239/Mmul_10_mac239_LP29_D83",
    "/projects/b1042/GoyalLab/egrody/extractedData/Figure1/singleCell/counts/Mmul_10_mac239/Mmul_10_mac239_invitro",
    "/projects/b1042/GoyalLab/egrody/extractedData/Figure1/singleCell/counts/Mmul_10_mac239/Mmul_10_mac239_P_acute_combined",
    "/projects/b1042/GoyalLab/egrody/extractedData/Figure1/singleCell/counts/Mmul_10_mac239/Mmul_10_mac239_W0",
    "/projects/b1042/GoyalLab/egrody/extractedData/Figure1/singleCell/counts/Mmul_10_mac239/Mmul_10_mac239_W2"
  ),
  stringsAsFactors = FALSE
)

## Copy each sample's cellranger outs/ into one flat counts directory, so that
## every reader can build a path as <counts.dir>/outs_<sample>/<matrix>/.
## `runs` needs a samples and a runs column, where runs is the cellranger run
## directory that contains outs/.
## `contents` limits what gets copied; the default takes only the two matrix
## directories, which is what Read10X and SeuratPipeline need. The BAMs are
## multiple GB per sample and the bamsort pipeline reads them from the original
## runs, so they are left in place. Pass contents = NULL to copy all of outs/.
collect_cellranger_outs <- function(runs = figure1_cellranger_runs,
                                    counts.dir = "/projects/b1042/GoyalLab/egrody/extractedData/Figure1/singleCell/counts/",
                                    contents = c("raw_feature_bc_matrix",
                                                 "filtered_feature_bc_matrix"),
                                    overwrite = FALSE) {
  dir.create(counts.dir, recursive = TRUE, showWarnings = FALSE)

  report <- data.frame(sample = character(), source = character(),
                       destination = character(), files = integer(),
                       status = character(), stringsAsFactors = FALSE)

  for (i in 1:nrow(runs)) {
    sample_name <- runs$samples[i]
    source_outs <- file.path(runs$runs[i], "outs")
    destination <- file.path(counts.dir, paste0("outs_", sample_name))

    # Copy the named subdirectories, or everything in outs/ when contents is NULL
    to_copy <- if (is.null(contents)) {
      list.files(source_outs, full.names = TRUE)
    } else {
      file.path(source_outs, contents)
    }
    to_copy <- to_copy[file.exists(to_copy)]

    # An empty matrix directory copies "successfully" while leaving the
    # destination unreadable, so treat a source with no files as a failure
    n_files <- length(list.files(to_copy, recursive = TRUE))

    if (!dir.exists(source_outs)) {
      status <- "missing source"
    } else if (n_files == 0) {
      status <- "source is empty"
    } else if (dir.exists(destination) && !overwrite) {
      status <- "already present, skipped"
    } else {
      dir.create(destination, recursive = TRUE, showWarnings = FALSE)
      copied <- file.copy(to_copy, destination, recursive = TRUE,
                          overwrite = overwrite, copy.date = TRUE)
      status <- if (all(copied)) "copied" else "copy failed"
    }

    cat(sprintf("%-8s %-24s %s\n", sample_name, status, source_outs))
    report <- rbind(report, data.frame(sample = sample_name,
                                       source = source_outs,
                                       destination = destination,
                                       files = n_files,
                                       status = status,
                                       stringsAsFactors = FALSE))
  }

  incomplete <- !report$status %in% c("copied", "already present, skipped")
  if (any(incomplete)) {
    warning("No counts copied for: ",
            paste(report$sample[incomplete], collapse = ", "))
  }
  report
}

# =============================================================================
# 10X viral library metrics
# =============================================================================

# Viral-specific sequencing saturation from a Cell Ranger molecule_info.h5.
#
# Cell Ranger's metrics_summary.csv reports one Sequencing Saturation for the
# whole transcriptome. A claim about viral reads needs the saturation of the
# viral feature alone, which requires the per-molecule read counts in
# molecule_info.h5. Saturation is Cell Ranger's own definition restricted to
# that one feature,
#
#     saturation = 1 - (viral UMIs) / (viral reads)
#
# where viral reads sum the `count` field over molecules assigned to the viral
# feature. It is reported twice: over all barcodes, and over the cell-called
# barcodes in barcode_info/pass_filter.
#
# Also reports the estimated viral library complexity C from the
# Lander-Waterman/Poisson relation
#
#     observed UMIs = C * (1 - exp(-reads / C))
#
# solved for C by bisection. The observed fraction UMIs/C says how much of the
# viral library the run recovered, which bounds how many additional molecules
# deeper sequencing could add.
#
# molecule_info.h5 holds one record per molecule -- tens of millions for a single
# library -- so feature_idx is read in full once and `count` and `barcode_idx`
# are then read only at the viral positions, rather than pulling three vectors
# of that length into memory.
tenx_viral_saturation <- function(molecule_info_h5, viral_feature = "mac239") {
  # rhdf5 is an optional Bioconductor package in packages.R, attached here
  # rather than at startup like the rest of the genomics stack.
  load_packages("rhdf5")

  if (!file.exists(molecule_info_h5)) {
    stop("molecule_info.h5 not found: ", molecule_info_h5)
  }

  # Resolve the viral feature to a 0-based index. Matching is exact and
  # case-sensitive on features/name, then on features/id. Erroring on absent or
  # ambiguous beats silently selecting a host gene with a similar name.
  feature_names <- as.character(rhdf5::h5read(molecule_info_h5, "features/name"))
  feature_ids   <- as.character(rhdf5::h5read(molecule_info_h5, "features/id"))
  hits <- which(feature_names == viral_feature)
  if (!length(hits)) hits <- which(feature_ids == viral_feature)
  if (!length(hits)) {
    stop("no feature named '", viral_feature, "'. The joint reference labels ",
         "the virus as a single feature; check features/name in the h5.")
  }
  if (length(hits) > 1) {
    stop("'", viral_feature, "' matches ", length(hits), " features (",
         paste(hits - 1L, collapse = ", "), "); expected exactly one.")
  }
  viral_idx <- hits[1] - 1L

  feature_idx <- as.vector(rhdf5::h5read(molecule_info_h5, "feature_idx"))
  viral_pos <- which(feature_idx == viral_idx)
  rm(feature_idx)
  if (!length(viral_pos)) {
    stop("feature '", viral_feature, "' carries no molecules in ", molecule_info_h5)
  }

  viral_counts   <- as.numeric(rhdf5::h5read(molecule_info_h5, "count",
                                             index = list(viral_pos)))
  viral_barcodes <- as.vector(rhdf5::h5read(molecule_info_h5, "barcode_idx",
                                            index = list(viral_pos)))

  # rhdf5 returns pass_filter transposed relative to its HDF5 shape, so the
  # barcode index of each called cell is the first ROW, not the first column.
  pass_filter <- rhdf5::h5read(molecule_info_h5, "barcode_info/pass_filter")
  cells <- unique(as.vector(pass_filter[1, ]))
  in_cell <- viral_barcodes %in% cells

  rbind(
    .tenx_saturation_row("all_barcodes", length(viral_pos), sum(viral_counts)),
    .tenx_saturation_row("cell_called_barcodes", sum(in_cell),
                         sum(viral_counts[in_cell]))
  )
}

# Estimate library complexity C by solving n_umis = C * (1 - exp(-n_reads / C)).
# Returns NA when every read was a distinct UMI: with no duplicates C is
# unbounded and the bisection has no bracket.
.tenx_estimate_complexity <- function(n_reads, n_umis) {
  if (n_umis >= n_reads) return(NA_real_)
  lo <- as.numeric(n_umis)
  hi <- as.numeric(n_umis) * 1e5
  for (i in seq_len(300)) {
    mid <- (lo + hi) / 2
    if (mid * (1 - exp(-n_reads / mid)) < n_umis) lo <- mid else hi <- mid
  }
  (lo + hi) / 2
}

# One output row. Rounding matches what the standalone script wrote, so a table
# produced here is comparable to one produced before the port.
.tenx_saturation_row <- function(label, n_umis, n_reads) {
  saturation <- if (n_reads) 1 - n_umis / n_reads else NA_real_
  complexity <- .tenx_estimate_complexity(n_reads, n_umis)
  data.frame(barcode_set          = label,
             viral_umis           = n_umis,
             viral_reads          = n_reads,
             reads_per_umi        = round(n_reads / max(1, n_umis), 4),
             viral_saturation     = round(saturation, 6),
             estimated_complexity = round(complexity, 1),
             fraction_observed    = round(n_umis / complexity, 6),
             stringsAsFactors     = FALSE)
}
