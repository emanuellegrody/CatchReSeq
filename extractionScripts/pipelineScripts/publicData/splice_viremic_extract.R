# Per-cell table of the infected CD4 T cells of the viremic public datasets
# (Collora, HIV-seq Primer) and the untreated SIV samples (Invitro, Pacute, W2),
# with splice status, viral burden, host library size and module scores. Run on
# a compute node; writes splice_viremic_cells.rds under extractedData.
chunk_code <- function(label) {
  start <- grep(paste0("^```\\{r ", label, "[,}]"), rmd)
  stopifnot(length(start) == 1)
  end <- start + which(grepl("^```\\s*$", rmd[(start + 1):length(rmd)]))[1]
  rmd[(start + 1):(end - 1)]
}
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

# Chunk bodies are taken from the notebook verbatim so this cannot drift from it.
rmd <- readLines(file.path(plot_dir, "Figure3.Rmd"), warn = FALSE)
pdf(NULL)

# public_load for min_hiv_reads, hiv_tcell_gate, public_paths and public_maps
code <- chunk_code("public_load")
code <- sub("^public_load_datasets <- .*$", "public_load_datasets <- character(0)", code)
eval(parse(text = code), envir = globalenv())

# Caches the Supp Fig 3 Panel D chunk reads, so they live under extractedData
# rather than beside this script.
spliceviremic.dir <- paste0(data.root, "Figure3/spliceViremic/")
dir.create(spliceviremic.dir, recursive = TRUE, showWarnings = FALSE)

# --- gene sets ----------------------------------------------------------------
art_hallmark_sets <- c("OXIDATIVE_PHOSPHORYLATION", "GLYCOLYSIS", "MTORC1_SIGNALING",
                       "INTERFERON_GAMMA_RESPONSE", "INTERFERON_ALPHA_RESPONSE",
                       "TNFA_SIGNALING_VIA_NFKB", "INFLAMMATORY_RESPONSE")
public_files <- c(
  Sys.glob(paste0(data.root, "Collora/DEG/Viremic_*_DEG_infected.csv")),
  Sys.glob(paste0(data.root, "HIVseq/DEG/Viremic_*_DEG_infected.csv")))

# The SIV splice+ signature: genes called splice+ in the Invitro Ca*Re contrast
# (SuppFigs.Rmd DEG: splice+ vs splice-, p_val_adj < 0.01 and |avg_log2FC| > 0.5).
# W2 and Pacute return no genes at that threshold, so Invitro is the only source.
siv_splice_deg <- read.csv(paste0(paste0(data.root, ""),
                                  "Figure1/combined/DEG/Invitro_DEG_splice.csv")) %>%
  filter(p_val_adj < 0.01, abs(avg_log2FC) > 0.5)
siv_splice_genes <- siv_splice_deg$gene[siv_splice_deg$avg_log2FC > 0]
cat("SIV Invitro splice+ genes:", paste(siv_splice_genes, collapse = ", "), "\n")

gene_sets <- c(cd4t_subset_markers,
               hallmark_gene_sets()[art_hallmark_sets],
               public_hiv_signature(public_files, fc = 0.5, pval = 0.01,
                                    min_tables = 2)["PUBLIC_HIV_UP"],
               list(SIV_SPLICE_UP = siv_splice_genes))

is_viral <- function(genes) genes == "HXB2" | grepl("mac239", genes, ignore.case = TRUE)

# Module scores plus the normalized expression of the SIV splice+ genes. Scores
# are computed over every cell of the object and only then restricted to
# `cells`: AddModuleScore bins genes by mean expression to draw its control
# genes, which a handful of infected cells cannot support.
score_cells <- function(obj, cells) {
  DefaultAssay(obj) <- "RNA"
  if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj, assay = "RNA")
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  host <- !is_viral(rownames(counts))
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
  present <- lapply(gene_sets, intersect, rownames(counts)[host])
  usable <- names(present)[lengths(present) >= 3]
  obj <- AddModuleScore(obj, features = present[usable], assay = "RNA",
                        name = "score", seed = 1, search = FALSE)
  scores <- obj@meta.data[cells, paste0("score", seq_along(usable)), drop = FALSE]
  colnames(scores) <- usable
  genes_present <- intersect(siv_splice_genes, rownames(obj))
  expression <- t(as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[genes_present, cells,
                                                                            drop = FALSE]))
  colnames(expression) <- paste0("gene_", genes_present)
  list(table = bind_cols(tibble(cellID   = cells,
                                nCount   = Matrix::colSums(counts[host, cells, drop = FALSE]),
                                nFeature = Matrix::colSums(counts[host, cells, drop = FALSE] > 0)),
                         as_tibble(scores), as_tibble(expression)),
       sets_found = lengths(present))
}

sets_found <- list()

# --- viremic public datasets --------------------------------------------------
public_cells <- purrr::map_dfr(c("Collora", "HIVseq"), function(dataset) {
  paths <- public_paths[[dataset]]
  libraries <- public_maps[[dataset]] %>% filter(condition == "Viremic")
  purrr::map_dfr(libraries$sample, function(sample) {
    obj <- load_human_sample(sample, paths$counts, paths$rds, qc = public_qc(dataset))
    obj$orig.ident <- sample
    # Each dataset under its own paper's HIV+ rule, as public_library_table and
    # the Figure3 chunks do: "collora" for Collora, "reads" for HIV-seq.
    obj <- attach_hxb2(RenameCells(obj, add.cell.id = sample), sample, paths$bamsort,
                       rule = public_hiv_rule(dataset))
    infected <- colnames(obj)[obj$HIV_status == "HIV+"]
    if (length(infected) < 2) { rm(obj); gc(); return(NULL) }
    md <- obj@meta.data[infected, ]
    scored <- score_cells(obj, infected)
    sets_found[[paste(dataset, sample)]] <<- scored$sets_found
    out <- scored$table %>%
      mutate(dataset = paste0(dataset, "_Viremic"), library = sample,
             participant = libraries$participant[match(sample, libraries$sample)],
             viral_UMIs = md$HXB2_US + md$HXB2_US.SS + md$HXB2_SS.MS + md$HXB2_SS +
                          md$HXB2_MS + md$HXB2_any,
             spliced = md$HXB2_SS.MS + md$HXB2_SS + md$HXB2_MS > 0)
    message(dataset, " ", sample, ": ", nrow(out), " infected, ", sum(out$spliced), " splice+")
    rm(obj); gc()
    out
  })
})

# --- untreated SIV samples ----------------------------------------------------
siv_combined.dir <- paste0(data.root, "Figure1/combined/dataframes/")
siv_rds.dir <- paste0(data.root, "Figure1/singleCell/Seurat/RDS/")
siv_cells <- purrr::map_dfr(c("Invitro", "Pacute", "W2"), function(sample_name) {
  combined <- read.csv(paste0(siv_combined.dir, sample_name, "_all.csv")) %>%
    mutate(cellID = as.character(cellID))
  viral_cols <- setdiff(colnames(combined), c("cellID", "UMAP1", "UMAP2", "cluster", "total"))
  splice_cols <- intersect(siv_splice_cols(combined), colnames(combined))
  combined <- combined %>%
    mutate(viral_UMIs = rowSums(across(all_of(viral_cols))),
           spliced = rowSums(across(all_of(splice_cols))) > 0) %>%
    filter(viral_UMIs > 0)
  obj <- readRDS(paste0(siv_rds.dir, sample_name, ".rds"))
  infected <- intersect(combined$cellID, colnames(obj))
  scored <- score_cells(obj, infected)
  sets_found[[sample_name]] <<- scored$sets_found
  m <- match(scored$table$cellID, combined$cellID)
  out <- scored$table %>%
    mutate(dataset = sample_name, library = sample_name, participant = sample_name,
           viral_UMIs = combined$viral_UMIs[m], spliced = combined$spliced[m])
  message(sample_name, ": ", nrow(out), " infected, ", sum(out$spliced), " splice+")
  rm(obj); gc()
  out
})

cells <- bind_rows(public_cells, siv_cells)
saveRDS(list(cells = cells, gene_sets = gene_sets, siv_splice_genes = siv_splice_genes,
             sets_found = sets_found), paste0(spliceviremic.dir, "splice_viremic_cells.rds"))

cat("\n=== infected cells per dataset ===\n")
print(as.data.frame(cells %>% group_by(dataset) %>%
  summarize(libraries = n_distinct(library), participants = n_distinct(participant),
            infected = n(), splice_pos = sum(spliced),
            median_viral_UMIs = median(viral_UMIs),
            median_viral_UMIs_splice_pos = median(viral_UMIs[spliced]),
            median_viral_UMIs_splice_neg = median(viral_UMIs[!spliced]))), row.names = FALSE)
cat("\n=== genes found per set (min across objects) ===\n")
print(do.call(pmin, sets_found))
