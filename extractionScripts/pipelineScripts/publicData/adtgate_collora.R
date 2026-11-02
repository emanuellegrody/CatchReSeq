# CD4+ T cell IDs gated on the TotalSeq-C protein counts inside the published
# GSE187515 matrices, using the threshold Collora et al. 2022 define themselves:
#
#   "cells having higher normalized read counts of barcoded surface antibody
#    staining, as measured by >90th percentile of the expression level of
#    isotype barcoded antibody controls, were defined as positive for the
#    respective surface protein expression."
#    -- STAR Methods, Surface protein expression analysis
#
# The three isotype controls they stained (mouse IgG1/IgG2a/IgG2b) are carried
# in the matrices, so the rule is applied as written: one threshold per library
# at the 90th percentile of the pooled isotype values, a marker is positive in a
# cell above it. Being referenced to each library's own background, the
# threshold needs no ADT depth guard -- unlike the HIVseq gate, whose panel has
# no isotype control and so uses an Otsu cut plus adt_depth_min.
#
# Normalization is CLR across features within each cell, NOT the per-feature CLR
# the HIVseq gate uses: per-feature CLR centres every feature on its own mean,
# which would put isotypes and real markers on the same centre and destroy the
# comparison this rule depends on.
#
# The gate is CD3+ AND CD4+, with NO protein CD8 exclusion. Per library, CD8
# sits at the isotype controls' own level, so its distribution overlaps the
# background the threshold is drawn from and no isotype percentile separates it.
# That is the expected reading of a negative selection: these cells were
# isolated as CD4+ T, so there are no CD8 T cells left and the CD8 ADT is
# background. The paper does not use protein CD8 to exclude either --
# it discards "clusters with expression of CD8A, KLRC1, or KLRC3 in more than
# 20% of cells" by gene expression at the cluster level, and uses the isotype
# rule only to call markers positive, which is how it is used here.
#
# Panels also differ by donor (STAR Methods, DNA-barcoded surface protein
# staining): donors 236/829/739/799 got the 24-marker panel, donors 640/910 a
# 13-marker panel. CD3 and CD4 are in both, so the gate is the same for all
# libraries; the panel difference is reported but does not change the gate.
#
# The matrices are raw (737,280 barcodes), so everything is computed after
# subsetting to the barcodes our own cellranger run called.
suppressMessages({library(Matrix); library(dplyr); library(purrr); library(tibble)})

matrices   <- "/projects/b1042/GoyalLab/egrody/publicData/Colloraetal2022/processed/matrices"
counts     <- "/projects/b1042/GoyalLab/egrody/publicData/Colloraetal2022/counts"
output.dir <- "/projects/b1042/GoyalLab/egrody/extractedData/Collora/ADT/"
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

isotypes      <- c("MouseIgG1", "MouseIgG2a", "MouseIgG2b")
isotype_pctl  <- 0.90   # the paper's ">90th percentile of isotype" rule
cd3_pos_min   <- 50     # library-level sanity check, see below

# Donors stained without CD8 (13-marker panel)
donors_no_cd8 <- c("640", "910")

# CLR within each cell, across features: x -> log1p(x / geometric mean(x))
clr_by_cell <- function(m) {
  logm <- log1p(m)
  gm   <- colMeans(logm)
  t(t(logm) - gm)
}

gate_one <- function(sample) {
  d    <- file.path(matrices, paste0(sample, "matrix"))
  feat <- read.delim(gzfile(file.path(d, "features.tsv.gz")), header = FALSE)
  bc   <- readLines(gzfile(file.path(d, "barcodes.tsv.gz")))
  m    <- readMM(gzfile(file.path(d, "matrix.mtx.gz")))

  is_adt <- feat$V3 == "Antibody Capture"
  adt <- m[is_adt, , drop = FALSE]
  rownames(adt) <- sub("_Total.*$", "", feat$V2[is_adt])
  colnames(adt) <- bc

  # Restrict to our own called cells before normalizing: on the raw matrix the
  # isotype percentile would be taken over empty droplets.
  bc_file <- file.path(counts, paste0("GRCh38_HXB2_", sample),
                       "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
  our_cells <- readLines(gzfile(bc_file))
  adt <- as.matrix(adt[, intersect(bc, our_cells), drop = FALSE])

  clr <- clr_by_cell(adt)
  thr <- quantile(as.vector(clr[isotypes, , drop = FALSE]), isotype_pctl,
                  na.rm = TRUE)
  pos <- function(marker) clr[marker, ] > thr

  participant <- sub("^[A-Za-z]+", "", sample)
  cd4t <- pos("CD3") & pos("CD4")

  cells <- tibble(sample      = sample,
                  participant = participant,
                  barcode     = colnames(adt)[cd4t]) %>%
    mutate(cellID = paste0(sample, "_", barcode),
           in_our_cellranger_cells = TRUE,
           adt_ok = 100 * mean(pos("CD3")) >= cd3_pos_min)

  # CD8 and the isotype medians are carried as diagnostics only: they are what
  # rules the protein CD8 exclusion out, so they should stay visible.
  stats <- tibble(sample = sample, participant = participant,
                  panel_has_CD8 = !participant %in% donors_no_cd8,
                  our_cells = length(our_cells), gated_cells = ncol(adt),
                  adt_umis_median = median(colSums(adt)),
                  isotype_threshold = round(unname(thr), 3),
                  CD8_median_raw = median(adt["CD8", ]),
                  isotype_median_raw = median(adt[isotypes, ]),
                  pct_CD3_pos = round(100 * mean(pos("CD3")), 1),
                  pct_CD4_pos = round(100 * mean(pos("CD4")), 1),
                  pct_CD8_pos_unused = round(100 * mean(pos("CD8")), 1),
                  CD4_T = sum(cd4t),
                  pct_CD4_T = round(100 * mean(cd4t), 1)) %>%
    mutate(adt_ok = pct_CD3_pos >= cd3_pos_min)

  rm(m, adt, clr); gc(FALSE)
  list(cells = cells, stats = stats)
}

samples <- sub("matrix$", "", basename(Sys.glob(file.path(matrices, "*matrix"))))
res   <- map(samples, gate_one)
cells <- bind_rows(map(res, "cells"))
stats <- bind_rows(map(res, "stats"))

cat("=== per library ===\n")
print(as.data.frame(stats %>% arrange(participant, sample)), row.names = FALSE)

# Adapted from the paper's own QC, which discarded clusters with <50% CD3
# surface staining as non-T. Applied per library here only to catch a library
# whose ADT failed outright; these are negatively selected CD4+ T cells, so a
# library well under 50% CD3+ would mean the staining, not the biology.
cat("\n=== libraries failing the >=", cd3_pos_min, "% CD3+ check ===\n")
failed <- stats %>% filter(!adt_ok)
if (nrow(failed) == 0) cat("none\n") else
  print(as.data.frame(failed %>% select(sample, pct_CD3_pos, adt_umis_median)),
        row.names = FALSE)

write.csv(cells, file.path(output.dir, "CD4T_cellIDs.csv"), row.names = FALSE)
write.csv(stats, file.path(output.dir, "CD4T_gate_summary.csv"), row.names = FALSE)
cat("\nwrote", nrow(cells), "cell IDs to",
    file.path(output.dir, "CD4T_cellIDs.csv"), "\n")
cat("DONE\n")
