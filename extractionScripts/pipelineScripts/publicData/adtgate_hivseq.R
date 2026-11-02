# CD4+ T cell IDs gated on the TotalSeq-C protein counts inside the published
# GSE266455 matrices, to be used in place of hiv_tcell_gate (gene-only CD3+CD4+,
# which loses half the cells to dropout).
#
# Gate: per marker, an Otsu threshold on CLR values, then
#   CD3+ CD4+ CD8- CD14- CD19- CD20-
# scored against the CD4 T cell counts the paper reports per participant and
# timepoint.
#
# Depth guard: PID145 and most of PID8026 are too shallow for any protein gate
# -- 15 markers over a few hundred ADT UMIs is sampling noise -- so their
# libraries are gated but marked adt_ok = FALSE.
# Nothing downstream should use a library whose flag is FALSE.
suppressMessages({library(Seurat); library(dplyr); library(purrr)})

matrices  <- "/projects/b1042/GoyalLab/egrody/publicData/HIVseq/matrices"
counts    <- "/projects/b1042/GoyalLab/egrody/publicData/HIVseq/counts"
output.dir <- "/projects/b1042/GoyalLab/egrody/extractedData/HIVseq/ADT/"
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

adt_depth_min <- 1000   # median ADT UMIs per cell required to trust the gate

published <- tibble(dir = basename(Sys.glob(file.path(matrices, "PID*")))) %>%
  mutate(PID       = sub("^PID(\\d+)W.*", "\\1", dir),
         week_pub  = sub("^PID\\d+W(\\d+)w.*", "\\1", dir),
         replicate = sub("^PID\\d+W\\d+w(\\d+).*", "\\1", dir),
         chemistry = ifelse(grepl("NOTHING$", dir), "NOTHING", "Primer"))

# Matched on participant + replicate + chemistry, not on the week token: the
# published matrices name PID1052 replicates 2/4/6 as W24, while their real
# suppressed timepoint is week 70 (Frouard et al., Nat Commun 2026, Methods:
# "following ART suppression (Week 24: PID8026; Week 45: PID8027; and Week 70:
# PID1052)"). W24 there is PID8026's week, so the week is taken from our own
# sample names (`week_ours`) and never from the published ones.
ours <- tibble(sample = sub("^GRCh38_HXB2_", "", basename(
  Sys.glob(file.path(counts, "GRCh38_HXB2_PID*"))))) %>%
  mutate(PID       = sub("^PID(\\d+)week.*", "\\1", sample),
         week_ours = sub("^PID\\d+week(\\d+)replicate.*", "\\1", sample),
         replicate = sub(".*replicate(\\d+).*", "\\1", sample),
         chemistry = ifelse(grepl("Primer$", sample), "Primer", "NOTHING"))
pairs <- inner_join(ours, published, by = c("PID", "replicate", "chemistry"))

otsu <- function(x, nbins = 256) {
  h <- hist(x, breaks = nbins, plot = FALSE)
  p <- h$counts / sum(h$counts); mids <- h$mids
  w0 <- cumsum(p); w1 <- 1 - w0
  m0 <- cumsum(p * mids) / ifelse(w0 == 0, NA, w0)
  m1 <- (sum(p * mids) - cumsum(p * mids)) / ifelse(w1 == 0, NA, w1)
  v  <- w0 * w1 * (m0 - m1)^2
  mids[which.max(replace(v, is.na(v), -Inf))]
}
clr_by_feature <- function(m) t(apply(m, 1, function(x) log1p(x) - mean(log1p(x))))

gate_one <- function(i) {
  row <- pairs[i, ]
  data <- Read10X(data.dir = file.path(matrices, row$dir))
  adt  <- as.matrix(data[["Antibody Capture"]])
  rownames(adt) <- sub("_TotalC.*$", "", rownames(adt))
  clr <- clr_by_feature(adt)
  thr <- apply(clr, 1, otsu)
  pos <- function(m) clr[m, ] > thr[[m]]
  cd4t <- pos("CD3") & pos("CD4") & !pos("CD8") &
          !pos("CD14") & !pos("CD19") & !pos("CD20")

  # Barcodes our own cellranger run called, so the list can be intersected
  # without loading the Seurat objects.
  bc_file <- file.path(counts, paste0("GRCh38_HXB2_", row$sample),
                       "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
  our_cells <- readLines(gzfile(bc_file))
  depth <- median(colSums(adt))

  cells <- tibble(sample      = row$sample,
                  published_matrix = row$dir,
                  PID         = row$PID,
                  week        = row$week_ours,
                  chemistry   = row$chemistry,
                  barcode     = colnames(adt)[cd4t]) %>%
    mutate(cellID = paste0(sample, "_", barcode),
           in_our_cellranger_cells = barcode %in% our_cells,
           adt_ok = depth >= adt_depth_min)

  stats <- tibble(sample = row$sample, published_matrix = row$dir,
                  PID = row$PID, week = row$week_ours,
                  chemistry = row$chemistry,
                  adt_umis_median = depth, adt_ok = depth >= adt_depth_min,
                  published_cells = ncol(adt), our_cells = length(our_cells),
                  CD4_T = sum(cd4t),
                  CD4_T_in_our_cells = sum(cells$in_our_cellranger_cells),
                  pct_CD4_T = round(100 * mean(cd4t), 1))
  rm(data, adt, clr); gc()
  list(cells = cells, stats = stats)
}

res   <- map(seq_len(nrow(pairs)), gate_one)
cells <- bind_rows(map(res, "cells"))
stats <- bind_rows(map(res, "stats"))

cat("=== per library ===\n")
print(as.data.frame(stats %>% arrange(PID, week, chemistry)))

targets <- tribble(
  ~PID,   ~group,     ~published_CD4_T,
  "1052", "W0 NOTHING",  4263, "1052", "W0 Primer",   8772,
  "8027", "W0 NOTHING",  8022, "8027", "W0 Primer",   9701,
  "8026", "W0 Primer",  35849, "145",  "W0 Primer",  19060,
  "1052", "ART",        15468, "8026", "ART",        42618,
  "8027", "ART",        17072)

cat("\n=== vs published CD4 T counts ===\n")
print(as.data.frame(stats %>%
  mutate(group = ifelse(week == "0", paste("W0", chemistry), "ART")) %>%
  group_by(PID, group, adt_ok) %>%
  summarize(libraries = n(), CD4_T = sum(CD4_T), .groups = "drop") %>%
  left_join(targets, by = c("PID", "group")) %>%
  mutate(pct_of_published = round(100 * CD4_T / published_CD4_T, 1)) %>%
  arrange(PID, group)))

write.csv(cells, file.path(output.dir, "CD4T_cellIDs.csv"), row.names = FALSE)
write.csv(stats, file.path(output.dir, "CD4T_gate_summary.csv"), row.names = FALSE)
cat("\nwrote", nrow(cells), "cell IDs to",
    file.path(output.dir, "CD4T_cellIDs.csv"), "\n")
cat("DONE\n")
