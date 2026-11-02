# Catch*ReSeq

Code to reproduce analysis and figures in the Catch*ReSeq manuscript (Grody et al., 2026)

## Layout

```
extractionScripts/
  pipelineScripts/      worker scripts
    cellranger/         cellranger count, one worker per library type
    bamsort/            BAM file read extraction
      coordinates/      the per-genome target region tables for stepTwo
    polya/              internal vs terminal priming (Figure2)
    isoform/            host exon ratios: BRIE2 (Figure4)
    cnmf/               cNMF factorization (Figure3)
    fish/               smFISH quantification (Figure2)
    memoryseq/          clone bulk RNA-seq and the heritable-gene analysis (Figure4)
    publicData/         re-analyzed public datasets: download, ADT gates,
                        label fixes, and downstream analyses
  submissionScripts/    SLURM wrappers and their sample CSVs, one dir per figure
plotScripts/
  Figure1-4.Rmd         one notebook per main figure
  SuppFigs.Rmd          every supplementary panel
  functions.R           analysis and plotting functions for notebooks
  packages.R            package dependency install and versions
  config.R              the four filesystem roots -- update upon download
  referenceData/        third-party reference files
```


## Configuring where the data lives

`plotScripts/config.R` defines four roots, each overridable by an environment
variable:

| variable | default | holds |
|---|---|---|
| `CARE_DATA_ROOT` | `…/extractedData/` | pipeline outputs, per figure |
| `CARE_PLOT_ROOT` | `…/plotData/` | the figure PDFs |
| `CARE_GENOME_ROOT` | `…/genomes/` | references and indices |
| `CARE_PUBLIC_ROOT` | `…/publicData/` | re-analysed published datasets |

The shell side takes `DATA_ROOT` and `GENOME_ROOT` for the same purpose.
Nothing else needs changing to run.

Source the files in this order:

```r
source("packages.R")
source("config.R")
source("functions.R")
```

## Reproducing a figure

1. Point the roots at your data, per the table above.
2. Open the notebook for that figure in `plotScripts/`. Each top-level section
   has an Extraction half and a Plotting half, and its prose names the
   prerequisite scripts. Extraction needs a cluster; plotting does not.

The extraction chain for the single-cell arms runs in order:

| step | writes |
|---|---|
| `stepOne_mkcounts.sh` | cellranger counts per library |
| `stepTwo_bamsort.sh` | per-region alignment tables from the BAM |
| `stepThree_process.sh` | UMI-collapsed read tables |
| `stepFour_minimums_CaRe.R` | detection-threshold filtered tables |
| `stepFive_markers.sh` | the per-cell marker tables the CD4 T gate reads |

`runAll.sh` in `submissionScripts/Figure1/` drives the whole chain for the SIV
arm and waits between steps.

Figure4 clonal expansion runs separately, driven by Nextflow:

| step | writes |
|---|---|
| `memoryseq_nfcore_rnaseq.sh` | nf-core/rnaseq gene TPM matrix |
| `memoryseq_qc.sh` | per-library QC, PDFs only |
| `memoryseq_variability.sh` | the over-dispersed-gene list Panel J reads |

## Dependencies

- R 4.4 with Seurat v5. `plotScripts/packages.R` installs everything and
  records resolved versions to `packages_versions.csv`; it attaches a core set
  and defers the Bioconductor genomics stack until a script asks for it.
- Python 3 with `pysam`, `numpy`, `pandas`, `scipy`, `scanpy`, `cnmf`,
  `gffutils`, `bigfish`, `tifffile`, `matplotlib`.
- cellranger (7.2 or 10.0; set `CELLRANGER_BIN`), STAR, samtools, seqkit,
  BRIE2.
- SLURM. The wrappers call `sbatch` and are written for an account/partition
  pair you will want to change at the top of each.

## bamsort script availability

The stepTwo read-extraction scripts, `bamsort_alignment_reads.py` and
`bamsort_alignment_reads_host.py`, are withheld and shared on request.

