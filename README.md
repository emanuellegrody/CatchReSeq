# Catch*ReSeq

Scripts to reproduce every figure in the Catch*ReSeq manuscript: targeted
capture and re-sequencing of viral transcripts from single-cell libraries,
applied to SIV in rhesus macaque CD4+ T cells, HIV in human donor samples, and
host splice isoforms in mouse liver.

The repository is self-contained. Nothing here calls a script outside it, and
every filesystem location is set in one place.

## Layout

```
extractionScripts/
  pipelineScripts/      workers -- the .py and .R that do the work, by stage
    cellranger/         cellranger count, one worker per library type
    bamsort/            read extraction and the stepTwo-stepFive chain
      coordinates/      the per-genome target region tables stepTwo reads
    polya/              internal vs terminal priming (Figure2)
    isoform/            host cassette-exon ratios: BRIE2 (Figure4)
    cnmf/               cNMF factorisation and its .h5ad input (Figure3)
    fish/               smFISH quantification (Figure2)
    memoryseq/          clone bulk RNA-seq and the heritable-gene call (Figure4)
    publicData/         re-analysed published datasets: download, ADT gates,
                        label fixes, and the analyses specific to them
  submissionScripts/    SLURM wrappers and their sample CSVs, one dir per figure
plotScripts/
  Figure1-4.Rmd         one notebook per main figure
  SuppFigs.Rmd          every supplementary panel
  functions.R           the analysis and plotting library the notebooks source
  packages.R            dependency install and version pinning
  config.R              the four filesystem roots -- the only thing to edit
  referenceData/        third-party reference files the notebooks read directly
```

`pipelineScripts/` is grouped by pipeline stage rather than by figure, because
a worker is usually shared across figures -- `stepThree_process_worker.sh` is
driven by ten wrappers spanning three figures and two public datasets.

A worker never hardcodes a path: `submissionScripts/<Figure>/*.sh` sets the
environment and submits, and the worker in `pipelineScripts/` reads it. Both
resolve each other from their own location, so the repository can be cloned
anywhere.

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
Nothing else needs changing to run this against your own copy of the data.

`config.R` must be sourced **between** `packages.R` and `functions.R`, because
`functions.R` builds default arguments from `data.root`:

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

The LD002 clone panel in Figure 4 is bulk RNA-seq and does not join that chain.
It runs on its own, and is the one stage driven by Nextflow rather than by
`sbatch` directly:

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

## Two scripts are not included

The stepTwo read-extraction scripts, `bamsort_alignment_reads.py` and
`bamsort_alignment_reads_host.py`, are withheld and shared on request. Every
other step is here, and the figures are built from stepThree/stepFour tables
rather than stepTwo output, so the plotting path is unaffected — see
[PRIVATE_SCRIPTS.md](PRIVATE_SCRIPTS.md) for exactly where the chain stops.

## Data

The sequencing data and the extracted intermediates are not in this repository;
see the data availability statement in the paper. The intermediates the notebooks
read are the stepThree and stepFour tables, which is what makes the plotting path
independent of the two withheld scripts.

[DATA_INPUTS.md](DATA_INPUTS.md) lists every raw input the extraction scripts
read -- which accessions are ours to deposit, which are published and only
cited, which references have to be rebuilt, and the arms whose FASTQs are no
longer on disk.
