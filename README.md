# Catch\*ReSeq

Scripts to reproduce every figure in the Catch\*ReSeq manuscript
(DOI: [10.64898/2026.11.12.763766](https://doi.org/10.64898/2026.11.12.763766)).

- **Demo:** [`demo/`](demo/) -- runs in about 5 seconds, no cluster needed.
- **What is not included:** [two read-extraction scripts](PRIVATE_SCRIPTS.md) -- reanalysis of the deposited data does not require these scripts.
- **License:** MIT, see [LICENSE](LICENSE).

---

## Layout

```
demo/                     a real data sample
environments/             Python dependency manifests, one per conda env
extractionScripts/
  pipelineScripts/        worker scripts
    cellranger/           cellranger count, one worker per library type
    bamsort/              read extraction and the stepTwo-stepFive pipeline
      coordinates/        the per-genome target region tables for stepTwo
    polya/                internal vs terminal polyA usage (Figure2)
    isoform/              host cassette-exon ratios: BRIE2 (Figure4)
    cnmf/                 cNMF factorization code and input (Figure3)
    fish/                 smFISH quantification (Figure2)
    memoryseq/            clone bulk RNAseq and heritable gene analyses (Figure4)
    publicData/           published datasets reanalysis scripts
  submissionScripts/      SLURM wrappers and sample CSVs, one dir per figure
plotScripts/
  Figure1-4.Rmd           one notebook per main figure
  SuppFigs.Rmd            every supplementary panel
  functions.R             all notebook analysis and plotting functions
  packages.R              dependency install and versions
  config.R                filesystem roots -- update
  referenceData/          third-party reference files
```

Please note that the stepTwo read-extraction scripts are withheld from this
public repository. For more information, see [PRIVATE_SCRIPTS.md](PRIVATE_SCRIPTS.md).

---

## 1. System requirements

### Operating system

Developed and run on Red Hat Enterprise Linux 8.10 (Linux 4.18.0-553).
The R plotting path and the demo are platform-independent; the extraction
path assumes a Linux cluster with SLURM.

### Software dependencies

| | version | notes |
|---|---|---|
| R | **4.4.0** | `module load R/4.4.0` throughout |
| R packages | 47, pinned | installed and version-locked by `plotScripts/packages.R`; the resolved set is in [`plotScripts/packages_versions.csv`](plotScripts/packages_versions.csv) -- Seurat 5.5.0, data.table 1.18.2.1, ggplot2 4.0.3, SingleCellExperiment 1.28.1, zellkonverter 1.16.0, harmony 2.0.2 |
| Python | **3.13.11** and **3.6.15** | two environments -- see [`environments/`](environments/). The published analysis ran in four; those manifests are kept there too |
| cellranger | **7.2.0** and **10.0.0** | chemistry-dependent; set `CELLRANGER_V3_BIN` / `CELLRANGER_V4_BIN`. Obtained separately from 10x Genomics under their own license |
| samtools | any recent | `module load samtools` |
| seqkit | any recent | FASTQ subsampling only |
| sra-toolkit | **3.2.1** | public-data download only |
| Nextflow + Singularity | Nextflow **25.10.4**, OpenJDK **21.0.10**, Singularity **3.8.1** | two arms. `memoryseq_nfcore_rnaseq.sh` runs nf-core/rnaseq 3.17.0 for the bulk clone data; `longread_scnanoseq.sh` runs nf-core/scnanoseq 1.0.0 for the Figure 2 long-read arm. Both need a conda env holding Nextflow -- see [`environments/README.md`](environments/README.md) |
| SLURM | any | extraction path only |

smFISH needs an environment of its own: big-fish 0.6.2 calls numpy aliases
removed in numpy 1.24 and imports a scikit-image module removed in 0.19.
Everything else -- read extraction, splice, polyA, cNMF and BRIE2 -- runs in one
python 3.13 environment. [`environments/README.md`](environments/README.md) says
which scripts activate which, and records the four environments the published
analysis actually used.

### Tested on

The published analysis was produced on:

- Red Hat Enterprise Linux 8.10, x86_64
- R 4.4.0 with the package versions in `packages_versions.csv`
- Python 3.9.18 (extraction), 3.13.11 (cNMF), 3.11.14 (BRIE2), 3.6.15 (smFISH)
  -- four environments, not the two the installation guide recommends
- cellranger 7.2.0 and 10.0.0
- SLURM on the Northwestern Quest cluster (`b1042` / `genomics`)

The demo additionally runs on a stock R 4.4 install with only `data.table`,
`dplyr` and `tidyr`.

### Non-standard hardware

**The extraction path needs a cluster.** The largest job allocation is
128 GB RAM, 24 h walltime, 16 cores.

Storage: roughly 1.2 TB for a full re-run -- about 646 GB of sequencing,
469 GB of imaging, 70 GB of references, plus intermediates.

---

## 2. Installation guide

```sh
git clone https://github.com/emanuellegrody/CatchReSeq.git
cd CatchReSeq
```

### R packages

```sh
cd plotScripts
Rscript -e 'source("packages.R")'
```

`packages.R` installs every missing package with correct versions.

### Python environments

Two: one for smFISH, one for everything else.

```sh
conda create -n care python=3.13.11 && conda activate care
pip install -r environments/care.txt

conda create -n bigfish_env python=3.6.15 && conda activate bigfish_env
pip install -r environments/bigfish.txt
```

The workers activate `cellranger`, `jupyter` and `brie2` by name, so either
build `care.txt` under each of those three names or change the `source activate`
line in the workers. [`environments/README.md`](environments/README.md) maps
names to scripts and lists the four environments the published analysis ran in,
should you want to reproduce those versions exactly.

### cellranger

Download 7.2.0 and 10.0.0 from 10x Genomics and point `CELLRANGER_V3_BIN` and
`CELLRANGER_V4_BIN` at the binaries. Needed only to re-run stepOne.

### HIVToolbox

The viral genome diversity analyses used for Catch*ReSeq SI 1E and SI 3F were 
produced with HIVToolbox, a third-party MATLAB package described in Li et al., 2015 
(DOI: [10.1186/s12977-015-0148-6](https://doi.org/10.1186/s12977-015-0148-6)) and recorded in 
[`plotScripts/referenceData/README.md`](plotScripts/referenceData/README.md).

### Typical install time

The R dependencies include the 47 packages in `packages_versions.csv` plus 
their transitive dependencies (224 total). On macOS or Windows,
wbudget 20-40 minutes. On Linux, installation will take 1-3 hours
on a single core; pass `Ncpus` to shorten it.

Installing a Python environment takes 1-3 minutes. cellranger is a ~10 GB
download.

---

## 3. Demo

```sh
bash demo/run_demo.sh
```

Runs stepThree and stepFour on the example stepTwo output in `demo/input/` 
and compares the result against the output published with the manuscript 
("D195" ART sample).

**Expected output**, written to a scratch directory and diffed against
`demo/expected_output/`:

| file | contents |
|---|---|
| `D195_UMI_read_counts_raw.csv` | 4,062 UMIs from 3,201 cells |
| `D195_UMI_read_counts_full.csv` | 2,333 UMIs from 1,745 called cells |
| `minimum/D195_CaRe_filtered.csv` | 52 cells past the detection thresholds |
| `minimum/D195_filters.txt` | the thresholds used |

ending in `== demo passed`.

**Expected run time: about 5 seconds** on a normal desktop, single core, once 
dependencies are installed.

[`demo/README.md`](demo/README.md) explains the data, the two steps, and how to
point the same steps at your own data.

---

## 4. Instructions for use

### Running the pipeline

`plotScripts/config.R` defines four roots, each overridable by an environment
variable. Nothing else needs changing.

| variable | default | holds |
|---|---|---|
| `CARE_DATA_ROOT` | `…/extractedData/` | pipeline outputs, per figure |
| `CARE_PLOT_ROOT` | `…/plotData/` | the figure PDFs |
| `CARE_GENOME_ROOT` | `…/genomes/` | references and indices |
| `CARE_PUBLIC_ROOT` | `…/publicData/` | re-analysed published datasets |

The shell side takes `DATA_ROOT` and `GENOME_ROOT` for the same purpose, and the
`publicData/` scripts read `CARE_DATA_ROOT` and `CARE_PUBLIC_ROOT` directly. One
script needs a fifth root: `hcr_positivity_overlay.py` reads the raw imaging
stacks, under `CARE_RAW_ROOT` (default `…/rawData/`).

Source files for plotting notebooks must be called in this order:

```r
source("packages.R")
source("config.R")
source("functions.R")
```

### Start to finish on a new machine

**1. Make one directory to hold the data.** Every root in the table above sits
inside it, so one block of `export`s configures the whole pipeline. Allow 2 GB
for the deposit. The optional public-data re-fetch in step 3 goes in the same
directory.

```sh
mkdir -p ~/care
```

**2. Clone the repository.** Anywhere -- it reads its own location and does not
need to sit beside the data. 15 MB.

```sh
git clone https://github.com/emanuellegrody/CatchReSeq.git ~/CatchReSeq
cd ~/CatchReSeq
```

**3. Download the deposit into the data directory.** It is in two places. Until
the manuscript is published, take it from the Google Drive folder:

<https://drive.google.com/drive/folders/1DAGRXQcvU5kC3DeqRspy2HYSW_7QHFoU?usp=sharing>

In a browser, use *Download* on the folder itself rather than on each
subfolder, so the per-figure tree is preserved; Drive serves it as a zip:

```sh
unzip ~/Downloads/<downloaded>.zip -d ~/care
```

From the command line, [`gdown`](https://github.com/wkentaro/gdown)
(`pip install gdown`) fetches the same folder without a browser:

```sh
cd ~/care && gdown --folder \
  'https://drive.google.com/drive/folders/1DAGRXQcvU5kC3DeqRspy2HYSW_7QHFoU'
```

From publication, the same files are the GEO supplementary archive:

```sh
tar -xzf GSE<accession>_CatchReSeq_processed.tar.gz -C ~/care
```

Whichever route you take, what matters is where it lands. The pipeline expects
`~/care/extractedData/Figure1/...`, so if the download leaves you holding a
directory whose immediate contents are `Figure1/`, `Figure3/` and `Figure4/`,
move it into place as `~/care/extractedData` before going on:

```sh
ls ~/care/extractedData/Figure1   # should list CaRe/ and singleCell/
```

Either route gives you the cellranger count matrices, the stepTwo alignment
tables, the splice correction tables and the `raw_cellIDs/` caches, laid out per
figure.

The three re-analysed published datasets are fetched from the original accessions;
they are needed only for the Figure 3 and Supp Fig 3 public-data panels:

```sh
cd ~/CatchReSeq/extractionScripts/submissionScripts/PublicData
DATA_ROOT=~/care bash runAll_publicData.sh download
```

This will need a cluster (the downloads are `--array` jobs). The accession
lists and the ENA manifest ship with the repository, beside each dataset's
`samples_manifest.csv` under `submissionScripts/PublicData/`, so the clone is
all the download stage needs.

**4. Install the dependencies.** Section 2 above. R and the three packages
`bamsort_process.R` needs are enough for steps 5 and 7; the full
`plotScripts/packages.R` set is needed to open a notebook, and the Python
environment only for the cNMF, BRIE2, polyA and smFISH panels.

**5. Check the install with the demo.**

```sh
bash demo/run_demo.sh
```

About 5 seconds, ending in `== demo passed`. If it does, stepThree and stepFour
will run on the deposit.

**6. Point the roots at the download.**

```sh
export CARE_DATA_ROOT=~/care/extractedData
export CARE_PLOT_ROOT=~/care/plotData
export CARE_PUBLIC_ROOT=~/care/publicData
export CARE_GENOME_ROOT=~/care/genomes
export DATA_ROOT=~/care
```

Note the asymmetry: the four `CARE_*` variables name their own directory, while
the shell-side `DATA_ROOT` names the directory *above* `extractedData`. The two
stepFour scripts read either one.

**7. Run stepThree and stepFour over the deposited tables.** The notebooks read
their output, not the deposited stepTwo tables, so this comes before opening
anything in `plotScripts/`.

On a SLURM cluster, `runAll.sh` takes the step to start from, so stepOne and the
withheld stepTwo are skipped:

```sh
cd extractionScripts/submissionScripts/Figure1 && mkdir -p logs
bash runAll.sh stepThree          # SIV Ca*Re: stepThree, waits, then stepFour
bash runAll_TenX.sh stepThree     # the matched 10x libraries
```

The other arms have no `runAll`; submit their wrappers directly, from their own
directory so the relative `logs/` resolves:

| arm | wrapper |
|---|---|
| Figure 3 Ca\*Re and 10x | `Figure3/stepThree_process_CaRe.sh`, `Figure3/stepThree_process.sh` |
| Figure 4 host targets | `Figure4/stepThree_process.sh` |
| Figure 4 SIV | `Figure4/stepThree_process_SIV.sh`, `stepThree_process_SIV_TenX.sh` |
| Collora et al. 2022 | `PublicData/HIV/stepThree_process.sh` |
| HIV-seq / Frouard et al. | `PublicData/HIVseq/stepThree_process.sh` |

Without SLURM, those wrappers cannot be used: each one builds an environment
block and hands it to `sbatch`, and the worker it submits begins with
`module purge`. Run the R script directly with the same block instead, one group
at a time -- `demo/run_demo.sh` is that call with real values in it, and the
`export` lines in any `stepThree_process.sh` are the values for the other
groups. stepFour needs no wrapper at all:

```sh
Rscript extractionScripts/pipelineScripts/bamsort/stepFour_minimums_CaRe.R
Rscript extractionScripts/pipelineScripts/bamsort/stepFour_minimums_TenX.R
```

**8. stepFive, only for the panels that use the CD4 T gate.**
`stepFive_markers.sh` reads the Seurat objects rather than the tables above,
and the deposit does not include them; regenerate them from the deposited
matrices with `SeuratPipeline()` first, per the Data section below.

**9. Open the notebook for the figure.** From `plotScripts/`, with the three
`source()` calls above at the top, then run the Plotting half of the section you
want. Panels land under `CARE_PLOT_ROOT`.

Going back further than stepTwo is a different proposition: it needs the SRA
reads (634 GB), the cellranger references built by `mkref_GRCh38_HXB2.sh` and
the `cellranger/` wrappers, cellranger itself from 10x, and the withheld
read-extraction script -- see [PRIVATE_SCRIPTS.md](PRIVATE_SCRIPTS.md).

### SLURM account and partition

Wrappers are written for the Northwestern Quest allocation `b1042` /
`genomics`. Scripts that submit with `sbatch --wrap` or build their own
`srun` line take `SLURM_ACCOUNT` and `SLURM_PARTITION`:

```sh
SLURM_ACCOUNT=myacct SLURM_PARTITION=mypart bash cNMF_workflow.sh
```

Scripts with their own `#SBATCH -A b1042 -p genomics` header would need to
be edited directly; run `grep -rn 'SBATCH -A\|SBATCH -p' extractionScripts/`
to find them. The `#SBATCH --output=` paths are relative, so `mkdir -p logs` 
in the directory you submit from.

### Reproducing a figure

1. Point the roots at your data, per the table above.
2. Open the notebook for that figure in `plotScripts/`. Each top-level section
   has an Extraction half and a Plotting half, and its prose names the
   prerequisite scripts. Extraction needs a cluster; plotting does not.

| notebook | covers |
|---|---|
| `Figure1.Rmd` | SIV arm: KLRB1 control, Ca\*Re vs 10X detection, Venn panels, UMAPs |
| `Figure2.Rmd` | polyA priming, splice maps, smFISH, long-read isoforms |
| `Figure3.Rmd` | cNMF patterns, viral load and splice state, ART comparison |
| `Figure4.Rmd` | host cassette exons, clonal expansion data, heritable genes, liver |
| `SuppFigs.Rmd` | every supplementary panel, grouped by `# Supp Fig N` |

Within each notebook, panels are named in the Markdown heading above the chunk
that draws them (e.g. `### Panel F: splice+ and viral transcript load`). A few
Supp Fig panels are produced from a main-figure notebook, and say so in the
heading.

The extraction pipeline runs in the following order:

| step | writes |
|---|---|
| `stepOne_mkcounts.sh` | cellranger counts per library |
| `stepTwo_bamsort.sh` | per-region alignment tables from the BAM |
| `stepThree_process.sh` | UMI-collapsed read tables |
| `stepFour_minimums_CaRe.R` | detection-threshold filtered tables |
| `stepFive_markers.sh` | the per-cell marker tables the CD4 T gate reads |

`runAll.sh` in `submissionScripts/Figure1/` runs the whole pipeline for the SIV
arm and waits between steps.

The clonal expansion data in Figure 4 is processed separately by Nextflow:

| step | writes |
|---|---|
| `memoryseq_nfcore_rnaseq.sh` | nf-core/rnaseq gene TPM matrix |
| `memoryseq_qc.sh` | per-library QC, PDFs only |
| `memoryseq_variability.sh` | heritable-gene list for plotting |

---

## Data

### What's included in the Google Drive link

The deposit is the [Google Drive folder](https://drive.google.com/drive/folders/1DAGRXQcvU5kC3DeqRspy2HYSW_7QHFoU?usp=sharing) until publication
and the GEO supplementary archive after it; section 4, step 3 has both download
routes. It holds the **cellranger count matrices**
and the **stepTwo alignment tables** -- the output of the withheld
read-extraction script, which is the earliest point the pipeline can be picked
up from without it. With those in place at `$CARE_DATA_ROOT`, the extraction
chain runs from stepThree onward using the scripts in this repository, and
every panel built on its output reproduces: the Ca\*Re detection and threshold
path, the per-sample UMI and category tables, and the Venn and detection
comparisons that rest on them.

Note the order this implies. The notebooks read stepThree and stepFour output,
not stepTwo, so run `stepThree_process.sh` and `stepFour_minimums_CaRe.R`
against the deposited tables before opening a notebook. `demo/` is that
sequence on one sample, with the expected output to check against.

Regenerating the Seurat objects from the deposited count matrices is possible
with `SeuratPipeline()` in `plotScripts/functions.R`. Note that clustering is
not bit-reproducible, so cluster labels may differ.

---

## License and citation

MIT, see [LICENSE](LICENSE). The data deposits carry their own terms; cellranger
is licensed separately by 10x Genomics.

If you use this code, please cite the Catch\*ReSeq manuscript (Grody et al.,
2026; (DOI: [10.64898/2026.11.12.763766](https://doi.org/10.64898/2026.11.12.763766))).
For questions about the withheld scripts, please contact the corresponding author at:
[emanuelle.grody@northwestern.edu](mailto:emanuelle.grody@northwestern.edu).
