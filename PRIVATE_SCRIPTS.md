# Scripts not included in this release

Please note that two essential Catch*ReSeq pipeline scripts are withheld from
the public repository and may be shared on request. These scripts perform isoform 
assignment at stepTwo. The pipeline will not run without them; the pipeline exits
with code 2 and a message pointing here if absent.

## The two scripts

### `extractionScripts/pipelineScripts/bamsort/bamsort_alignment_reads.py`

Catch*ReSeq reads-extraction script for viral targets. Given an aligned BAM
and one genomic window, it counts unique (UMI, cell barcode) combinations in
that window and writes

    UMI, cellID, count, primary_count

one row per molecule per window. `count` includes every read overlapping the
window; `primary_count` charges a read that straddles a boundary to a single
window to avoid duplication.

### `extractionScripts/pipelineScripts/bamsort/bamsort_alignment_reads_host.py`

Catch*ReSeq reads-extraction script for host targets. At a target exon, it
separates reads that contain the exon from reads that skip it, and writes the
inclusion and exclusion counts to two files with the columns

    UMI, cellID, count

## Affected pipeline step

The withheld scripts affect only **stepTwo** of the extraction chain:

| step | included? |
|---|---|
| `stepOne_mkcounts.sh` — cellranger | yes |
| `stepTwo_bamsort*.sh` — **the withheld scripts** | **no** |
| `stepThree_process.sh` — UMI collapse | yes |
| `stepFour_minimums.R` — detection thresholds | yes |
| `stepFive_markers.sh` — marker tables for the CD4 T gate | yes |
| `plotScripts/*.Rmd` — every figure | yes |

## Unaffected processes

**Re-running everything downstream of stepTwo** including figure plotting.
The stepTwo outputs are deposited on [Google Drive](https://drive.google.com/drive/folders/1VZqFvq4PPePjBebWIqZQCKbCONigd2-u?usp=sharing) (see the README).

**Re-running cellranger**, from the deposited sequencing data.

## Requesting withheld scripts

Please contact the corresponding author at:
[emanuelle.grody@northwestern.edu](mailto:emanuelle.grody@northwestern.edu).
