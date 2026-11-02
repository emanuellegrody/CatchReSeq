# Scripts not included in this release

Two scripts are withheld from the public repository and shared on request. This
file records what they are, exactly where the pipeline stops without them, and
what you can still do.

Every script that invokes either one checks for it first and exits with code 2
and a message pointing here, rather than failing obscurely — at submit time, not
buried in a job log.

## The two scripts

### `extractionScripts/pipelineScripts/bamsort/bamsort_alignment_reads.py`

The read-extraction step for the viral arms. Given an aligned BAM and one
genomic window, it counts unique (UMI, cell barcode) combinations in that
window and writes

    UMI, cellID, count, primary_count

one row per molecule per window. `count` includes every read overlapping the
window; `primary_count` charges a read that straddles a boundary to a single
window, so a sum across windows counts it once.

### `extractionScripts/pipelineScripts/bamsort/bamsort_alignment_reads_host.py`

The same step for the host splice targets. At a target exon it separates reads
that contain the exon from reads that skip it, and writes the inclusion and
exclusion counts to two files with the columns

    UMI, cellID, count

## Where the pipeline stops

This is **stepTwo** of the extraction chain. Everything on either side of it is
included:

| step | included? |
|---|---|
| `stepOne_mkcounts.sh` — cellranger | yes |
| `stepTwo_bamsort*.sh` — **the withheld scripts** | **no** |
| `stepThree_process.sh` — UMI collapse | yes |
| `stepFour_minimums.R` — detection thresholds | yes |
| `stepFive_markers.sh` — marker tables for the CD4 T gate | yes |
| `plotScripts/*.Rmd` — every figure | yes |

Nine submission chains reach stepTwo: Figure1 (three arms), Figure3 (two), and
the Giron, HIV and HIV-seq public-data chains use the viral script; Figure4's
host-target arm and the liver arm use the host script.

## What you can still do

**Make the figures.** The notebooks read stepThree and stepFour outputs, not
stepTwo outputs. Those tables are part of the data record deposited with the
paper, not of this repository — see its data availability statement. With them
in place at `$CARE_DATA_ROOT` (see the README), every chunk in
`plotScripts/*.Rmd` runs as written.

**Re-run everything downstream of stepTwo**, from those same deposited tables.

**Re-run cellranger**, from the deposited sequencing data.

What you cannot do without the scripts is re-derive the stepTwo tables from a
BAM. Note that this also requires the BAMs, which are sequencing-derived and not
in this repository either, so the step is only reachable by someone who already
holds that data.

## Requesting them

Open an issue, or email the corresponding author of the paper, saying which
dataset you want to run. Drop the file you are sent at

    extractionScripts/pipelineScripts/<name>.py

and re-run; the wrappers resolve it from their own location and nothing else
needs changing. Both paths are already listed in `.gitignore`, so neither can be
committed back by accident.
