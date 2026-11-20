# Demo

Runs stepThree and stepFour extraction steps on the dataframes generated from stepTwo 
and checks the result against the output reported in the manuscript.

```sh
bash demo/run_demo.sh
```

Because the files are small, you do not need to use the cluster. Pass a directory
argument to keep the intermediate tree: `bash demo/run_demo.sh /tmp/caredemo`.

## What the data is

The ART-suppressed sample from Figure 1: rhesus macaque CD4+ T cells on
suppressive antiretroviral therapy, five Catch\*ReSeq amplicon libraries (`env`,
`SSenv`, `pol`, `D1`, `tat` aka D4) across the five SIVmac239 genome windows.

| `input/` | what it is | written by |
|---|---|---|
| `alignment/` (25 CSV) | per-window read tables, `UMI, cellID, count, primary_count` | stepTwo |
| `alignment/*_alignment_summary.txt` (5) | per-library window summary | stepTwo |
| `splice/` (6 CSV) | splice-site corrections, `CB, UB` | `bamsort_splice_sites.py` |
| `raw_cellIDs/D195_raw_cellIDs.csv` | cellranger-called barcodes, one `cellID` column | stepOne |

387 KB in total; the whole `demo/` directory, with the expected output, is 694 KB.

The `alignment/` tables are the output of `bamsort_alignment_reads.py`, which is
withheld from this release and shared on request (see
[PRIVATE_SCRIPTS.md](../PRIVATE_SCRIPTS.md)). Genomic window coordinates are in
[`bamsort_coordinates_SIVmac239.csv`](../extractionScripts/pipelineScripts/bamsort/coordinates/bamsort_coordinates_SIVmac239.csv).

Because `raw_cellIDs` are supplied precomputed, no cellranger matrix is read and
Seurat is never loaded. The demo needs only `data.table`, `dplyr` and `tidyr`.

## What it does

| step | script | what it does |
|---|---|---|
| stepThree | `bamsort/bamsort_process.R` | collapses reads to one splice category per UMI, applies the splice correction, resolves multi-mapped UMIs, filters to called cells |
| stepFour | `bamsort/stepFour_minimums_CaRe.R` | applies the per-sample detection thresholds and writes the per-cell table the figures read |

## Expected output

Written to `<work>/Figure1/CaRe/ART/bamsort/reads/`, then compared against
`expected_output/`:

| file | contents |
|---|---|
| `D195_UMI_read_counts_raw.csv` | 4,062 UMIs from 3,201 cells, unfiltered |
| `D195_UMI_read_counts_full.csv` | 2,333 UMIs from 1,745 called cells |
| `minimum/D195_CaRe_filtered.csv` | 52 cells passing `min_reads = 2`, `min_umi = 3` |
| `minimum/D195_filters.txt` | the thresholds used, for the figures to verify against |

Console output ends with:

```
OK       D195_UMI_read_counts_raw.csv
OK       D195_UMI_read_counts_full.csv
OK       minimum/D195_CaRe_filtered.csv
OK       minimum/D195_filters.txt

== demo passed
```

Three `LTR_D1` tables are empty and produce a `File is empty` warning. That is
correct: no library in this sample yielded a read in that window, and the
summary files record it as `Rows 0`.

## Expected run time

**About 5 seconds** on a normal desktop, after R and the three packages are
installed.

To install R and the dependencies, see the installation section of the
[top-level README](../README.md).

## Running it on a different sample

`run_demo.sh` stages its bundled input into the standard `extractedData/` layout
and then sets the same environment variables the SLURM wrappers set. To run the
same two steps on your own data, use
`extractionScripts/submissionScripts/Figure1/stepThree_process.sh` and
`runAll.sh stepFour` instead, which do this across all samples and groups.
