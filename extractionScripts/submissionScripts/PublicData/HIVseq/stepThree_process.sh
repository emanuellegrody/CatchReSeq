#!/bin/bash
# stepThree_process.sh  (HIVseq)
# Combines the stepTwo per-region alignment tables into one category-per-UMI
# table per sample: <PROCESS_OUTPUT_DIR>/<sample>_UMI_read_counts_full.csv
# (filtered to raw cellIDs) and _UMI_read_counts_raw.csv (unfiltered).
#
# Pairs with stepTwo_bamsort.sh.
#
# Reads samples_manifest.csv and submits ONE job covering all samples. Run AFTER
# stepTwo_bamsort.sh has finished, since this step gathers across every region
# of every sample. Edit the paths below, then run:
#   bash stepThree_process.sh
#
# Detection thresholds are NOT applied here. Call set_minimums() on the
# _UMI_read_counts_full.csv from the analysis scripts so thresholds can be
# varied without re-running this step — which is what the public datasets need,
# since their sequencing depth differs from the CaRe runs.

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can
# be cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../../pipelineScripts"

# NOTE: not run for this dataset. Figure3.Rmd reads the stepTwo alignment
# and splice tables directly (hxb2_per_cell), so no per-UMI category table
# is built here; these paths follow the convention the rest of the dataset
# uses rather than pointing at existing output.
### UPDATE ALL PATHS
export R_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_process.R"
export FUNCTIONS_R="${PIPELINE_DIR}/bamsort/bamsort_functions.R"
# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data bash stepThree_process.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

export RAW_CELLIDS_OUTPUT_DIR="${DATA_ROOT}/extractedData/HIVseq/raw_cellIDs"
export GENOME_PREFIX="GRCh38_HXB2"
export MIN_UMI_COUNT="100"
# Parent dir of the cellranger runs, i.e. <folder>/<GENOME_PREFIX>_<sample>/outs/
COUNTS_ROOT="${DATA_ROOT}/publicData/HIVseq/counts/"
export PROCESS_INPUT_DIRS="${DATA_ROOT}/extractedData/HIVseq/bamsort/"
export PROCESS_OUTPUT_DIR="${DATA_ROOT}/extractedData/HIVseq/bamsort/reads/"
# stepTwo_bamsort.sh extracts splice junctions but does not run the
# bamsort_splice_sites.py filter, so there are no *_splicesites_filtered.csv to
# apply. Set this once that filter is run for these samples.
# export PROCESS_CORRECTION_DIR=".../bamsort/splice/combined/"

### DO NOT UPDATE BELOW HERE
SAMPLES_CSV="${SCRIPT_DIR}/samples_manifest.csv"
WORKER="${PIPELINE_DIR}/bamsort/stepThree_process_worker.sh"

for f in "$SAMPLES_CSV" "$WORKER"; do
    [ -f "$f" ] || { echo "Error: file not found: $f"; exit 1; }
done

# One entry per manifest row. The cellranger run dir is named after the sample
# itself here, so the dataset and folder lists are the sample name and a
# constant.
SAMPLES=""; DATASETS=""; FOLDERS=""
while IFS=, read -r SAMPLE _rest; do
    SAMPLE="$(echo "$SAMPLE" | tr -d '[:space:]\r')"
    [ -z "$SAMPLE" ] && continue

    SAMPLES="${SAMPLES:+$SAMPLES|}${SAMPLE}"
    DATASETS="${DATASETS:+$DATASETS|}${SAMPLE}"
    FOLDERS="${FOLDERS:+$FOLDERS|}${COUNTS_ROOT}"
done < <(tail -n +2 "$SAMPLES_CSV")

if [ -z "$SAMPLES" ]; then
    echo "Error: no samples read from $SAMPLES_CSV"
    exit 1
fi

export PROCESS_SAMPLES="$SAMPLES"
export RAW_CELLIDS_SAMPLES="$SAMPLES"
export RAW_CELLIDS_DATASETS="$DATASETS"
export RAW_CELLIDS_FOLDERS="$FOLDERS"

mkdir -p logs

echo "Submitting stepThree for: ${SAMPLES//|/, }"
sbatch \
    --job-name="stepThree-HIVseq" \
    --output="logs/stepThree_process_%j.txt" \
    --error="logs/stepThree_process_%j.err" \
    --export=ALL \
    "$WORKER"
