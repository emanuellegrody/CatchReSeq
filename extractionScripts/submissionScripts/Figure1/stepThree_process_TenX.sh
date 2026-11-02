#!/bin/bash
# stepThree_process_TenX.sh  (Figure1 TenX)
# Combines the stepTwo_bamsort_TenX.sh per-region alignment tables into one
# category-per-UMI table per sample:
#   <TenX_ROOT>/bamsort/reads/<sample>_UMI_read_counts_full.csv  (raw_cellID filtered)
#   <TenX_ROOT>/bamsort/reads/<sample>_UMI_read_counts_raw.csv   (unfiltered)
#
# Same worker and same R script as the Ca*Re stepThree_process.sh.
#
# Run AFTER every stepTwo_bamsort_TenX.sh job has finished. Follow with 
# stepFour_minimums_TenX.R

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can be
# cloned anywhere
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

# Filesystem root for every input read and output written below. Override to run
# the pipeline against a different copy of the data:
#   DATA_ROOT=/my/data bash <this script>
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

### UPDATE ALL PATHS
TenX_ROOT="${DATA_ROOT}/extractedData/Figure1/singleCell"
export R_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_process.R"
export FUNCTIONS_R="${PIPELINE_DIR}/bamsort/bamsort_functions.R"
export RAW_CELLIDS_OUTPUT_DIR="${TenX_ROOT}/raw_cellIDs"
export MIN_UMI_COUNT="100"

# Every sample already has a cached <sample>_raw_cellIDs.csv in
# RAW_CELLIDS_OUTPUT_DIR, so bamsort_process.R loads the barcode lists from disk.
export GENOME_PREFIX="Mmul_10_mac239"

export PROCESS_INPUT_DIRS="${TenX_ROOT}/bamsort/alignment/"
export PROCESS_OUTPUT_DIR="${TenX_ROOT}/bamsort/reads/"
export PROCESS_CORRECTION_DIR="${TenX_ROOT}/bamsort/splice/combined/"

### DO NOT UPDATE BELOW HERE
SAMPLES_CSV="${SCRIPT_DIR}/stepTwo_bamsort_TenX.csv"
WORKER="${PIPELINE_DIR}/bamsort/stepThree_process_worker.sh"

for f in "$SAMPLES_CSV" "$WORKER"; do
    [ -f "$f" ] || { echo "Error: file not found: $f"; exit 1; }
done

mkdir -p logs

# Build the |-separated parallel lists the R script expects, one entry per row.
SAMPLES=""; DATASETS=""; FOLDERS=""
while IFS=, read -r SAMPLE _COUNTS_DIR RAW_CELLIDS_DATASET RAW_CELLIDS_FOLDER; do
    SAMPLE="$(echo "$SAMPLE" | tr -d '[:space:]\r')"
    [ -z "$SAMPLE" ] && continue
    RAW_CELLIDS_DATASET="$(echo "$RAW_CELLIDS_DATASET" | tr -d '[:space:]\r')"
    RAW_CELLIDS_FOLDER="$(echo "$RAW_CELLIDS_FOLDER" | tr -d '[:space:]\r')"

    SAMPLES="${SAMPLES:+$SAMPLES|}${SAMPLE}"
    DATASETS="${DATASETS:+$DATASETS|}${RAW_CELLIDS_DATASET}"
    FOLDERS="${FOLDERS:+$FOLDERS|}${RAW_CELLIDS_FOLDER}"
done < <(tail -n +2 "$SAMPLES_CSV")

if [ -z "$SAMPLES" ]; then
    echo "Error: no samples read from $SAMPLES_CSV"
    exit 1
fi

export PROCESS_SAMPLES="$SAMPLES"
export RAW_CELLIDS_SAMPLES="$SAMPLES"
export RAW_CELLIDS_DATASETS="$DATASETS"
export RAW_CELLIDS_FOLDERS="$FOLDERS"

echo "Submitting TenX stepThree: ${SAMPLES//|/, }"
sbatch \
    --job-name="tenxThree" \
    --output="logs/stepThree_process_TenX_%j.txt" \
    --error="logs/stepThree_process_TenX_%j.err" \
    --export=ALL \
    "$WORKER"

echo "TenX stepThree job submitted."
