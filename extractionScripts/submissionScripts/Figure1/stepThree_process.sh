#!/bin/bash
# stepThree_process.sh  (Figure1 Ca*Re)
# Combines the stepTwo per-region alignment tables into one category-per-UMI
# table per sample: <group>/bamsort/reads/<sample>_UMI_read_counts_full.csv
# (filtered to raw cellIDs) and _UMI_read_counts_raw.csv (unfiltered).
#
# Reads the same stepTwo_bamsort.csv as stepTwo_bamsort.sh and submits one job
# per group, each covering that group's samples and reading/writing only inside
# its own <CARE_ROOT>/<group>/ tree. Edit the paths below, then run:
#   bash stepThree_process.sh
#
# Run after all stepTwo_bamsort.sh jobs of a group have finished
#
# raw_cellIDs are derived from the single-cell TenX runs, not from the Ca*Re
# libraries, so they stay in one shared dir across all groups.
#

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
CARE_ROOT="${DATA_ROOT}/extractedData/Figure1/CaRe"
export R_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_process.R"
export FUNCTIONS_R="${PIPELINE_DIR}/bamsort/bamsort_functions.R"
export RAW_CELLIDS_OUTPUT_DIR="${DATA_ROOT}/extractedData/Figure1/singleCell/raw_cellIDs"
export GENOME_PREFIX="Mmul_10_mac239"
export MIN_UMI_COUNT="100"

### DO NOT UPDATE BELOW HERE
SAMPLES_CSV="${SCRIPT_DIR}/stepTwo_bamsort.csv"
WORKER="${PIPELINE_DIR}/bamsort/stepThree_process_worker.sh"

for f in "$SAMPLES_CSV" "$WORKER"; do
    [ -f "$f" ] || { echo "Error: file not found: $f"; exit 1; }
done

mkdir -p logs

# Groups in the order they first appear in the CSV
GROUP_LIST=$(tail -n +2 "$SAMPLES_CSV" | cut -d, -f2 | tr -d '[:blank:]\r' | awk 'NF && !seen[$0]++')
if [ -z "$GROUP_LIST" ]; then
    echo "Error: no groups read from $SAMPLES_CSV"
    exit 1
fi

for GROUP in $GROUP_LIST; do
    # Build the |-separated parallel lists the R script expects, one entry per
    # sample row of this group.
    SAMPLES=""; DATASETS=""; FOLDERS=""
    while IFS=, read -r SAMPLE ROW_GROUP _ALIGNMENT_BAM_IDS RAW_CELLIDS_DATASET RAW_CELLIDS_FOLDER; do
        SAMPLE="$(echo "$SAMPLE" | tr -d '[:space:]\r')"
        ROW_GROUP="$(echo "$ROW_GROUP" | tr -d '[:space:]\r')"
        [ -z "$SAMPLE" ] && continue
        [ "$ROW_GROUP" = "$GROUP" ] || continue
        RAW_CELLIDS_DATASET="$(echo "$RAW_CELLIDS_DATASET" | tr -d '[:space:]\r')"
        RAW_CELLIDS_FOLDER="$(echo "$RAW_CELLIDS_FOLDER" | tr -d '[:space:]\r')"

        SAMPLES="${SAMPLES:+$SAMPLES|}${SAMPLE}"
        DATASETS="${DATASETS:+$DATASETS|}${RAW_CELLIDS_DATASET}"
        FOLDERS="${FOLDERS:+$FOLDERS|}${RAW_CELLIDS_FOLDER}"
    done < <(tail -n +2 "$SAMPLES_CSV")

    if [ -z "$SAMPLES" ]; then
        echo "Error: no samples read for group $GROUP"
        exit 1
    fi

    GROUP_DIR="${CARE_ROOT}/${GROUP}"
    export PROCESS_INPUT_DIRS="${GROUP_DIR}/bamsort/alignment/"
    export PROCESS_OUTPUT_DIR="${GROUP_DIR}/bamsort/reads/"
    export PROCESS_CORRECTION_DIR="${GROUP_DIR}/bamsort/splice/combined/"
    export PROCESS_SAMPLES="$SAMPLES"
    export RAW_CELLIDS_SAMPLES="$SAMPLES"
    export RAW_CELLIDS_DATASETS="$DATASETS"
    export RAW_CELLIDS_FOLDERS="$FOLDERS"

    echo "Submitting stepThree for group ${GROUP}: ${SAMPLES//|/, }"
    sbatch \
        --job-name="stepThree-${GROUP}" \
        --output="logs/stepThree_process_${GROUP}_%j.txt" \
        --error="logs/stepThree_process_${GROUP}_%j.err" \
        --export=ALL \
        "$WORKER"
done

echo "All stepThree jobs submitted."
