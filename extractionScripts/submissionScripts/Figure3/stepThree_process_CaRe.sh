#!/bin/bash
# stepThree_process_CaRe.sh  (EGS026 CaRe HIV clade amplicons, Figure3)
# Combines the stepTwo per-region alignment tables into one category-per-UMI
# table per clade: <PROCESS_OUTPUT_DIR>/<sample>_UMI_read_counts_full.csv
# (filtered to raw cellIDs) and _UMI_read_counts_raw.csv (unfiltered).
#
# Named apart from stepThree_process.sh in this directory, which processes the
# single-cell HXB2 bamsort of the 10X pools; this one is the CaRe amplicon
# pipeline and pairs with stepTwo_bamsort_CaRe.sh.
#
# Reads the same stepTwo_bamsort_CaRe.csv and submits ONE job covering all
# clades. Run AFTER all stepTwo_bamsort_CaRe.sh jobs have finished, since this
# step gathers across every region of every clade.
#
# The sample token carries a trailing underscore ("A_", not "A"). process_bamsort
# selects a sample's alignment CSVs by substring match on the filename, and a bare
# "A" also matches the A1/A7 coordinate targets inside other clades' files
# (C_env_bamsort_alignment_A1_D4.csv); "A_" matches only clade A's libraries.
# The same token keys raw_cellIDs and the splice-correction filename.
#
# Detection thresholds are NOT applied here. Call set_minimums() on the
# _UMI_read_counts_full.csv from the analysis scripts so thresholds can be
# varied without re-running this step.

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can
# be cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

### UPDATE ALL PATHS
export R_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_process.R"
export FUNCTIONS_R="${PIPELINE_DIR}/bamsort/bamsort_functions.R"
# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data bash stepThree_process_CaRe.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

export RAW_CELLIDS_OUTPUT_DIR="${DATA_ROOT}/extractedData/Figure3/singleCell/raw_cellIDs"
export GENOME_PREFIX="GRCh38_HXB2"
export MIN_UMI_COUNT="100"
export PROCESS_INPUT_DIRS="${DATA_ROOT}/extractedData/Figure3/CaRe/bamsort/alignment/"
export PROCESS_OUTPUT_DIR="${DATA_ROOT}/extractedData/Figure3/CaRe/bamsort/reads/"
export PROCESS_CORRECTION_DIR="${DATA_ROOT}/extractedData/Figure3/CaRe/bamsort/splice/combined/"

### DO NOT UPDATE BELOW HERE
SAMPLES_CSV="${SCRIPT_DIR}/stepTwo_bamsort_CaRe.csv"
WORKER="${PIPELINE_DIR}/bamsort/stepThree_process_worker.sh"

for f in "$SAMPLES_CSV" "$WORKER"; do
    [ -f "$f" ] || { echo "Error: file not found: $f"; exit 1; }
done

# Build the |-separated parallel lists the R script expects, one entry per
# sample row of the CSV.
SAMPLES=""; DATASETS=""; FOLDERS=""
while IFS=, read -r SAMPLE _ALIGNMENT_BAM_IDS RAW_CELLIDS_DATASET RAW_CELLIDS_FOLDER; do
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

mkdir -p logs

echo "Submitting stepThree for: ${SAMPLES//|/, }"
sbatch \
    --job-name="stepThree-Figure3-CaRe" \
    --output="logs/stepThree_process_CaRe_%j.txt" \
    --error="logs/stepThree_process_CaRe_%j.err" \
    --export=ALL \
    "$WORKER"
