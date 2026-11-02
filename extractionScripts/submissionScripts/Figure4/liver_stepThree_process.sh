#!/bin/bash
# liver_stepThree_process.sh  (host targets, Figure4 liver section)
# Combines the stepTwo per-region alignment tables into one target-per-UMI
# table per sample: <PROCESS_OUTPUT_DIR>/<sample>_UMI_read_counts_full.csv
# (filtered to raw cellIDs) and _UMI_read_counts_raw.csv (unfiltered).
#
# Host equivalent of the viral stepThree: no splice-category resolution, so
# each UMI keeps the target its alignment file came from.
#
# Reads the same liver_stepTwo_bamsort.csv as liver_stepTwo_bamsort.sh and submits ONE job
# covering all samples. Run AFTER all liver_stepTwo_bamsort.sh jobs have finished,
# since this step gathers across every target of every sample.
#
# Detection thresholds are NOT applied here. Call set_minimums(mode = "host")
# on the _UMI_read_counts_full.csv from the analysis scripts so thresholds can
# be varied without re-running this step.

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can
# be cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

# Filesystem roots for every input read and output written below. Override to
# run against a different copy of the data:
#   DATA_ROOT=/my/data bash liver_stepThree_process.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

### UPDATE ALL PATHS
export R_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_process_host.R"
export FUNCTIONS_R="${PIPELINE_DIR}/bamsort/bamsort_functions.R"
export RAW_CELLIDS_OUTPUT_DIR="${DATA_ROOT}/extractedData/Figure4/liver/singleCell/raw_cellIDs"
export GENOME_PREFIX="GRCm39"
export MIN_UMI_COUNT="100"
export PROCESS_INPUT_DIRS="${DATA_ROOT}/extractedData/Figure4/liver/CaRe/stepTwo/alignment/"
export PROCESS_OUTPUT_DIR="${DATA_ROOT}/extractedData/Figure4/liver/CaRe/stepTwo/processed/"
# Optional |-separated allowed targets; unset accepts every target found
export PROCESS_TARGETS="Slk|Dgkd|Lsm14b|Kras|Nf2"

### DO NOT UPDATE BELOW HERE
SAMPLES_CSV="${SCRIPT_DIR}/liver_stepTwo_bamsort.csv"
WORKER="${PIPELINE_DIR}/bamsort/stepThree_process_worker.sh"

for f in "$SAMPLES_CSV" "$WORKER"; do
    [ -f "$f" ] || { echo "Error: file not found: $f"; exit 1; }
done

# Build the |-separated parallel lists the R script expects, one entry per
# sample row of the CSV. raw_cellIDs_sample repeats across rows that share a
# biological sample (e.g. every WT_* target maps to "Normal"); the repeats are
# kept so the lists stay row-aligned, and bamsort_process_host.R loads the
# cached cellID list after generating it once.
SAMPLES=""; DATASETS=""; FOLDERS=""; CELLID_SAMPLES=""; NICKNAMES=""
while IFS=, read -r SAMPLE RAW_CELLIDS_DATASET RAW_CELLIDS_FOLDER RAW_CELLIDS_SAMPLE PROCESS_NICKNAME; do
    SAMPLE="$(echo "$SAMPLE" | tr -d '[:space:]\r')"
    [ -z "$SAMPLE" ] && continue
    RAW_CELLIDS_DATASET="$(echo "${RAW_CELLIDS_DATASET:-}" | tr -d '[:space:]\r')"
    RAW_CELLIDS_FOLDER="$(echo "${RAW_CELLIDS_FOLDER:-}" | tr -d '[:space:]\r')"
    RAW_CELLIDS_SAMPLE="$(echo "${RAW_CELLIDS_SAMPLE:-}" | tr -d '[:space:]\r')"
    PROCESS_NICKNAME="$(echo "${PROCESS_NICKNAME:-}" | tr -d '[:space:]\r')"

    # The CSV path is relative to DATA_ROOT; an absolute path is taken as given.
    # bamsort_process_host.R pastes this straight into the Read10X path, and the
    # job's working directory is not DATA_ROOT, so it has to be absolute here.
    case "$RAW_CELLIDS_FOLDER" in /*) ;; *) RAW_CELLIDS_FOLDER="${DATA_ROOT}/$RAW_CELLIDS_FOLDER" ;; esac

    if [ -z "$RAW_CELLIDS_DATASET" ] || [ -z "$RAW_CELLIDS_FOLDER" ]; then
        echo "Error: row for sample '$SAMPLE' missing raw_cellIDs_dataset and/or raw_cellIDs_folder"
        echo "       expected CSV header: sample,raw_cellIDs_dataset,raw_cellIDs_folder[,raw_cellIDs_sample][,process_nickname]"
        exit 1
    fi

    SAMPLES="${SAMPLES:+$SAMPLES|}${SAMPLE}"
    DATASETS="${DATASETS:+$DATASETS|}${RAW_CELLIDS_DATASET}"
    FOLDERS="${FOLDERS:+$FOLDERS|}${RAW_CELLIDS_FOLDER}"
    CELLID_SAMPLES="${CELLID_SAMPLES:+$CELLID_SAMPLES|}${RAW_CELLIDS_SAMPLE:-$SAMPLE}"
    NICKNAMES="${NICKNAMES:+$NICKNAMES|}${PROCESS_NICKNAME:-$SAMPLE}"
done < <(tail -n +2 "$SAMPLES_CSV")

if [ -z "$SAMPLES" ]; then
    echo "Error: no samples read from $SAMPLES_CSV"
    exit 1
fi

export PROCESS_SAMPLES="$SAMPLES"
export PROCESS_NICKNAMES="$NICKNAMES"
export RAW_CELLIDS_SAMPLES="$CELLID_SAMPLES"
export RAW_CELLIDS_DATASETS="$DATASETS"
export RAW_CELLIDS_FOLDERS="$FOLDERS"

mkdir -p logs

echo "Submitting stepThree for: ${SAMPLES//|/, }"
sbatch \
    --job-name="stepThree-liver" \
    --output="logs/stepThree_process_%j.txt" \
    --error="logs/stepThree_process_%j.err" \
    --export=ALL \
    "$WORKER"
