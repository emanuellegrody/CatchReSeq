#!/bin/bash
# stepTwo_bamsort.sh (host targets, EGS024)
# Reads stepTwo_bamsort.csv and submits one bamsort sbatch job per sample.
# Edit the paths below + the CSV, then run: bash stepTwo_bamsort.sh
#
# stepTwo_bamsort.csv columns:
#   sample                — SAMPLE name (matches the BAM filename stem and is
#                           used as the sample key throughout the pipeline)
#   raw_cellIDs_dataset   — 10X dataset suffix; combined with GENOME_PREFIX and
#                           raw_cellIDs_folder to locate the cellranger output
#                           (e.g. "GFP_GEX" → "<folder>/<GENOME_PREFIX>_<dataset>/outs/raw_feature_bc_matrix/")
#   raw_cellIDs_folder    — parent directory containing the cellranger output,
#                           relative to DATA_ROOT
#                           folder (must end with a trailing slash)
#   raw_cellIDs_sample    — OPTIONAL. Filename stem for the cached cellID list
#                           (<RAW_CELLIDS_OUTPUT_DIR>/<raw_cellIDs_sample>_raw_cellIDs.csv).
#                           Defaults to SAMPLE. Set to a shared biological-sample
#                           name (e.g. "GFP_invitro") so multiple CaRe-target
#                           workers reuse one cached list rather than each
#                           regenerating an identical copy.
#   process_nickname      — OPTIONAL. Pattern used by process_bamsort_host to
#                           regex-match alignment files (".*<nickname>.*_alignment_.*\\.csv$").
#                           Defaults to SAMPLE, which matches files produced by
#                           this worker; override only when aggregating files
#                           that don't carry the SAMPLE name in their filename.

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can be
# cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"
PLOT_DIR="${SCRIPT_DIR}/../../../plotScripts"

# Filesystem root for every input read and output written below. Override to run
# against a different copy of the data:
#   DATA_ROOT=/my/data bash stepTwo_bamsort.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

### UPDATE ALL PATHS (shared across samples)
# --- Step 1 (alignment) ---
export COORDS_CSV="${PIPELINE_DIR}/bamsort/coordinates/bamsort_coordinates_host_Mmul10.csv"
export INPUT_ROOT="${DATA_ROOT}/extractedData/Figure4/CaRe/counts/concat/Mmul_10_mac239_GFP_"
export INPUT_END="/outs/possorted_genome_bam.bam"
export OUTPUT_ROOT="${DATA_ROOT}/extractedData/Figure4/CaRe/stepTwo/alignment/"
export ALIGNMENT_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_alignment_reads_host.py"
# Fail at submit time rather than inside each queued job: the worker checks this
# too, but by then the error is buried in a SLURM log.
if [ ! -f "$ALIGNMENT_SCRIPT" ]; then
    echo "Error: $(basename "$ALIGNMENT_SCRIPT") is not present."
    echo
    echo "It is withheld from the public release of this repository and is"
    echo "available on request. See PRIVATE_SCRIPTS.md for what it does, its"
    echo "interface, and how to ask for it. Drop the file in"
    echo "extractionScripts/pipelineScripts/ and re-run; nothing else changes."
    exit 2
fi


# --- Step 2 (R processing: raw_cellIDs + process_bamsort) ---
export R_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_process_host.R"
export FUNCTIONS_R="${PLOT_DIR}/functions.R"
export RAW_CELLIDS_OUTPUT_DIR="${DATA_ROOT}/extractedData/Figure4/singleCell/raw_cellIDs"
export GENOME_PREFIX="Mmul_10_mac239"
export MIN_UMI_COUNT="100"
export PROCESS_INPUT_DIRS="${OUTPUT_ROOT}"
export PROCESS_OUTPUT_DIR="${DATA_ROOT}/extractedData/Figure4/CaRe/stepTwo/processed/"
# Optional |-separated allowed targets
export PROCESS_TARGETS="PTPRC|BCL2L1|IL7R"

### DO NOT UPDATE BELOW HERE
SAMPLES_CSV="${SCRIPT_DIR}/stepTwo_bamsort.csv"
WORKER="${PIPELINE_DIR}/bamsort/stepTwo_bamsort_worker_host.sh"

if [ ! -f "$SAMPLES_CSV" ]; then
    echo "Error: samples CSV not found: $SAMPLES_CSV"
    exit 1
fi
if [ ! -f "$WORKER" ]; then
    echo "Error: worker script not found: $WORKER"
    exit 1
fi
if [ ! -f "$COORDS_CSV" ]; then
    echo "Error: coords CSV not found: $COORDS_CSV"
    exit 1
fi

mkdir -p logs

tail -n +2 "$SAMPLES_CSV" | while IFS=, read -r SAMPLE RAW_CELLIDS_DATASET RAW_CELLIDS_FOLDER RAW_CELLIDS_SAMPLE PROCESS_NICKNAME; do
    SAMPLE="$(echo "$SAMPLE" | tr -d '[:space:]\r')"
    RAW_CELLIDS_DATASET="$(echo "${RAW_CELLIDS_DATASET:-}" | tr -d '[:space:]\r')"
    RAW_CELLIDS_FOLDER="$(echo "${RAW_CELLIDS_FOLDER:-}" | tr -d '[:space:]\r')"
    # The CSV path is relative to DATA_ROOT; an absolute path is taken as given.
    case "$RAW_CELLIDS_FOLDER" in /*) ;; *) RAW_CELLIDS_FOLDER="${DATA_ROOT}/$RAW_CELLIDS_FOLDER" ;; esac
    RAW_CELLIDS_SAMPLE="$(echo "${RAW_CELLIDS_SAMPLE:-}" | tr -d '[:space:]\r')"
    PROCESS_NICKNAME="$(echo "${PROCESS_NICKNAME:-}" | tr -d '[:space:]\r')"
    [ -z "$SAMPLE" ] && continue

    if [ -z "$RAW_CELLIDS_DATASET" ] || [ -z "$RAW_CELLIDS_FOLDER" ]; then
        echo "Error: row for sample '$SAMPLE' missing raw_cellIDs_dataset and/or raw_cellIDs_folder"
        echo "       expected CSV header: sample,raw_cellIDs_dataset,raw_cellIDs_folder[,raw_cellIDs_sample][,process_nickname]"
        exit 1
    fi

    export SAMPLE RAW_CELLIDS_DATASET RAW_CELLIDS_FOLDER RAW_CELLIDS_SAMPLE PROCESS_NICKNAME

    echo "Submitting bamsort: $SAMPLE"
    sbatch \
        --job-name="bamsort-${SAMPLE}" \
        --output="logs/stepTwo_bamsort_${SAMPLE}_%j.txt" \
        --error="logs/stepTwo_bamsort_${SAMPLE}_%j.err" \
        --export=ALL \
        "$WORKER"
done

echo "All bamsort jobs submitted."
