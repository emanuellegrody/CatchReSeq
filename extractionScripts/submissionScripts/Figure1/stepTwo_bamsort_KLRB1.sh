#!/bin/bash
# stepTwo_bamsort_KLRB1.sh  (Figure1 KLRB1 panels A-D)
# Reads stepTwo_bamsort_KLRB1.csv and submits bamsort sbatch job
# Run after stepOne_mkcounts_KLRB1.sh for the Ca*Re rows
#
# There is no stepThree or stepFour here. Figure1.Rmd collapses reads per UMI
# and applies the detection threshold in the chunks under "# KLRB1".
#
# stepTwo_bamsort_KLRB1.csv columns (paths are relative to DATA_ROOT):
#   sample         — label for the job name and logs
#   output_prefix  — filename stem, giving <prefix>_bamsort_alignment_KLRB1.csv
#   bam            — the cellranger BAM to read
#   output_root    — directory to write into (must end with a trailing slash)
#
# Edit the paths below + the CSV, then run: bash stepTwo_bamsort_KLRB1.sh

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can be
# cloned anywhere
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

# Filesystem root for every input read and output written below. Override to run
# the pipeline against a different copy of the data:
#   DATA_ROOT=/my/data bash stepTwo_bamsort_KLRB1.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

### UPDATE ALL PATHS (shared across samples)
export COORDS_CSV="${PIPELINE_DIR}/bamsort/coordinates/bamsort_coordinates_KLRB1.csv"
export ALIGNMENT_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_alignment_reads.py"

### DO NOT UPDATE BELOW HERE
SAMPLES_CSV="${SCRIPT_DIR}/stepTwo_bamsort_KLRB1.csv"
WORKER="${PIPELINE_DIR}/bamsort/stepTwo_bamsort_worker_KLRB1.sh"

for f in "$SAMPLES_CSV" "$WORKER" "$COORDS_CSV" "$ALIGNMENT_SCRIPT"; do
    if [ ! -f "$f" ]; then
        case "$(basename "$f")" in
            bamsort_alignment_reads.py|bamsort_alignment_reads_host.py)
                echo "Error: $(basename "$f") is not present."
                echo
                echo "It is withheld from the public release of this repository and is"
                echo "available on request. See PRIVATE_SCRIPTS.md for what it does, its"
                echo "interface, and how to request access to it. Drop the file in"
                echo "extractionScripts/pipelineScripts/ and re-run; nothing else changes."
                exit 2
                ;;
        esac
        echo "Error: required file not found: $f"
        exit 1
    fi
done

mkdir -p logs

tail -n +2 "$SAMPLES_CSV" | while IFS=, read -r SAMPLE OUTPUT_PREFIX BAM_REL OUTPUT_REL; do
    SAMPLE="$(echo "$SAMPLE" | tr -d '[:space:]\r')"
    OUTPUT_PREFIX="$(echo "${OUTPUT_PREFIX:-}" | tr -d '[:space:]\r')"
    BAM_REL="$(echo "${BAM_REL:-}" | tr -d '[:space:]\r')"
    OUTPUT_REL="$(echo "${OUTPUT_REL:-}" | tr -d '[:space:]\r')"
    [ -z "$SAMPLE" ] && continue

    if [ -z "$OUTPUT_PREFIX" ] || [ -z "$BAM_REL" ] || [ -z "$OUTPUT_REL" ]; then
        echo "Error: row for sample '$SAMPLE' is missing a field"
        echo "       expected CSV header: sample,output_prefix,bam,output_root"
        exit 1
    fi

    BAM_FILE="${DATA_ROOT}/${BAM_REL}"
    OUTPUT_ROOT="${DATA_ROOT}/${OUTPUT_REL}"

    if [ ! -f "$BAM_FILE" ]; then
        echo "Error: BAM not found for '$SAMPLE': $BAM_FILE"
        echo "       Ca*Re rows need stepOne_mkcounts_KLRB1.sh to have run first."
        exit 1
    fi

    export SAMPLE OUTPUT_PREFIX BAM_FILE OUTPUT_ROOT

    echo "Submitting bamsort: $SAMPLE"
    sbatch \
        --job-name="bamsortK-${SAMPLE}" \
        --output="logs/stepTwo_bamsort_KLRB1_${SAMPLE}_%j.txt" \
        --error="logs/stepTwo_bamsort_KLRB1_${SAMPLE}_%j.err" \
        --export=ALL \
        "$WORKER"
done

echo "All KLRB1 bamsort jobs submitted."
