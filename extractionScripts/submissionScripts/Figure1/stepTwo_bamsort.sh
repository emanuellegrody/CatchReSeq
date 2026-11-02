#!/bin/bash
# stepTwo_bamsort.sh  (Figure1 Ca*Re)
# Reads stepTwo_bamsort.csv and submits one bamsort sbatch job per sample:
# per-region alignment tables + splice junctions, over all five amplicon
# libraries of that sample (env, SSenv, pol, D1, D4).
# Edit the paths below + the CSV, then run: bash stepTwo_bamsort.sh
#
# Run AFTER stepOne_mkcounts.sh: both the alignment BAMs and the splice BAMs are
# the cellranger runs it produces under <CARE_ROOT>/<group>/counts/.
#
# Follow with stepThree_process.sh, which combines the per-region tables into
# one category-per-UMI table per sample.
#
# Every path is derived per row from the CSV group column, so each group writes
# only inside its own <CARE_ROOT>/<group>/ tree.

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can be
# cloned anywhere
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

# Filesystem root for every input read and output written below. Override to run
# the pipeline against a different copy of the data:
#   DATA_ROOT=/my/data bash <this script>
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

### UPDATE ALL PATHS (shared across samples)
CARE_ROOT="${DATA_ROOT}/extractedData/Figure1/CaRe"

# --- Step 1 (alignment) ---
export COORDS_CSV="${PIPELINE_DIR}/bamsort/coordinates/bamsort_coordinates_SIVmac239.csv"
export INPUT_END="/outs/possorted_genome_bam.bam"
export ALIGNMENT_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_alignment_reads.py"
if [ ! -f "$ALIGNMENT_SCRIPT" ]; then
    echo "Error: $(basename "$ALIGNMENT_SCRIPT") is not present."
    echo
    echo "It is withheld from the public release of this repository and is"
    echo "available on request. See PRIVATE_SCRIPTS.md for what it does, its"
    echo "interface, and how to request access to it. Drop the file in"
    echo "extractionScripts/pipelineScripts/ and re-run; nothing else changes."
    exit 2
fi

GENOME_PREFIX="Mmul_10_mac239"

# --- Step 2 (splice extraction + filter) ---
export SPLICE_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_splice.py"
export SPLICE_FILTER_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_splice_sites.py"
export SPLICE_INPUT_END="/outs/possorted_genome_bam.bam"
export SPLICE_TARGETS="_D1_PL|_env_PL|_pol_PL|_SSenv_PL|_tat"
export SPLICE_FILTER_DONOR="D4"
export SPLICE_MIN_INTRON="500"

### DO NOT UPDATE BELOW HERE
SAMPLES_CSV="${SCRIPT_DIR}/stepTwo_bamsort.csv"
WORKER="${PIPELINE_DIR}/bamsort/stepTwo_bamsort_worker.sh"

if [ ! -f "$SAMPLES_CSV" ]; then
    echo "Error: samples CSV not found: $SAMPLES_CSV"
    exit 1
fi
if [ ! -f "$WORKER" ]; then
    echo "Error: worker script not found: $WORKER"
    exit 1
fi

mkdir -p logs

tail -n +2 "$SAMPLES_CSV" | while IFS=, read -r SAMPLE GROUP ALIGNMENT_BAM_IDS RAW_CELLIDS_DATASET RAW_CELLIDS_FOLDER; do
    SAMPLE="$(echo "$SAMPLE" | tr -d '[:space:]')"
    [ -z "$SAMPLE" ] && continue
    GROUP="$(echo "$GROUP" | tr -d '[:space:]')"
    RAW_CELLIDS_FOLDER="$(echo "${RAW_CELLIDS_FOLDER:-}" | tr -d '[:space:]\r')"
    # The CSV path is relative to DATA_ROOT; an absolute path is taken as given.
    case "$RAW_CELLIDS_FOLDER" in /*) ;; *) RAW_CELLIDS_FOLDER="${DATA_ROOT}/$RAW_CELLIDS_FOLDER" ;; esac
    if [ -z "$GROUP" ]; then
        echo "Error: empty group for sample $SAMPLE"
        exit 1
    fi

    GROUP_DIR="${CARE_ROOT}/${GROUP}"
    export INPUT_ROOT="${GROUP_DIR}/counts/${GENOME_PREFIX}_"
    export OUTPUT_ROOT="${GROUP_DIR}/bamsort/alignment/"
    export SPLICE_INPUT_ROOT="${GROUP_DIR}/counts/${GENOME_PREFIX}_"
    export SPLICE_OUTPUT_DIR="${GROUP_DIR}/bamsort/splice"
    # outputDir used for stepThree loading
    export PROCESS_CORRECTION_DIR="${GROUP_DIR}/bamsort/splice/combined/"

    ALIGNMENT_BAM_ID="${ALIGNMENT_BAM_IDS%%|*}"

    export SAMPLE ALIGNMENT_BAM_ID ALIGNMENT_BAM_IDS RAW_CELLIDS_DATASET RAW_CELLIDS_FOLDER

    echo "Submitting bamsort: $SAMPLE  group=${GROUP}  (${ALIGNMENT_BAM_IDS//|/, })"
    sbatch \
        --job-name="bamsort-${SAMPLE}" \
        --output="logs/stepTwo_bamsort_${SAMPLE}_%j.txt" \
        --error="logs/stepTwo_bamsort_${SAMPLE}_%j.err" \
        --export=ALL \
        "$WORKER"
done

echo "All bamsort jobs submitted."
