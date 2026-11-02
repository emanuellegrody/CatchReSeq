#!/bin/bash
# stepThree_process_SIV_TenX.sh  (Figure4 / EGS024 SIV single-cell TenX library)
# Same worker and same R script as stepThree_process_SIV.sh; only the
# input/output dirs and the sample name differ, so the two arms are categorised
# by identical code.
#
# SAMPLE is GEX_GFP_invitro, not GFP_invitro. bamsort/alignment/ holds four
# generations of this library under four prefixes:
#   GEX_GFP_invitro_*        Figure4/singleCell/bamsort   <- current
#   GFP_invitro_*            EGS024/invitro/bamsort
#   GFP_GEX_*                EGS024/TenX/bamsort
#   GFP_invitro_GFP_GEX_*    an earlier pass
# process_bamsort globs `.*<sample>.*`, so a GFP_invitro run would sweep three
# of the four together and double-count. Naming the sample GEX_GFP_invitro
# selects the current set alone, at the cost of that prefix appearing in the
# output filenames -- which is what the Figure4.Rmd SIV chunk reads.
#
# PROCESS_CORRECTION_DIR is what makes MS calls possible for TenX at all. Without
# it, a TenX library can only ever produce US, US-SS and any: extract_category
# reads the isoform from the target primer in the filename, which a TenX library
# does not have, so A1_D4 and A7_LTR fall to any-D4/any-A7. The splice filter
# promotes exactly those two to MS when the junction is confirmed at donor D4.
# The TenX worker writes one <sample>_splicesites_filtered.csv per sample, so no
# combine step is needed here.
#
# Run AFTER every stepTwo job has finished, and BEFORE stepFour minimums.

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can
# be cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

### UPDATE ALL PATHS
# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data bash stepThree_process_SIV_TenX.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

TenX_ROOT="${DATA_ROOT}/extractedData/Figure4/singleCell"
export R_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_process.R"
export FUNCTIONS_R="${PIPELINE_DIR}/bamsort/bamsort_functions.R"
export RAW_CELLIDS_OUTPUT_DIR="${TenX_ROOT}/raw_cellIDs"
export GENOME_PREFIX="Mmul_10_mac239_GFP"
export MIN_UMI_COUNT="100"

SAMPLE="GEX_GFP_invitro"
RAW_CELLIDS_DATASET="GEX"
RAW_CELLIDS_FOLDER="${TenX_ROOT}/counts/concat/"

export PROCESS_INPUT_DIRS="${TenX_ROOT}/bamsort/alignment/"
export PROCESS_OUTPUT_DIR="${TenX_ROOT}/bamsort/reads/"
export PROCESS_CORRECTION_DIR="${TenX_ROOT}/bamsort/splice/combined/"

### DO NOT UPDATE BELOW HERE
WORKER="${PIPELINE_DIR}/bamsort/stepThree_process_worker.sh"

[ -f "$WORKER" ]   || { echo "Error: file not found: $WORKER"; exit 1; }
[ -f "$R_SCRIPT" ] || { echo "Error: file not found: $R_SCRIPT"; exit 1; }
[ -d "$PROCESS_INPUT_DIRS" ] || { echo "Error: no alignment dir: $PROCESS_INPUT_DIRS"; exit 1; }

# Preflight: see stepThree_process_SIV.sh. The cache is keyed by SAMPLE, and
# both libraries were called on the one Mmul_10_mac239_GFP_GEX cellranger run,
# so the TenX_-prefixed name is seeded from the GFP_invitro cache rather than
# re-derived from the 10X matrix.
CACHE="${RAW_CELLIDS_OUTPUT_DIR}/${SAMPLE}_raw_cellIDs.csv"
BASE_CACHE="${RAW_CELLIDS_OUTPUT_DIR}/GFP_invitro_raw_cellIDs.csv"
if [ ! -f "$CACHE" ]; then
    if [ -f "$BASE_CACHE" ]; then
        cp "$BASE_CACHE" "$CACHE"
        echo "Seeded $CACHE from $BASE_CACHE ($(($(wc -l < "$CACHE") - 1)) barcodes)"
    else
        echo "Error: no raw cellID cache at $CACHE and none to seed from at $BASE_CACHE"
        exit 1
    fi
fi

N_ALIGN=$(ls "${PROCESS_INPUT_DIRS}${SAMPLE}"_bamsort_alignment_*.csv 2>/dev/null | wc -l)
[ "$N_ALIGN" -gt 0 ] || { echo "Error: no alignment CSVs matching $SAMPLE in $PROCESS_INPUT_DIRS"; exit 1; }

CORRECTION_FILE="${PROCESS_CORRECTION_DIR%/}/${SAMPLE}_splicesites_filtered.csv"
if [ -f "$CORRECTION_FILE" ]; then
    echo "Splice corrections: $CORRECTION_FILE ($(($(wc -l < "$CORRECTION_FILE") - 1)) records)"
else
    echo "Warning: no $CORRECTION_FILE -- TenX will produce no MS calls"
fi

export PROCESS_SAMPLES="$SAMPLE"
export RAW_CELLIDS_SAMPLES="$SAMPLE"
export RAW_CELLIDS_DATASETS="$RAW_CELLIDS_DATASET"
export RAW_CELLIDS_FOLDERS="$RAW_CELLIDS_FOLDER"

echo "Submitting SIV TenX stepThree for ${SAMPLE} (${N_ALIGN} alignment CSVs)"
sbatch \
    --job-name="f4sivtenx" \
    --output="logs/stepThree_process_SIV_TenX_%j.txt" \
    --error="logs/stepThree_process_SIV_TenX_%j.err" \
    --export=ALL \
    "$WORKER"

echo "SIV TenX stepThree job submitted."
