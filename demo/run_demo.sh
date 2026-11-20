#!/bin/bash
# run_demo.sh -- the Catch*ReSeq demo.
#
# Runs stepThree (UMI collapse and category assignment) and stepFour (detection
# thresholds) on the real D195 ART sample bundled in input/, then compares the
# result against expected_output/.
#
# No SLURM, no `module`, no cluster. R plus data.table, dplyr and tidyr is the
# whole requirement. See README.md in this directory for timings and for what
# each output column means.
#
#   bash run_demo.sh [work_dir]
#
# work_dir defaults to a temporary directory that is removed on exit; pass one
# to keep the intermediate tree for inspection.

set -euo pipefail

DEMO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${DEMO_DIR}/.." && pwd)"
PIPELINE_DIR="${REPO_DIR}/extractionScripts/pipelineScripts"

if [ $# -ge 1 ]; then
    WORK_DIR="$1"
    mkdir -p "$WORK_DIR"
    KEEP_WORK_DIR=1
else
    WORK_DIR="$(mktemp -d)"
    KEEP_WORK_DIR=0
    trap 'rm -rf "$WORK_DIR"' EXIT
fi

say() { printf '\n== %s\n' "$*"; }

# The pipeline addresses its inputs through the standard extractedData layout,
# so the demo stages its bundled copies into that shape rather than teaching the
# scripts a second set of paths.
ART_DIR="${WORK_DIR}/Figure1/CaRe/ART/bamsort"
mkdir -p "${ART_DIR}/alignment" "${ART_DIR}/splice/combined" \
         "${ART_DIR}/reads" "${WORK_DIR}/Figure1/singleCell/raw_cellIDs"
cp "${DEMO_DIR}"/input/alignment/*                "${ART_DIR}/alignment/"
cp "${DEMO_DIR}"/input/splice/*                   "${ART_DIR}/splice/combined/"
cp "${DEMO_DIR}"/input/raw_cellIDs/*              "${WORK_DIR}/Figure1/singleCell/raw_cellIDs/"

say "stepThree: collapsing reads to one category per UMI"
# raw_cellIDs are supplied precomputed, so the cellranger matrix is never read
# and Seurat is never needed. RAW_CELLIDS_DATASETS and RAW_CELLIDS_FOLDERS name
# where that matrix would have been found, and are unused on this path.
FUNCTIONS_R="${PIPELINE_DIR}/bamsort/bamsort_functions.R" \
RAW_CELLIDS_OUTPUT_DIR="${WORK_DIR}/Figure1/singleCell/raw_cellIDs" \
RAW_CELLIDS_DATASETS="LP29_D195" \
RAW_CELLIDS_SAMPLES="D195" \
RAW_CELLIDS_FOLDERS="${WORK_DIR}/Figure1/singleCell/counts/" \
PROCESS_SAMPLES="D195" \
PROCESS_INPUT_DIRS="${ART_DIR}/alignment/" \
PROCESS_OUTPUT_DIR="${ART_DIR}/reads/" \
PROCESS_CORRECTION_DIR="${ART_DIR}/splice/combined/" \
GENOME_PREFIX="Mmul_10_mac239" \
MIN_UMI_COUNT="100" \
    Rscript "${PIPELINE_DIR}/bamsort/bamsort_process.R"

say "stepFour: applying detection thresholds"
CARE_DATA_ROOT="${WORK_DIR}" \
    Rscript "${PIPELINE_DIR}/bamsort/stepFour_minimums_CaRe.R"

say "comparing against expected_output/"
STATUS=0
compare() {
    local produced="$1" expected="$2" label="$3"
    if [ ! -f "$produced" ]; then
        echo "MISSING  ${label}"
        STATUS=1
        return
    fi
    # filters.txt records the date it was written, so that one line is expected
    # to differ and is excluded from the comparison.
    if diff -q <(grep -v '^processed_date' "$produced") \
               <(grep -v '^processed_date' "$expected") >/dev/null; then
        echo "OK       ${label}"
    else
        echo "DIFFERS  ${label}"
        diff <(grep -v '^processed_date' "$expected") \
             <(grep -v '^processed_date' "$produced") | head -20
        STATUS=1
    fi
}

READS_DIR="${ART_DIR}/reads"
compare "${READS_DIR}/D195_UMI_read_counts_raw.csv" \
        "${DEMO_DIR}/expected_output/D195_UMI_read_counts_raw.csv"  "D195_UMI_read_counts_raw.csv"
compare "${READS_DIR}/D195_UMI_read_counts_full.csv" \
        "${DEMO_DIR}/expected_output/D195_UMI_read_counts_full.csv" "D195_UMI_read_counts_full.csv"
compare "${READS_DIR}/minimum/D195_CaRe_filtered.csv" \
        "${DEMO_DIR}/expected_output/minimum/D195_CaRe_filtered.csv" "minimum/D195_CaRe_filtered.csv"
compare "${READS_DIR}/minimum/D195_filters.txt" \
        "${DEMO_DIR}/expected_output/minimum/D195_filters.txt"       "minimum/D195_filters.txt"

if [ "$KEEP_WORK_DIR" = "1" ]; then
    say "outputs left in ${READS_DIR}"
fi

if [ "$STATUS" = "0" ]; then
    say "demo passed"
else
    say "demo FAILED -- see the differences above"
fi
exit "$STATUS"
