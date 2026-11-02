#!/bin/bash
# stepThree_process_SIV.sh  (Figure4 / EGS024 SIV Ca*Re amplicon library)
# Combines the stepTwo per-region alignment tables into one category-per-UMI
# table:
#   <CARE_ROOT>/bamsort/reads/GFP_invitro_UMI_read_counts_full.csv  (resolved)
#   <CARE_ROOT>/bamsort/reads/GFP_invitro_UMI_read_counts_raw.csv   (unresolved)
#
# The host targets have their own stepThree_process.sh in this directory; that
# one runs bamsort_process_host.R over PTPRC/BCL2L1/IL7R in the stepTwo/ tree.
# This is the viral arm and reads the bamsort/ tree.
#
# The alignment CSVs of both arms share bamsort/alignment/, but process_bamsort
# globs `.*<sample>.*_bamsort_alignment_.*`, and PTPRC/IL7R/BCL2L1 carry no
# sample name, so a GFP_invitro run picks up the 25 viral files only.
#
# RAW_CELLIDS_SAMPLES must name the same string as PROCESS_SAMPLES.
# process_bamsort gates the whole resolution block on
# `sample %in% names(raw_cellIDs)` and falls through to an unfiltered,
# UNRESOLVED table on a cat() warning when it does not. The preflight below
# turns that silent fallback into an error.
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
#   DATA_ROOT=/my/data bash stepThree_process_SIV.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

CARE_ROOT="${DATA_ROOT}/extractedData/Figure4/CaRe"
export R_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_process.R"
export FUNCTIONS_R="${PIPELINE_DIR}/bamsort/bamsort_functions.R"
export RAW_CELLIDS_OUTPUT_DIR="${DATA_ROOT}/extractedData/Figure4/singleCell/raw_cellIDs"
# <RAW_CELLIDS_FOLDER><GENOME_PREFIX>_<RAW_CELLIDS_DATASET>/outs/raw_feature_bc_matrix/
# is only rebuilt if the cache file is deleted; the cellranger run is
# Mmul_10_mac239_GFP_GEX.
export GENOME_PREFIX="Mmul_10_mac239_GFP"
export MIN_UMI_COUNT="100"

SAMPLE="GFP_invitro"
RAW_CELLIDS_DATASET="GEX"
RAW_CELLIDS_FOLDER="${DATA_ROOT}/extractedData/Figure4/singleCell/counts/concat/"

# Targets as they appear after the sample name in the alignment and splice
# filenames, i.e. GFP_invitro<TARGET>_bamsort_alignment_<region>.csv.
SPLICE_TARGETS=("_D1" "_nef" "_pol" "_SSenv" "_tat")

export PROCESS_INPUT_DIRS="${CARE_ROOT}/bamsort/alignment/"
export PROCESS_OUTPUT_DIR="${CARE_ROOT}/bamsort/reads/"
export PROCESS_CORRECTION_DIR="${CARE_ROOT}/bamsort/splice/combined/"

### DO NOT UPDATE BELOW HERE
WORKER="${PIPELINE_DIR}/bamsort/stepThree_process_worker.sh"

[ -f "$WORKER" ]      || { echo "Error: file not found: $WORKER"; exit 1; }
[ -f "$R_SCRIPT" ]    || { echo "Error: file not found: $R_SCRIPT"; exit 1; }
[ -d "$PROCESS_INPUT_DIRS" ] || { echo "Error: no alignment dir: $PROCESS_INPUT_DIRS"; exit 1; }

# Preflight: without this cache file process_bamsort skips resolve_multimap,
# apply_umi_corrections, the any-*/SS-MS* consolidation and the cell filter, and
# writes _full as a copy of _raw without failing.
CACHE="${RAW_CELLIDS_OUTPUT_DIR}/${SAMPLE}_raw_cellIDs.csv"
if [ ! -f "$CACHE" ]; then
    echo "Error: no raw cellID cache at $CACHE"
    echo "  Delete-and-rebuild path needs ${RAW_CELLIDS_FOLDER}${GENOME_PREFIX}_${RAW_CELLIDS_DATASET}/outs/raw_feature_bc_matrix/"
    exit 1
fi

N_ALIGN=$(ls "${PROCESS_INPUT_DIRS}"*"${SAMPLE}"*_bamsort_alignment_*.csv 2>/dev/null | wc -l)
[ "$N_ALIGN" -gt 0 ] || { echo "Error: no alignment CSVs matching $SAMPLE in $PROCESS_INPUT_DIRS"; exit 1; }

mkdir -p logs

# apply_umi_corrections reads a single <sample>_splicesites_filtered.csv, while
# bamsort_splice_sites.py writes one per target. Newer stepTwo workers
# concatenate them at the end of their own run; the stepTwo behind this tree
# predates that, so the combine happens here. Targets with no junctions (nef,
# pol) leave no filtered file and are skipped, as in the worker.
COMBINED_FILE="${PROCESS_CORRECTION_DIR%/}/${SAMPLE}_splicesites_filtered.csv"
mkdir -p "$PROCESS_CORRECTION_DIR"
: > "$COMBINED_FILE"
wrote_header=0
for TARGET_NAME in "${SPLICE_TARGETS[@]}"; do
    TARGET_FILE="${PROCESS_CORRECTION_DIR%/}/${SAMPLE}${TARGET_NAME}_splicesites_filtered.csv"
    [ -f "$TARGET_FILE" ] || continue
    if [ "$wrote_header" -eq 0 ]; then
        cat "$TARGET_FILE" >> "$COMBINED_FILE"
        wrote_header=1
    else
        tail -n +2 "$TARGET_FILE" >> "$COMBINED_FILE"
    fi
done
if [ "$wrote_header" -eq 0 ]; then
    echo "Warning: no per-target filtered splice files for $SAMPLE; removing empty $COMBINED_FILE"
    rm -f "$COMBINED_FILE"
else
    echo "Combined splice corrections: $COMBINED_FILE ($(($(wc -l < "$COMBINED_FILE") - 1)) records)"
fi

export PROCESS_SAMPLES="$SAMPLE"
export RAW_CELLIDS_SAMPLES="$SAMPLE"
export RAW_CELLIDS_DATASETS="$RAW_CELLIDS_DATASET"
export RAW_CELLIDS_FOLDERS="$RAW_CELLIDS_FOLDER"

echo "Submitting SIV Ca*Re stepThree for ${SAMPLE} (${N_ALIGN} alignment CSVs)"
sbatch \
    --job-name="f4sivThree" \
    --output="logs/stepThree_process_SIV_%j.txt" \
    --error="logs/stepThree_process_SIV_%j.err" \
    --export=ALL \
    "$WORKER"

echo "SIV Ca*Re stepThree job submitted."
