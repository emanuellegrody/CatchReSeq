#!/bin/bash
# stepTwo_bamsort_CaRe.sh  (EGS026 CaRe HIV clade amplicons, Figure3)
# Reads stepTwo_bamsort_CaRe.csv and submits one bamsort sbatch job per clade.
# Edit the paths below + the CSV, then run: bash stepTwo_bamsort_CaRe.sh
#
# Named apart from stepTwo_bamsort.sh in this directory, which is the single-cell
# HXB2 bamsort for the 10X pools; this one is the CaRe amplicon pipeline
# (alignment -> splice extract+filter -> R processing) and shares the Figure1
# worker.
#
# Run AFTER stepOne_mkcounts.sh (amplicon BAMs) and stepOne_cellranger.sh
# (single-cell raw_feature_bc_matrix, the source of raw_cellIDs) have finished.
# Follow with stepThree_process_CaRe.sh.
#
# The sample token carries a trailing underscore ("A_", not "A"), matching
# stepThree_process_CaRe.sh: process_bamsort selects a sample's alignment CSVs
# by substring match on the filename, and a bare "A" also matches the A1/A7
# coordinate targets inside other clades' files (C_env_bamsort_alignment_A1_D4.csv);
# "A_" matches only clade A's libraries. The same token keys raw_cellIDs and the
# splice-correction filename, so it is used consistently rather than overridden
# per step.

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can
# be cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

### UPDATE ALL PATHS (shared across samples)
# --- Step 1 (alignment) ---
export COORDS_CSV="${PIPELINE_DIR}/bamsort/coordinates/bamsort_coordinates_HXB2.csv"
export VIRAL_CONTIG="HXB2"
# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data bash stepTwo_bamsort_CaRe.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

export INPUT_ROOT="${DATA_ROOT}/extractedData/Figure3/CaRe/counts/GRCh38_HXB2_"
export INPUT_END="/outs/possorted_genome_bam.bam"
export OUTPUT_ROOT="${DATA_ROOT}/extractedData/Figure3/CaRe/bamsort/alignment/"
export ALIGNMENT_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_alignment_reads.py"
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


# --- Step 2 (splice extraction + filter) ---
# Splice BAMs are this same EGS026 amplicon set; targets carry no leading
# underscore because the sample token already ends in one.
export SPLICE_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_splice.py"
export SPLICE_FILTER_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_splice_sites.py"
export SPLICE_INPUT_ROOT="${DATA_ROOT}/extractedData/Figure3/CaRe/counts/GRCh38_HXB2_"
export SPLICE_INPUT_END="/outs/possorted_genome_bam.bam"
export SPLICE_OUTPUT_DIR="${DATA_ROOT}/extractedData/Figure3/CaRe/bamsort/splice"
export SPLICE_TARGETS="D1|D4|pol|env|nef"
export SPLICE_FILTER_DONOR="D4"
# Label junctions against the HXB2 splice sites, not the mac239 default.
export SPLICE_SITES="HXB2"

# PROCESS_CORRECTION_DIR is where the splice filter writes its confirmed
# junctions; stepThree_process_CaRe.sh reads them back from the same path.
export PROCESS_CORRECTION_DIR="${DATA_ROOT}/extractedData/Figure3/CaRe/bamsort/splice/combined/"

### DO NOT UPDATE BELOW HERE
SAMPLES_CSV="${SCRIPT_DIR}/stepTwo_bamsort_CaRe.csv"
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

tail -n +2 "$SAMPLES_CSV" | while IFS=, read -r SAMPLE ALIGNMENT_BAM_IDS RAW_CELLIDS_DATASET RAW_CELLIDS_FOLDER; do
    SAMPLE="$(echo "$SAMPLE" | tr -d '[:space:]')"
    [ -z "$SAMPLE" ] && continue

    RAW_CELLIDS_FOLDER="$(echo "${RAW_CELLIDS_FOLDER:-}" | tr -d '[:space:]\r')"
    # The CSV path is relative to DATA_ROOT; an absolute path is taken as given.
    case "$RAW_CELLIDS_FOLDER" in /*) ;; *) RAW_CELLIDS_FOLDER="${DATA_ROOT}/$RAW_CELLIDS_FOLDER" ;; esac

    # The worker still checks ALIGNMENT_BAM_ID; it runs over the full |-list.
    ALIGNMENT_BAM_ID="${ALIGNMENT_BAM_IDS%%|*}"

    export SAMPLE ALIGNMENT_BAM_ID ALIGNMENT_BAM_IDS RAW_CELLIDS_DATASET RAW_CELLIDS_FOLDER

    echo "Submitting bamsort: $SAMPLE  (${ALIGNMENT_BAM_IDS//|/, })"
    sbatch \
        --job-name="bamsort-${SAMPLE}" \
        --output="logs/stepTwo_bamsort_${SAMPLE}_%j.txt" \
        --error="logs/stepTwo_bamsort_${SAMPLE}_%j.err" \
        --export=ALL \
        "$WORKER"
done

echo "All bamsort jobs submitted."
