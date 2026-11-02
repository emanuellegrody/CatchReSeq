#!/bin/bash
# stepTwo_bamsort_TenX.sh  (Figure1 single-cell TenX libraries)
# Reads stepTwo_bamsort_TenX.csv and submits one bamsort sbatch job per sample:
# per-region alignment tables + splice junctions + the canonical-donor filter,
# over the one whole-transcriptome cellranger BAM of that sample.
#
# This is the TenX arm of the same comparison stepTwo_bamsort.sh runs on the
# Ca*Re amplicon libraries. It exists as a separate wrapper because a TenX sample
# is a single library with no target suffix, which the amplicon worker cannot
# express: it builds every path as <root><sample><target><end> over a non-empty
# target list. Everything downstream is shared.
#
# The TenX BAMs are not in one tree — they are the original cellranger runs from
# EGS002, EGS004 and EGS009 — so each row of the CSV carries its own counts
# dir, relative to DATA_ROOT, rather than a shared root. All of them align to
# the same mac239 contig (LN:9602, cellranger 7.2.0) as the Ca*Re BAMs, which
# is what makes the region coordinates in bamsort_coordinates_TenX.csv valid
# for both arms.
#
# SPLICE_MIN_INTRON and SPLICE_FILTER_DONOR must match stepTwo_bamsort.sh, or
# the TenX and Ca*Re junction sets are not comparable.
#
# Edit the paths below + the CSV, then run: bash stepTwo_bamsort_TenX.sh
#
# Follow with stepThree_process_TenX.sh, which combines the per-region tables
# into one category-per-UMI table per sample and applies the splice correction.

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can be
# cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

# Filesystem root for every input read and output written below. Override to run
# the pipeline against a different copy of the data:
#   DATA_ROOT=/my/data bash <this script>
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

### UPDATE ALL PATHS (shared across samples)
TenX_ROOT="${DATA_ROOT}/extractedData/Figure1/singleCell"
PIPELINE_DIR="${PIPELINE_DIR}"

export COORDS_CSV="${PIPELINE_DIR}/bamsort/coordinates/bamsort_coordinates_SIVmac239.csv"
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

export SPLICE_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_splice.py"
export SPLICE_FILTER_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_splice_sites.py"

export VIRAL_CONTIG="mac239"
export SPLICE_FILTER_DONOR="D4"
export SPLICE_MIN_INTRON="500"

export OUTPUT_ROOT="${TenX_ROOT}/bamsort/alignment/"
export SPLICE_OUTPUT_DIR="${TenX_ROOT}/bamsort/splice"
export PROCESS_CORRECTION_DIR="${TenX_ROOT}/bamsort/splice/combined/"

### DO NOT UPDATE BELOW HERE
SAMPLES_CSV="${SCRIPT_DIR}/stepTwo_bamsort_TenX.csv"
WORKER="${PIPELINE_DIR}/bamsort/stepTwo_bamsort_worker_TenX.sh"

for f in "$SAMPLES_CSV" "$WORKER" "$COORDS_CSV"; do
    [ -f "$f" ] || { echo "Error: file not found: $f"; exit 1; }
done

mkdir -p logs

tail -n +2 "$SAMPLES_CSV" | while IFS=, read -r SAMPLE COUNTS_DIR _RAW_CELLIDS_DATASET _RAW_CELLIDS_FOLDER; do
    SAMPLE="$(echo "$SAMPLE" | tr -d '[:space:]\r')"
    COUNTS_DIR="$(echo "$COUNTS_DIR" | tr -d '[:space:]\r')"
    # The CSV path is relative to DATA_ROOT; an absolute path is taken as given.
    case "$COUNTS_DIR" in /*) ;; *) COUNTS_DIR="${DATA_ROOT}/$COUNTS_DIR" ;; esac
    [ -z "$SAMPLE" ] && continue

    if [ -z "$COUNTS_DIR" ]; then
        echo "Error: row for sample '$SAMPLE' has no counts_dir"
        exit 1
    fi

    export ALIGNMENT_BAM="${COUNTS_DIR}/outs/possorted_genome_bam.bam"
    if [ ! -f "$ALIGNMENT_BAM" ]; then
        echo "Error: BAM not found for $SAMPLE: $ALIGNMENT_BAM"
        exit 1
    fi

    export SAMPLE

    echo "Submitting TenX bamsort: $SAMPLE  ($(basename "$COUNTS_DIR"))"
    sbatch \
        --job-name="tenxbam-${SAMPLE}" \
        --output="logs/stepTwo_bamsort_TenX_${SAMPLE}_%j.txt" \
        --error="logs/stepTwo_bamsort_TenX_${SAMPLE}_%j.err" \
        --export=ALL \
        "$WORKER"
done

echo "All TenX bamsort jobs submitted."
