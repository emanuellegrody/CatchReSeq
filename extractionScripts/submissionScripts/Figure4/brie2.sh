#!/bin/bash
# brie2.sh  (Figure4 BRIE2)
# Reads brie2.csv and submits one BRIE2 sbatch job per row: brie-count,
# brie-quant, then the per-cell isoform CSV the .Rmd reads.
#
# Replaces the two near-identical brie2_all.sh / brie2_each.sh pairs that lived
# under submissionScripts/EGS024 (macaque host targets) and EGL002 (mouse liver
# targets). The two differed only in the GFF, the sample ids, the output naming
# and whether the barcode-filtering step ran, so all four are now columns in
# brie2.csv.
#
# Usage:
#   bash brie2.sh                # every row
#   bash brie2.sh invitro        # one row
#   bash brie2.sh liver_WT liver_KO
#
# Requires the `brie2` conda environment (brie-count / brie-quant on PATH in it).

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can be
# cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

# Filesystem roots for every input read and output written below. Override to
# run against a different copy of the data:
#   GENOME_ROOT=/my/genomes bash brie2.sh
GENOME_ROOT="${GENOME_ROOT:-/projects/b1042/GoyalLab/egrody/genomes}"

SAMPLES_CSV="${SCRIPT_DIR}/brie2.csv"
WORKER="${PIPELINE_DIR}/isoform/brie2_worker.sh"
FILTER="${PIPELINE_DIR}/isoform/filter_count_matrices.py"

for f in "$SAMPLES_CSV" "$WORKER"; do
    [ -f "$f" ] || { echo "Error: required file not found: $f"; exit 1; }
done

mkdir -p logs

# DATA_ROOT already names extractedData, so the CSV's paths are relative to its
# parent; strip the trailing component once rather than in every row.
ROOT="$(dirname "$DATA_ROOT")"

submitted=0
while IFS=, read -r row counts_dir sample out_name gff output_dir barcode_csv export_csv events; do
    row="$(echo "$row" | tr -d '[:space:]\r')"
    [ -z "$row" ] && continue
    case "$row" in \#*|row) continue ;; esac
    if [ $# -gt 0 ]; then
        want=0
        for a in "$@"; do [ "$a" = "$row" ] && want=1; done
        [ $want -eq 1 ] || continue
    fi

    COUNTS="${ROOT}/$(echo "$counts_dir" | tr -d '[:space:]\r')"
    SAMPLE="$(echo "$sample" | tr -d '[:space:]\r')"
    OUT_NAME="$(echo "$out_name" | tr -d '[:space:]\r')"
    GFF="${GENOME_ROOT}/$(echo "$gff" | tr -d '[:space:]\r')"
    OUTDIR="${ROOT}/$(echo "$output_dir" | tr -d '[:space:]\r')"
    BARCODE_CSV="$(echo "${barcode_csv:-}" | tr -d '[:space:]\r')"
    EXPORT_CSV="$(echo "$export_csv" | tr -d '[:space:]\r')"
    EVENTS="$(echo "${events:-}" | tr -d '\r')"

    [ -f "$GFF" ] || { echo "Error [$row]: GFF not found: $GFF"; exit 1; }
    [ -d "${COUNTS}/${SAMPLE}" ] || { echo "Error [$row]: counts dir not found: ${COUNTS}/${SAMPLE}"; exit 1; }
    [ -f "${COUNTS}/${SAMPLE}/outs/possorted_genome_bam.bam" ] \
        || { echo "Error [$row]: BAM not found under ${COUNTS}/${SAMPLE}/outs/"; exit 1; }

    mkdir -p "${OUTDIR}/${OUT_NAME}"

    # Restrict BRIE2 to a target-gene cell list when the row names one. The
    # worker picks up filtered_barcodes.tsv from the output dir if it is there.
    if [ -n "$BARCODE_CSV" ]; then
        CSV_PATH="${ROOT}/${BARCODE_CSV}"
        [ -f "$CSV_PATH" ] || { echo "Error [$row]: barcode CSV not found: $CSV_PATH"; exit 1; }
        [ -f "$FILTER" ]   || { echo "Error [$row]: $FILTER not found"; exit 1; }
        echo "[$row] filtering barcodes from $(basename "$CSV_PATH")"
        python "$FILTER" "$CSV_PATH" "$OUT_NAME" "$COUNTS" "$OUTDIR"
    fi

    echo "Submitting BRIE2: $row  ($SAMPLE -> ${OUTDIR}/${OUT_NAME})"
    sbatch \
        --job-name="brie2-${row}" \
        --output="logs/brie2_${row}_%j.txt" \
        --error="logs/brie2_${row}_%j.err" \
        --export=ALL \
        "$WORKER" "$SAMPLE" "$COUNTS" "$OUTDIR" "$GFF" "$OUT_NAME" "$EXPORT_CSV" $EVENTS
    submitted=$((submitted + 1))
done < <(grep -v '^#' "$SAMPLES_CSV")

if [ "$submitted" -eq 0 ]; then
    echo "No rows selected. Rows in brie2.csv:"
    grep -v '^#' "$SAMPLES_CSV" | tail -n +2 | cut -d, -f1 | sed 's/^/  /'
    exit 1
fi
echo "Submitted $submitted BRIE2 job(s)."
