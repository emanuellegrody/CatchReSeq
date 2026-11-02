#!/bin/bash
# stepOne_mkcounts_KLRB1.sh  (Figure1 KLRB1 Ca*Re )
# Reads stepOne_mkcounts_KLRB1.csv and submits cellranger count sbatch job.
#
# Counts land at
#   <DATA_ROOT>/<OUTPUT_REL>/Mmul_10_mac239_<Sample_Name>/outs/
# which is where stepTwo_bamsort_KLRB1.csv expects the Ca*Re BAMs. The TenX arm
# needs no stepOne: it reads the original EGS002 cellranger runs.
#
# stepOne_mkcounts_KLRB1.csv columns:
#   Sample_Name  — cellranger --sample, matching the FASTQ filename prefix
#   Fastq_Path   — dir holding the FASTQs, relative to DATA_ROOT
#
# Edit the paths below + the CSV, then run: bash stepOne_mkcounts_KLRB1.sh

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can be
# cloned anywhere
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

# Filesystem root for every input read and output written below. Override to run
# the pipeline against a different copy of the data:
#   DATA_ROOT=/my/data bash stepOne_mkcounts_KLRB1.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

### UPDATE ALL PATHS (shared across samples)
TRANSCRIPTOME_PATH="${DATA_ROOT}/genomes/Mmul_10_mac239"
OUTPUT_REL="extractedData/Figure1/CaRe/counts/Mmul_10_mac239"
CHEMISTRY="v3"    # 'v3' (cellranger 7.2.0) or 'v4' (cellranger 10.0.0)
OUTS_ONLY="FALSE" # 'TRUE' keeps only [genome]_[sample]/outs/ after cellranger count

### DO NOT UPDATE BELOW HERE
SAMPLES_CSV="${SCRIPT_DIR}/stepOne_mkcounts_KLRB1.csv"
WORKER="${PIPELINE_DIR}/cellranger/stepOne_mkcounts_worker.sh"

if [ "$CHEMISTRY" != "v3" ] && [ "$CHEMISTRY" != "v4" ]; then
    echo "Error: CHEMISTRY must be 'v3' or 'v4', got '$CHEMISTRY'"
    exit 1
fi
if [ "$OUTS_ONLY" != "TRUE" ] && [ "$OUTS_ONLY" != "FALSE" ]; then
    echo "Error: OUTS_ONLY must be 'TRUE' or 'FALSE', got '$OUTS_ONLY'"
    exit 1
fi
for f in "$SAMPLES_CSV" "$WORKER"; do
    if [ ! -f "$f" ]; then
        echo "Error: required file not found: $f"
        exit 1
    fi
done
if [ ! -d "$TRANSCRIPTOME_PATH" ]; then
    echo "Error: transcriptome not found: $TRANSCRIPTOME_PATH"
    exit 1
fi

OUTPUT_PATH="${DATA_ROOT}/${OUTPUT_REL}"
mkdir -p logs

tail -n +2 "$SAMPLES_CSV" | while IFS=, read -r SAMPLE FASTQ_REL; do
    SAMPLE="$(echo "$SAMPLE" | tr -d '[:space:]\r')"
    FASTQ_REL="$(echo "${FASTQ_REL:-}" | tr -d '[:space:]\r')"
    [ -z "$SAMPLE" ] && continue

    if [ -z "$FASTQ_REL" ]; then
        echo "Error: no Fastq_Path for sample $SAMPLE"
        echo "       expected CSV header: Sample_Name,Fastq_Path"
        exit 1
    fi

    FASTQ_PATH="${DATA_ROOT}/${FASTQ_REL}"
    if [ ! -d "$FASTQ_PATH" ]; then
        echo "Error: FASTQ dir not found for '$SAMPLE': $FASTQ_PATH"
        exit 1
    fi

    echo "Submitting mkcounts: $SAMPLE  fastqs=${FASTQ_PATH}"
    sbatch \
        --job-name="mkcntK-${SAMPLE}" \
        --output="logs/stepOne_mkcounts_KLRB1_${SAMPLE}_%j.txt" \
        --error="logs/stepOne_mkcounts_KLRB1_${SAMPLE}_%j.err" \
        --export="ALL,SAMPLE=${SAMPLE},FASTQ_PATH=${FASTQ_PATH},TRANSCRIPTOME_PATH=${TRANSCRIPTOME_PATH},OUTPUT_PATH=${OUTPUT_PATH},CHEMISTRY=${CHEMISTRY},OUTS_ONLY=${OUTS_ONLY}" \
        "$WORKER"
done

echo "All KLRB1 mkcounts jobs submitted."
