#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=70G
#SBATCH -t 10:00:00
#SBATCH --verbose

# Per-sample worker: runs cellranger count for one sample.
# Submitted by Figure1/stepOne_mkcounts.sh; not run directly.
#
# Required env vars (from wrapper):
#   SAMPLE              sample name (matches FASTQ prefix)
#   FASTQ_PATH          dir containing FASTQs
#   TRANSCRIPTOME_PATH  cellranger reference dir
#   OUTPUT_PATH         dir to cd into before cellranger count
#   CHEMISTRY           'v3' (cellranger 7.2.0) or 'v4' (cellranger 10.0.0)
# Optional env vars:
#   PROJECT             BCL Convert project folder name; passed as --project
#   OUTS_ONLY           'TRUE' to keep only outs/ inside the run dir after cellranger count; default FALSE
#   CELLRANGER_V3_BIN   cellranger 7.2.0 binary, used when CHEMISTRY=v3
#   CELLRANGER_V4_BIN   cellranger 10.0.0 binary, used when CHEMISTRY=v4

set -euo pipefail

for v in SAMPLE FASTQ_PATH TRANSCRIPTOME_PATH OUTPUT_PATH CHEMISTRY; do
    if [ -z "${!v:-}" ]; then
        echo "Error: required env var $v not set"
        exit 1
    fi
done

case "$CHEMISTRY" in
    v3)
        cellranger_bin="${CELLRANGER_V3_BIN:-/home/egy2296/packages/cellranger-7.2.0/cellranger}"
        bam_flag=""
        chemistry_flag="--chemistry=SC3Pv3"
        ;;
    v4)
        cellranger_bin="${CELLRANGER_V4_BIN:-/home/egy2296/packages/cellranger-10.0.0/cellranger}"
        bam_flag="--create-bam=true"
        chemistry_flag="--chemistry=SC3Pv4"
        ;;
    *)
        echo "Error: CHEMISTRY must be 'v3' or 'v4', got '$CHEMISTRY'"
        exit 1
        ;;
esac

project_flag=""
if [ -n "${PROJECT:-}" ]; then
    project_flag="--project=$PROJECT"
fi

source activate cellranger

mkdir -p "$OUTPUT_PATH"
cd "$OUTPUT_PATH"
genome_name=$(basename "$TRANSCRIPTOME_PATH")

"$cellranger_bin" count \
    --id="${genome_name}_${SAMPLE}" \
    --fastqs="$FASTQ_PATH" \
    --sample="$SAMPLE" \
    --transcriptome="$TRANSCRIPTOME_PATH" \
    --nosecondary \
    $chemistry_flag $bam_flag $project_flag

if [ "${OUTS_ONLY:-FALSE}" = "TRUE" ]; then
    run_dir="${OUTPUT_PATH}/${genome_name}_${SAMPLE}"
    echo "OUTS_ONLY=TRUE: removing non-outs entries under $run_dir"
    find "$run_dir" -mindepth 1 -maxdepth 1 ! -name outs -exec rm -rf {} +
fi
