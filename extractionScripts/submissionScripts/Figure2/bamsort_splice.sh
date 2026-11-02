#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort_splice
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=5G
#SBATCH -t 4:00:00
#SBATCH --output=logs/bamsort_splice_%j.txt
#SBATCH --error=logs/bamsort_splice_%j.err
#SBATCH --verbose

# Script locations are resolved from this script's own path so the repo can
# be cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

# bamsort_splice.sh  (Figure2 Panels C and J: TenX splice junctions)
# Writes <output_dir>/<sample>_splicesites.csv, the raw per-junction tables
# Figure2.Rmd reads for the TenX arm of the donor->acceptor maps (Panel C) and the
# Panel J recovery rate.
#
# Copied into CatchReSeq from extractionScripts/submissionScripts/Figure1/singleCell/
# bamsort_splice_invitro.sh + bamsort_splice_20251102.sh and EGS009/
# bamsort_splice_20251102.sh, which are the same script run once per sample. The
# three TenX samples are combined here because they are one figure's inputs, and
# they span two datasets, so INPUT_DIR is per-sample rather than global.
#
# The Ca*Re arm of the same panels comes from submissionScripts/Figure1/ stepTwo.
# MIN_INTRON must match SPLICE_MIN_INTRON there or the two arms are not
# comparable: the default of 1000 hides the D3 acceptors (D3-A3 528, D3-A4 600,
# D3-A5 741) that 500 recovers.
#
# Run: sbatch bamsort_splice.sh

module purge
source activate cellranger
module load python
module load samtools

SPLICE_SCRIPT="${PIPELINE_DIR}/bamsort/bamsort_splice.py"
MIN_INTRON=500
VIRAL_CONTIG="mac239"

# sample : BAM path : output dir. P_acute_combined sits under Mmul_10_mac239v5_,
# not Mmul_10_mac239_, but carries the same mac239 contig (LN:9602) as the rest.
# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data bash bamsort_splice.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

EXTRACTED="${DATA_ROOT}/extractedData"
SAMPLES=(
    "invitro|${EXTRACTED}/Figure1/singleCell/counts/Mmul_10_mac239/Mmul_10_mac239_invitro|${EXTRACTED}/Figure1/singleCell/bamsort/splice"
    "W2|${EXTRACTED}/Figure1/singleCell/counts/Mmul_10_mac239/Mmul_10_mac239_W2|${EXTRACTED}/Figure1/singleCell/bamsort/splice"
    "P_acute_combined|${EXTRACTED}/Figure1/singleCell/counts/Mmul_10_mac239v5_P_acute_combined|${EXTRACTED}/Figure1/singleCell/bamsort/splice"
)

for ENTRY in "${SAMPLES[@]}"; do
    IFS='|' read -r SAMPLE_NAME COUNTS_DIR OUTPUT_DIR <<< "$ENTRY"

    INPUT_BAM="${COUNTS_DIR}/outs/possorted_genome_bam.bam"
    if [ ! -f "$INPUT_BAM" ]; then
        echo "Warning: Input BAM file not found: $INPUT_BAM"
        continue
    fi

    mkdir -p "${OUTPUT_DIR}"
    OUTPUT_CSV="${OUTPUT_DIR}/${SAMPLE_NAME}_splicesites.csv"

    echo "Processing sample: $SAMPLE_NAME"
    echo "Input BAM: $INPUT_BAM"
    echo "Output CSV: $OUTPUT_CSV"
    echo "Min intron: $MIN_INTRON"

    python3 "$SPLICE_SCRIPT" \
        "$INPUT_BAM" \
        "$OUTPUT_CSV" \
        --reference "$VIRAL_CONTIG" \
        --min-intron "$MIN_INTRON"
done
