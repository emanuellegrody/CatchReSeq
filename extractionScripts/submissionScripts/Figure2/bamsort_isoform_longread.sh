#!/bin/bash
# bamsort_isoform_longread.sh  (Figure2 Panel D, long-read isoforms)
# Calls isoforms per long-read molecule for the ONT cDNA library and
# writes cDNA_isoforms.csv, which Figure2.Rmd reads for the short-read vs
# long-read comparison and the spliced-isoform fraction.
#
# Run: sbatch bamsort_isoform_longread.sh
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort_isoform
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=1G
#SBATCH -t 0:10:00
#SBATCH --output=logs/bamsort_isoform_%j.txt
#SBATCH --error=logs/bamsort_isoform_%j.err
#SBATCH --verbose

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can be
# cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

# Filesystem root for every input read and output written below. Override to run
# against a different copy of the data:
#   DATA_ROOT=/my/data bash bamsort_isoform_longread.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

module purge
source activate cellranger
module load python
module load samtools

# Base directory for input BAM files
INPUT_DIR="${DATA_ROOT}/extractedData/Figure2/longRead/nf-core/cDNA/bam/dedup"
OUTPUT_DIR="${DATA_ROOT}/extractedData/Figure2/longRead/bamsort/isoform"

# List of sample names (space-separated)
SAMPLES=("cDNA")

mkdir -p "${OUTPUT_DIR}"

for SAMPLE_NAME in "${SAMPLES[@]}"; do
    # Input BAM file path
    INPUT_BAM="${INPUT_DIR}/${SAMPLE_NAME}.dedup.sorted.bam"
    if [ ! -f "$INPUT_BAM" ]; then
        echo "Warning: Input BAM file not found: $INPUT_BAM"
        continue
    fi

    # Output CSV filename
    OUTPUT_FILENAME="${SAMPLE_NAME}_isoforms.csv"
    OUTPUT_CSV="${OUTPUT_DIR}/${OUTPUT_FILENAME}"
    
    echo "Processing sample: $SAMPLE_NAME"
    echo "Input BAM: $INPUT_BAM"
    echo "Output CSV: $OUTPUT_CSV"
    
    # The donor/acceptor coordinates in bamsort_isoform_longread.py are written
    # in the 9602bp mac239.fa used everywhere else in this project. This BAM was
    # aligned by nf-core against a 10279bp mac239, whose coordinates run 553bp
    # higher (measured from junction-spanning reads: the dominant gap at
    # 6597-8802 is canonical D4/6043-A7/8249). Without the offset the D4-A7
    # region test fails and multiply-spliced molecules are called "any".
    # --merge-supplementary pools each read's primary and supplementary blocks,
    # since minimap2 splits some spliced long reads across both.
    python3 "${PIPELINE_DIR}/bamsort/bamsort_isoform_longread.py" \
        "$INPUT_BAM" \
        "$OUTPUT_CSV" \
        --coordinate-offset 553 \
        --merge-supplementary
done
