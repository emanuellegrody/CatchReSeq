#!/bin/bash
# stepOne_mkcounts.sh  (EGS026 CaRe HIV clade amplicons, Figure3)
# Reads stepOne_mkcounts.csv and submits one cellranger count sbatch job per row.
# Covers the 15 HIV clade A/C/D amplicon libraries in fastq_concat/dualIndex
# (both sequencing runs concatenated per sample), aligned
# to the joint human GRCh38 + HXB2 reference so the bamsort step can pull HXB2
# reads with CB/UB tags. The SIV amplicons in that same FASTQ dir belong to the
# Mmul_10_mac239 sample set and are not covered here.
# Edit the paths below + the CSV, then run: bash stepOne_mkcounts.sh

set -euo pipefail

### UPDATE ALL PATHS (shared across samples)
# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data GENOME_ROOT=/my/genomes bash stepOne_mkcounts.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"
GENOME_ROOT="${GENOME_ROOT:-/projects/b1042/GoyalLab/egrody/genomes}"

FASTQ_PATH="${DATA_ROOT}/rawData/EGS026/fastq_concat/dualIndex"
TRANSCRIPTOME_PATH="${GENOME_ROOT}/GRCh38_HXB2"
OUTPUT_PATH="${DATA_ROOT}/extractedData/Figure3/CaRe/counts"
CHEMISTRY="v3"   # 'v3' (cellranger 7.2.0) or 'v4' (cellranger 10.0.0)
PROJECT=""       # leave "" if FASTQ_PATH has only one project
OUTS_ONLY="FALSE" # 'TRUE' keeps only [genome]_[sample]/outs/ after cellranger count; 'FALSE' keeps the full run dir

### DO NOT UPDATE BELOW HERE
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLES_CSV="${SCRIPT_DIR}/stepOne_mkcounts.csv"
WORKER="${SCRIPT_DIR}/../../pipelineScripts/cellranger/stepOne_mkcounts_worker.sh"

if [ "$CHEMISTRY" != "v3" ] && [ "$CHEMISTRY" != "v4" ]; then
    echo "Error: CHEMISTRY must be 'v3' or 'v4', got '$CHEMISTRY'"
    exit 1
fi
if [ "$OUTS_ONLY" != "TRUE" ] && [ "$OUTS_ONLY" != "FALSE" ]; then
    echo "Error: OUTS_ONLY must be 'TRUE' or 'FALSE', got '$OUTS_ONLY'"
    exit 1
fi
if [ ! -f "$SAMPLES_CSV" ]; then
    echo "Error: samples CSV not found: $SAMPLES_CSV"
    exit 1
fi
if [ ! -f "$WORKER" ]; then
    echo "Error: worker script not found: $WORKER"
    exit 1
fi

mkdir -p logs

# Find column indices (matches mkcounts_all.sh schema)
get_col_index() {
    awk -F ',' -v col="$1" 'NR==1 {for (i=1; i<=NF; i++) if ($i == col) print i; exit}' "$SAMPLES_CSV"
}
sample_name_index=$(get_col_index "Sample_Name")
sample_project_index=$(get_col_index "Sample_Project")
if [ -z "$sample_name_index" ]; then
    echo "Error: 'Sample_Name' column not found in $SAMPLES_CSV"
    exit 1
fi

tail -n +2 "$SAMPLES_CSV" | while IFS= read -r line; do
    [ -z "$line" ] && continue
    IFS=',' read -r -a fields <<< "$line"
    SAMPLE="$(echo "${fields[sample_name_index-1]}" | tr -d '[:space:]')"
    [ -z "$SAMPLE" ] && continue

    # Per-row Sample_Project takes priority, else fall back to wrapper-level PROJECT
    ROW_PROJECT=""
    if [ -n "$sample_project_index" ]; then
        ROW_PROJECT="$(echo "${fields[sample_project_index-1]}" | tr -d '[:space:]')"
    fi
    [ -z "$ROW_PROJECT" ] && ROW_PROJECT="$PROJECT"

    echo "Submitting mkcounts: $SAMPLE  project=${ROW_PROJECT:-none}"
    sbatch \
        --job-name="mkcounts-${SAMPLE}" \
        --output="logs/stepOne_mkcounts_${SAMPLE}_%j.txt" \
        --error="logs/stepOne_mkcounts_${SAMPLE}_%j.err" \
        --export="ALL,SAMPLE=${SAMPLE},FASTQ_PATH=${FASTQ_PATH},TRANSCRIPTOME_PATH=${TRANSCRIPTOME_PATH},OUTPUT_PATH=${OUTPUT_PATH},CHEMISTRY=${CHEMISTRY},PROJECT=${ROW_PROJECT},OUTS_ONLY=${OUTS_ONLY}" \
        "$WORKER"
done

echo "All mkcounts jobs submitted."
