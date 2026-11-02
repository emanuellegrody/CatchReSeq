#!/bin/bash
# stepOne_mkcounts.sh  (Figure1 SIV Ca*Re)
# Reads stepOne_mkcounts.csv and submits one cellranger count sbatch job per row.
#
# Counts land under <OUTPUT_ROOT>/<Group>/counts/, with the group taken from the
# CSV Group column:
#   viremia     Pacute, W2
#   ART         D195
#   invitro     Invitro
#   uninfected  D0, W0
#
# Edit the paths below + the CSV, then run: bash stepOne_mkcounts.sh

set -euo pipefail

# Filesystem root for every input read and output written below. Override to run
# the pipeline against a different copy of the data:
#   DATA_ROOT=/my/data bash <this script>
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

### UPDATE ALL PATHS (shared across samples)
TRANSCRIPTOME_PATH="${DATA_ROOT}/genomes/Mmul_10_mac239"
OUTPUT_ROOT="${DATA_ROOT}/extractedData/Figure1/CaRe"
FASTQ_PATH=""    # per-row Fastq_Path takes priority; this is the fallback
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
fastq_path_index=$(get_col_index "Fastq_Path")
group_index=$(get_col_index "Group")
if [ -z "$sample_name_index" ]; then
    echo "Error: 'Sample_Name' column not found in $SAMPLES_CSV"
    exit 1
fi
if [ -z "$group_index" ]; then
    echo "Error: 'Group' column not found in $SAMPLES_CSV"
    exit 1
fi

tail -n +2 "$SAMPLES_CSV" | while IFS= read -r line; do
    [ -z "$line" ] && continue
    IFS=',' read -r -a fields <<< "$line"
    SAMPLE="$(echo "${fields[sample_name_index-1]}" | tr -d '[:space:]')"
    [ -z "$SAMPLE" ] && continue

    GROUP="$(echo "${fields[group_index-1]}" | tr -d '[:space:]')"
    if [ -z "$GROUP" ]; then
        echo "Error: empty Group for sample $SAMPLE"
        exit 1
    fi
    OUTPUT_PATH="${OUTPUT_ROOT}/${GROUP}/counts"

    # Per-row Fastq_Path takes priority, else fall back to wrapper-level FASTQ_PATH
    ROW_FASTQ=""
    if [ -n "$fastq_path_index" ]; then
        ROW_FASTQ="$(echo "${fields[fastq_path_index-1]}" | tr -d '[:space:]')"
    fi
    [ -z "$ROW_FASTQ" ] && ROW_FASTQ="$FASTQ_PATH"
    # The CSV path is relative to DATA_ROOT; an absolute path is taken as given.
    case "$ROW_FASTQ" in /*) ;; *) ROW_FASTQ="${DATA_ROOT}/$ROW_FASTQ" ;; esac
    if [ -z "$ROW_FASTQ" ]; then
        echo "Error: no FASTQ path for sample $SAMPLE"
        exit 1
    fi

    # Per-row Sample_Project takes priority, else fall back to wrapper-level PROJECT
    ROW_PROJECT=""
    if [ -n "$sample_project_index" ]; then
        ROW_PROJECT="$(echo "${fields[sample_project_index-1]}" | tr -d '[:space:]')"
    fi
    [ -z "$ROW_PROJECT" ] && ROW_PROJECT="$PROJECT"

    echo "Submitting mkcounts: $SAMPLE  group=${GROUP}  fastqs=${ROW_FASTQ}"
    sbatch \
        --job-name="mkcounts-${SAMPLE}" \
        --output="logs/stepOne_mkcounts_${SAMPLE}_%j.txt" \
        --error="logs/stepOne_mkcounts_${SAMPLE}_%j.err" \
        --export="ALL,SAMPLE=${SAMPLE},FASTQ_PATH=${ROW_FASTQ},TRANSCRIPTOME_PATH=${TRANSCRIPTOME_PATH},OUTPUT_PATH=${OUTPUT_PATH},CHEMISTRY=${CHEMISTRY},PROJECT=${ROW_PROJECT},OUTS_ONLY=${OUTS_ONLY}" \
        "$WORKER"
done

echo "All mkcounts jobs submitted."
