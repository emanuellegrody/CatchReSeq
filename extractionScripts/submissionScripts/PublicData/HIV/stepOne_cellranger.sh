#!/bin/bash
# stepOne_cellranger.sh
# Submit one cellranger job per TenX sample (10X 5' v2, vs joint human GRCh38 +
# HXB2) for the Collora et al. 2022 HIV data. Aligning to the combined reference
# maps host + HIV in one pass, so HIV reads get CB/UB/xf tags directly.
# Reads samples_manifest.csv (SAMPLE,GSM,RUNS),
# where RUNS is a '|'-separated list of SRA runs (lanes) for that sample, all
# combined into one cellranger run.
#
# Only the 30 scRNA-seq (TenX) samples are listed; ADT/VDJ libraries are excluded
# (no HIV cDNA). Manifest was built from ENA metadata for PRJNA777601.
# Edit paths below, then: bash stepOne_cellranger.sh

set -euo pipefail

### UPDATE PATHS
# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data GENOME_ROOT=/my/genomes bash stepOne_cellranger.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"
GENOME_ROOT="${GENOME_ROOT:-/projects/b1042/GoyalLab/egrody/genomes}"

export FASTQ_SRC_DIR="${DATA_ROOT}/publicData/Colloraetal2022/fastq"
export TRANSCRIPTOME="${GENOME_ROOT}/GRCh38_HXB2"
export OUTPUT_PATH="${DATA_ROOT}/publicData/Colloraetal2022/counts"
export OUTS_ONLY="FALSE"

### DO NOT UPDATE BELOW
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${SCRIPT_DIR}/samples_manifest.csv"
WORKER="${SCRIPT_DIR}/../../../pipelineScripts/cellranger/hiv_cellranger_worker.sh"

[ -f "$MANIFEST" ] || { echo "Missing $MANIFEST"; exit 1; }
[ -f "$WORKER" ]   || { echo "Missing $WORKER"; exit 1; }
mkdir -p "${SCRIPT_DIR}/logs"

# Optional filter: pass sample names as args to submit only those; else submit all.
WANTED=" $* "

tail -n +2 "$MANIFEST" | while IFS=, read -r SAMPLE GSM RUNS; do
    SAMPLE="$(echo "$SAMPLE" | tr -d '[:space:]')"
    [ -z "$SAMPLE" ] && continue
    if [ "$#" -gt 0 ] && [[ "$WANTED" != *" $SAMPLE "* ]]; then continue; fi

    export SAMPLE RUNS
    echo "Submitting cellranger: $SAMPLE ($GSM)  lanes: ${RUNS//|/ }"
    sbatch \
        --job-name="cr-${SAMPLE}" \
        --output="${SCRIPT_DIR}/logs/stepOne_${SAMPLE}_%j.txt" \
        --error="${SCRIPT_DIR}/logs/stepOne_${SAMPLE}_%j.err" \
        --export=ALL \
        "$WORKER"
done
echo "Done submitting cellranger jobs."
