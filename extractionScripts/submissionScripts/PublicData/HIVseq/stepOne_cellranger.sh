#!/bin/bash
# stepOne_cellranger.sh  (HIVseq / PRJNA1107214)
# Submit one cellranger job per TenX sample (10X 5' v2, vs joint human GRCh38 + HXB2).
# Reads samples_manifest.csv (SAMPLE,GSM,RUNS,ORIGINAL_GEO_SAMPLE); RUNS is
# '|'-separated SRA runs. SAMPLE is the corrected library name and is what the
# cellranger run dir is named after; ORIGINAL_GEO_SAMPLE is the name the GSM was
# filed under in GEO, which differs for six PID8027 rows whose SRA runs GEO
# attached to the wrong library (see the PID8027 section of Figure3.Rmd).
# (complete lanes only). The worker auto-detects the barcode read (_1 or _3) and
# cDNA read (_2 or _4) per run. Only the 16 samples with a downloaded barcode
# read are listed; 8 GSMs (PID145*, PID1052*) are missing the barcode read and
# need re-download first.
# Edit paths below, then: bash stepOne_cellranger.sh  [SAMPLE ...]

set -euo pipefail

### UPDATE PATHS
# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data GENOME_ROOT=/my/genomes bash stepOne_cellranger.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"
GENOME_ROOT="${GENOME_ROOT:-/projects/b1042/GoyalLab/egrody/genomes}"

export FASTQ_SRC_DIR="${DATA_ROOT}/publicData/HIVseq/fastq"
export TRANSCRIPTOME="${GENOME_ROOT}/GRCh38_HXB2"
export OUTPUT_PATH="${DATA_ROOT}/publicData/HIVseq/counts"
export OUTS_ONLY="FALSE"

### DO NOT UPDATE BELOW
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${SCRIPT_DIR}/samples_manifest.csv"
WORKER="${SCRIPT_DIR}/../../../pipelineScripts/cellranger/hiv_cellranger_worker.sh"

[ -f "$MANIFEST" ] || { echo "Missing $MANIFEST"; exit 1; }
[ -f "$WORKER" ]   || { echo "Missing $WORKER"; exit 1; }
mkdir -p "${SCRIPT_DIR}/logs"

WANTED=" $* "
tail -n +2 "$MANIFEST" | while IFS=, read -r SAMPLE GSM RUNS ORIGINAL_GEO_SAMPLE; do
    SAMPLE="$(echo "$SAMPLE" | tr -d '[:space:]')"; [ -z "$SAMPLE" ] && continue
    RUNS="$(echo "$RUNS" | tr -d '[:space:]')"
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
