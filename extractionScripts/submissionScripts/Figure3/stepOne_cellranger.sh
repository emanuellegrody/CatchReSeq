#!/bin/bash
# stepOne_cellranger.sh  (EGS026 re-sequenced 10X pools, Figure3)
# Submit one feature-barcode cellranger job per pool (10X 3' v3 vs joint human
# GRCh38 + HXB2), producing a matrix with genes + HXB2 + CITE-seq ADT + donor
# hashtags. The Gene Expression library comes from the EGS026 fastq_concat/tenx
# bcl-convert output, with both sequencing runs concatenated per sample; the
# Antibody Capture library is not part of EGS026 and is
# taken from the matching Giron et al. 2024 (GSE254483) SRA run for that pool.
# Reads samples_manifest.csv (SAMPLE,GSM,TenX_FASTQ_SAMPLE,FB_RUN).
# Edit paths below, then: bash stepOne_cellranger.sh  [SAMPLE ...]

set -euo pipefail

### UPDATE PATHS
# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data GENOME_ROOT=/my/genomes bash stepOne_cellranger.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"
GENOME_ROOT="${GENOME_ROOT:-/projects/b1042/GoyalLab/egrody/genomes}"

export TenX_FASTQ_DIR="${DATA_ROOT}/rawData/EGS026/fastq_concat/tenx"
export FB_FASTQ_SRC_DIR="${DATA_ROOT}/publicData/Gironetal2024/fastq"
export TRANSCRIPTOME="${GENOME_ROOT}/GRCh38_HXB2"
export FEATURE_REF="${DATA_ROOT}/publicData/Gironetal2024/feature_ref/feature_reference.csv"
export OUTPUT_PATH="${DATA_ROOT}/extractedData/Figure3/singleCell/counts"
export OUTS_ONLY="FALSE"

### DO NOT UPDATE BELOW
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${SCRIPT_DIR}/samples_manifest.csv"
WORKER="${SCRIPT_DIR}/../../pipelineScripts/cellranger/featurebarcode_cellranger_worker.sh"

[ -f "$MANIFEST" ] || { echo "Missing $MANIFEST"; exit 1; }
[ -f "$WORKER" ]   || { echo "Missing $WORKER"; exit 1; }
mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "$OUTPUT_PATH"

WANTED=" $* "
tail -n +2 "$MANIFEST" | while IFS=, read -r SAMPLE GSM TenX_FASTQ_SAMPLE FB_RUN; do
    SAMPLE="$(echo "$SAMPLE" | tr -d '[:space:]')"; [ -z "$SAMPLE" ] && continue
    TenX_FASTQ_SAMPLE="$(echo "$TenX_FASTQ_SAMPLE" | tr -d '[:space:]')"
    FB_RUN="$(echo "$FB_RUN" | tr -d '[:space:]')"
    if [ "$#" -gt 0 ] && [[ "$WANTED" != *" $SAMPLE "* ]]; then continue; fi
    export SAMPLE TenX_FASTQ_SAMPLE FB_RUN
    echo "Submitting cellranger: $SAMPLE ($GSM)  TenX=$TenX_FASTQ_SAMPLE  FB=$FB_RUN"
    sbatch \
        --job-name="cr-EGS026-${SAMPLE}" \
        --output="${SCRIPT_DIR}/logs/stepOne_${SAMPLE}_%j.txt" \
        --error="${SCRIPT_DIR}/logs/stepOne_${SAMPLE}_%j.err" \
        --export=ALL \
        "$WORKER"
done
echo "Done submitting cellranger jobs."
