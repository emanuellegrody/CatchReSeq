#!/bin/bash
# memoryseq_qc.sh  (Figure4 / LD002 MemorySeq clone panel)
# Submits pipelineScripts/memoryseq/memoryseq_qc.R, which checks the LD002 TPM
# matrix for the between-library artefacts the heritable-gene call would read as
# heritable variation. No figure panel depends on it; run it once after
# memoryseq_nfcore_rnaseq.sh and before trusting the hit list.

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can
# be cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

R_SCRIPT="${PIPELINE_DIR}/memoryseq/memoryseq_qc.R"

# Filesystem roots for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data PLOT_ROOT=/my/plots bash memoryseq_qc.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"
PLOT_ROOT="${PLOT_ROOT:-${DATA_ROOT}/plotData}"

export TPM_MATRIX="${DATA_ROOT}/extractedData/LD002/nf-core/star_salmon/salmon.merged.gene_tpm.tsv"
export PLOT_DIR="${PLOT_ROOT}/Figure4/memoryseq_qc/"
export SAMPLE_PATTERN="^Clone_"

[ -f "$R_SCRIPT" ] || { echo "Error: file not found: $R_SCRIPT"; exit 1; }
[ -f "$TPM_MATRIX" ] || {
    echo "Error: TPM matrix not found: $TPM_MATRIX"
    echo "       Run memoryseq_nfcore_rnaseq.sh first."
    exit 1
}

mkdir -p logs

sbatch \
    --job-name="ld002qc" \
    --account=b1042 \
    --partition=genomics \
    --time=00:20:00 \
    --mem=8G \
    --nodes=1 \
    --ntasks=1 \
    --output="logs/memoryseq_qc_%j.txt" \
    --error="logs/memoryseq_qc_%j.err" \
    --wrap="module purge && module load R/4.4.0 && Rscript ${R_SCRIPT}"

echo "memoryseq QC job submitted."
echo "Output: ${PLOT_DIR}"
