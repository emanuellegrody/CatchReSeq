#!/bin/bash
# =============================================================================
# cNMF Workflow - Phase 2: Postprocess with User-Selected K
# =============================================================================
# Run this after Phase 1 completes and you have inspected the K selection plot.
#
# Usage:
#   ./cNMF_postprocess.sh <sample_name> <K> [density]
#
# Examples:
#   ./cNMF_postprocess.sh D195 7         # K=7, density=2.0 (default)
#   ./cNMF_postprocess.sh D195 8 0.1    # K=8, density=0.1
#   ./cNMF_postprocess.sh Pacute 6 0.1
#   CNMF_BASE_DIR=<dir>/subpattern ./cNMF_postprocess.sh Pacute_cleaned 6 0.1
#
# The K selection plot is at:
#   {CNMF_OUTPUT_DIR}/{sample_name}/{sample_name}.k_selection.png
# =============================================================================

set -e

# =============================================================================
# DATASET CONFIGURATION
# =============================================================================

# Script locations are resolved from this script's own path so the repo can be
# cloned anywhere; nothing below reaches outside it.
WRAPPER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${WRAPPER_DIR}/../../pipelineScripts"

# Filesystem root for every input read and output written. Override to run
# against a different copy of the data:
#   DATA_ROOT=/my/data ./cNMF_postprocess.sh <run_name> <K> [density]
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody/extractedData}"

BASE_OUTPUT_DIR="${CNMF_BASE_DIR:-${DATA_ROOT}/Figure1/singleCell/cNMF}"

# =============================================================================
# Parse arguments
# =============================================================================
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "ERROR: Sample name and K value required."
    echo ""
    echo "Usage: ./cNMF_postprocess.sh <sample_name> <K> [density]"
    echo ""
    echo "Examples:"
    echo "  ./cNMF_postprocess.sh D195 8 0.1"
    echo "  ./cNMF_postprocess.sh Pacute 7 0.05"
    exit 1
fi

export CNMF_RUN_NAME="$1"
export CNMF_OUTPUT_DIR="${BASE_OUTPUT_DIR}"
SELECTED_K=$2
DENSITY=${3:-2.0}

# =============================================================================
# SLURM CONFIGURATION
# =============================================================================
ACCOUNT="b1042"
PARTITION="genomics"
SCRIPT_DIR="${PIPELINE_DIR}/cnmf"
LOG_DIR="${WRAPPER_DIR}/logs"
CONDA_ENV="jupyter"

# =============================================================================
# Validate that Phase 1 completed
# =============================================================================
PICKLE_PATH="${CNMF_OUTPUT_DIR}/${CNMF_RUN_NAME}_cnmf_object.pkl"
if [ ! -f "${PICKLE_PATH}" ]; then
    echo "ERROR: cNMF object not found at ${PICKLE_PATH}"
    echo "Phase 1 may not have completed. Run cNMF_workflow.sh first."
    exit 1
fi

mkdir -p "${LOG_DIR}"

echo "=============================================="
echo "cNMF Workflow - Phase 2: Postprocess"
echo "=============================================="
echo "Time: $(date)"
echo "  Run name:    ${CNMF_RUN_NAME}"
echo "  Selected K:  ${SELECTED_K}"
echo "  Density:     ${DENSITY}"
echo "=============================================="

# =============================================================================
# Submit postprocess job
# =============================================================================
JOB_POST=$(sbatch --parsable <<POST_SCRIPT
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -p ${PARTITION}
#SBATCH --job-name=cnPA_post
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=2G
#SBATCH -t 0:15:00
#SBATCH --output=${LOG_DIR}/cnmf_post_${CNMF_RUN_NAME}_%j.out
#SBATCH --error=${LOG_DIR}/cnmf_post_${CNMF_RUN_NAME}_%j.err

export CNMF_OUTPUT_DIR="${CNMF_OUTPUT_DIR}"
export CNMF_RUN_NAME="${CNMF_RUN_NAME}"

module purge
source activate ${CONDA_ENV}
cd ${SCRIPT_DIR}

echo "=== Postprocess (K=${SELECTED_K}, density=${DENSITY}) ==="
python3 cNMF.py postprocess --k ${SELECTED_K} --density ${DENSITY}

echo ""
echo "=============================================="
echo "cNMF Workflow Complete: \$(date)"
echo "=============================================="
echo "Output files in: ${CNMF_OUTPUT_DIR}"
echo "  - ${CNMF_RUN_NAME}_usage_K${SELECTED_K}.csv"
echo "  - ${CNMF_RUN_NAME}_spectra_K${SELECTED_K}.csv"
echo "  - ${CNMF_RUN_NAME}_spectra_tpm_K${SELECTED_K}.csv"
echo "  - ${CNMF_RUN_NAME}_top_genes_K${SELECTED_K}.csv"
echo ""
echo "Next: Run cNMF_analysis.R"
POST_SCRIPT
)

echo "Submitted postprocess job: ${JOB_POST}"
echo "Monitor with: squeue -u \$USER"
echo "=============================================="