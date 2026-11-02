#!/bin/bash
# =============================================================================
# cNMF Workflow - Phase 1: Factorize and Generate K Selection Plot
# =============================================================================
# Submits prep -> factorize -> combine -> k_select as a SLURM job chain.
# Stops after producing the K selection plot so the user can choose K.
#
# Usage:
#   ./cNMF_workflow.sh <sample>
#
# <sample> picks a row from cNMF_workflow.csv, which carries the input .h5ad,
# the run name and the K that was chosen for that sample. One script serves
# every sample; edit the CSV, not this file.
#
# After this completes, inspect the K selection plot at:
#   {CNMF_OUTPUT_DIR}/{CNMF_RUN_NAME}/{CNMF_RUN_NAME}.k_selection.png
#
# Then run Phase 2:
#   ./cNMF_postprocess.sh <K> [density]
#
# Steps executed:
#   1. preprocess  - Filter genes, ensure raw counts
#   2. prepare     - Initialize cNMF, select HVGs
#   3. factorize   - Parallel NMF (array job)
#   4. combine     - Merge worker results
#   5. k_select    - Generate K selection plot and stability metrics
# =============================================================================

set -e  # Exit on error


# Script locations are resolved from this script's own path so the repo can be
# cloned anywhere; nothing below reaches outside it.
WRAPPER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${WRAPPER_DIR}/../../pipelineScripts"

# Filesystem root for every input read and output written. Override to run
# against a different copy of the data:
#   DATA_ROOT=/my/data ./cNMF_workflow.sh <sample>
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody/extractedData}"

# =============================================================================
# DATASET CONFIGURATION - from cNMF_workflow.csv, picked by the <sample> argument
# =============================================================================
SAMPLE="${1:-}"
if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 <sample>"
    echo "Samples in cNMF_workflow.csv:"
    grep -v '^#' "${WRAPPER_DIR}/cNMF_workflow.csv" | tail -n +2 | cut -d, -f1 | sed 's/^/  /'
    exit 1
fi

SAMPLES_CSV="${WRAPPER_DIR}/cNMF_workflow.csv"
[ -f "$SAMPLES_CSV" ] || { echo "Error: not found: $SAMPLES_CSV"; exit 1; }

ROW=$(grep -v '^#' "$SAMPLES_CSV" | tail -n +2 | awk -F, -v s="$SAMPLE" '$1 == s')
[ -n "$ROW" ] || { echo "Error: sample '$SAMPLE' not in $SAMPLES_CSV"; exit 1; }

export CNMF_RUN_NAME="$(echo "$ROW" | cut -d, -f2)"
CNMF_H5AD_REL="$(echo "$ROW" | cut -d, -f3)"
CNMF_K="$(echo "$ROW" | cut -d, -f4)"

export CNMF_INPUT_H5AD="${DATA_ROOT}/${CNMF_H5AD_REL}"
export CNMF_OUTPUT_DIR="${CNMF_BASE_DIR:-${DATA_ROOT}/Figure1/singleCell/cNMF}"

if [ ! -f "$CNMF_INPUT_H5AD" ]; then
    echo "Error: input not found: $CNMF_INPUT_H5AD"
    echo "       Write it with pipelineScripts/cnmf/extract_h5ad.R first."
    exit 1
fi

# cNMF parameters
export CNMF_K_RANGE="5,6,7,8,9,10"  # Comma-separated K values to test
export CNMF_N_ITER=200              # NMF iterations per K
export CNMF_N_HVG=2000              # Number of highly variable genes
export CNMF_SEED=36                 # Random seed for reproducibility
export CNMF_N_WORKERS=4             # Parallel workers (must match array size)

# =============================================================================
# SLURM CONFIGURATION
# =============================================================================
ACCOUNT="b1042"
PARTITION="genomics"
SCRIPT_DIR="${PIPELINE_DIR}/cnmf"
LOG_DIR="${WRAPPER_DIR}/logs"
CONDA_ENV="jupyter"

# =============================================================================
# Setup
# =============================================================================
mkdir -p "${LOG_DIR}"
cd "${SCRIPT_DIR}"

echo "=============================================="
echo "cNMF Workflow - Phase 1"
echo "=============================================="
echo "Time: $(date)"
echo ""
echo "Dataset configuration:"
echo "  Input:       ${CNMF_INPUT_H5AD}"
echo "  Output dir:  ${CNMF_OUTPUT_DIR}"
echo "  Run name:    ${CNMF_RUN_NAME}"
echo "  Chosen K:    ${CNMF_K} (pass to cNMF_postprocess.sh)"
echo ""
echo "cNMF parameters:"
echo "  K range:     ${CNMF_K_RANGE}"
echo "  Iterations:  ${CNMF_N_ITER}"
echo "  HVGs:        ${CNMF_N_HVG}"
echo "  Workers:     ${CNMF_N_WORKERS}"
echo "=============================================="

# =============================================================================
# Job 1: Preprocess + Prepare (sequential, single job)
# =============================================================================
JOB_PREP=$(sbatch --parsable <<PREP_SCRIPT
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -p ${PARTITION}
#SBATCH --job-name=cnmf_prep
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 0:30:00
#SBATCH --output=${LOG_DIR}/cnmf_prep_${CNMF_RUN_NAME}_%j.out
#SBATCH --error=${LOG_DIR}/cnmf_prep_${CNMF_RUN_NAME}_%j.err

export CNMF_INPUT_H5AD="${CNMF_INPUT_H5AD}"
export CNMF_OUTPUT_DIR="${CNMF_OUTPUT_DIR}"
export CNMF_RUN_NAME="${CNMF_RUN_NAME}"
export CNMF_K_RANGE="${CNMF_K_RANGE}"
export CNMF_N_ITER="${CNMF_N_ITER}"
export CNMF_N_HVG="${CNMF_N_HVG}"
export CNMF_SEED="${CNMF_SEED}"
export CNMF_N_WORKERS="${CNMF_N_WORKERS}"

module purge
source activate ${CONDA_ENV}
cd ${SCRIPT_DIR}

echo "=== Configuration ==="
python3 cNMF.py config

echo ""
echo "=== Preprocess ==="
python3 cNMF.py preprocess

echo ""
echo "=== Prepare ==="
python3 cNMF.py prepare

echo ""
echo "Prep complete: \$(date)"
PREP_SCRIPT
)
echo "Submitted prep job: ${JOB_PREP}"

# =============================================================================
# Job 2: Factorize (array job, depends on prep)
# =============================================================================
JOB_FACTOR=$(sbatch --parsable --dependency=afterok:${JOB_PREP} --array=0-$((CNMF_N_WORKERS-1)) <<FACTOR_SCRIPT
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -p ${PARTITION}
#SBATCH --job-name=cnmf_fact
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=2G
#SBATCH -t 1:00:00
#SBATCH --output=${LOG_DIR}/cnmf_factor_${CNMF_RUN_NAME}_%A_%a.out
#SBATCH --error=${LOG_DIR}/cnmf_factor_${CNMF_RUN_NAME}_%A_%a.err

export CNMF_OUTPUT_DIR="${CNMF_OUTPUT_DIR}"
export CNMF_RUN_NAME="${CNMF_RUN_NAME}"
export CNMF_N_WORKERS="${CNMF_N_WORKERS}"

module purge
source activate ${CONDA_ENV}
cd ${SCRIPT_DIR}

echo "Worker \${SLURM_ARRAY_TASK_ID} starting: \$(date)"
python3 cNMF.py factorize
echo "Worker \${SLURM_ARRAY_TASK_ID} complete: \$(date)"
FACTOR_SCRIPT
)
echo "Submitted factorize array job: ${JOB_FACTOR}"

# =============================================================================
# Job 3: Combine + K selection plot (depends on all factorize tasks)
# =============================================================================
JOB_KSELECT=$(sbatch --parsable --dependency=afterok:${JOB_FACTOR} <<KSELECT_SCRIPT
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -p ${PARTITION}
#SBATCH --job-name=cnmf_ksel
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=2G
#SBATCH -t 0:15:00
#SBATCH --output=${LOG_DIR}/cnmf_kselect_${CNMF_RUN_NAME}_%j.out
#SBATCH --error=${LOG_DIR}/cnmf_kselect_${CNMF_RUN_NAME}_%j.err

export CNMF_OUTPUT_DIR="${CNMF_OUTPUT_DIR}"
export CNMF_RUN_NAME="${CNMF_RUN_NAME}"
export CNMF_K_RANGE="${CNMF_K_RANGE}"

module purge
source activate ${CONDA_ENV}
cd ${SCRIPT_DIR}

echo "=== Combine ==="
python3 cNMF.py combine

echo ""
echo "=== K Selection ==="
python3 cNMF.py k_select

echo ""
echo "Phase 1 complete: \$(date)"
KSELECT_SCRIPT
)
echo "Submitted combine + k_select job: ${JOB_KSELECT}"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=============================================="
echo "Phase 1 Job Chain Submitted"
echo "=============================================="
echo "  Prep job:      ${JOB_PREP}"
echo "  Factorize job: ${JOB_FACTOR} (array 0-$((CNMF_N_WORKERS-1)))"
echo "  K select job:  ${JOB_KSELECT}"
echo ""
echo "Monitor with: squeue -u \$USER"
echo ""
echo "When complete, inspect the K selection plot:"
echo "  ${CNMF_OUTPUT_DIR}/${CNMF_RUN_NAME}/${CNMF_RUN_NAME}.k_selection.png"
echo ""
echo "Then run Phase 2:"
echo "  ./cNMF_postprocess.sh <K> [density]"
echo "=============================================="
