#!/bin/bash
# =============================================================================
# cNMF Workflow - Phase 1: Factorize and Generate K Selection Plot
# =============================================================================
# Runs on cleaned count matrices exported from Figure3.Rmd cNMF section.
# Submits prep -> factorize -> combine -> k_select as a SLURM job chain
# for each sample sequentially.
#
# Usage:
#   ./cNMF_subpattern.sh [sample ...]        (default: every row in the CSV)
#
# After this completes, inspect the K selection plot at:
#   {CNMF_OUTPUT_DIR}/{CNMF_RUN_NAME}/{CNMF_RUN_NAME}.k_selection.png
#
# Then run Phase 2:
#   ./cNMF_postprocess.sh <K> [density]
#
# Steps executed (per sample):
#   0. convert    - Convert cleaned TSV to h5ad
#   1. preprocess - Filter genes, ensure raw counts
#   2. prepare    - Initialize cNMF, select HVGs
#   3. factorize  - Parallel NMF (array job)
#   4. combine    - Merge worker results
#   5. k_select   - Generate K selection plot and stability metrics
# =============================================================================

set -e  # Exit on error

# =============================================================================
# DATASET CONFIGURATION
# =============================================================================

# Script locations are resolved from this script's own path so the repo can be
# cloned anywhere; nothing below reaches outside it.
WRAPPER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${WRAPPER_DIR}/../../pipelineScripts"

# Filesystem root for every input read and output written. Override to run
# against a different copy of the data:
#   DATA_ROOT=/my/data ./cNMF_subpattern.sh [sample ...]
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody/extractedData}"

SUBPATTERN_DIR="${CNMF_BASE_DIR:-${DATA_ROOT}/Figure1/singleCell/cNMF}/subpattern"
CNMF_INPUT_DIR="${SUBPATTERN_DIR}"
CNMF_OUTPUT_DIR="${SUBPATTERN_DIR}"

# Samples and their cleaned run names come from cNMF_subpattern.csv. Naming a
# sample on the command line runs just that one.
SAMPLES_CSV="${WRAPPER_DIR}/cNMF_subpattern.csv"
[ -f "$SAMPLES_CSV" ] || { echo "Error: not found: $SAMPLES_CSV"; exit 1; }

SAMPLES=()
RUN_NAMES=()
while IFS=, read -r csv_sample csv_run csv_k; do
    [ -z "$csv_sample" ] && continue
    case "$csv_sample" in \#*) continue ;; esac
    if [ $# -gt 0 ]; then
        match=0
        for want in "$@"; do [ "$want" = "$csv_sample" ] && match=1; done
        [ $match -eq 1 ] || continue
    fi
    SAMPLES+=("$csv_sample")
    RUN_NAMES+=("$csv_run")
done < <(grep -v '^#' "$SAMPLES_CSV" | tail -n +2)

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "Error: no samples selected. Rows in $SAMPLES_CSV:"
    grep -v '^#' "$SAMPLES_CSV" | tail -n +2 | cut -d, -f1 | sed 's/^/  /'
    exit 1
fi
echo "Samples: ${SAMPLES[*]}"

# cNMF parameters
export CNMF_K_RANGE="5,6,7,8,9"
export CNMF_N_ITER=200
export CNMF_N_HVG=2000
export CNMF_SEED=36
export CNMF_N_WORKERS=4

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

echo "=============================================="
echo "cNMF Workflow - Phase 1 (cleaned inputs)"
echo "=============================================="
echo "Time: $(date)"
echo ""

LAST_JOB=""

for idx in "${!SAMPLES[@]}"; do
  SAMPLE="${SAMPLES[$idx]}"
  RUN_NAME="${RUN_NAMES[$idx]}"
  INPUT_TSV="${CNMF_INPUT_DIR}/${SAMPLE}_counts_cleaned.tsv"
  INPUT_H5AD="${CNMF_INPUT_DIR}/${SAMPLE}_counts_cleaned.h5ad"

  echo ""
  echo "----------------------------------------------"
  echo "Sample: ${SAMPLE}  Run: ${RUN_NAME}"
  echo "  Input TSV:  ${INPUT_TSV}"
  echo "  Output dir: ${CNMF_OUTPUT_DIR}"
  echo "----------------------------------------------"

  # Build dependency flag if chaining after previous sample
  DEP_FLAG=""
  if [ -n "${LAST_JOB}" ]; then
    DEP_FLAG="--dependency=afterok:${LAST_JOB}"
  fi

  # ===========================================================================
  # Job 1: Convert TSV to h5ad + Preprocess + Prepare
  # ===========================================================================
  JOB_PREP=$(sbatch --parsable ${DEP_FLAG} <<PREP_SCRIPT
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -p ${PARTITION}
#SBATCH --job-name=cnPAc_prep
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 0:30:00
#SBATCH --output=${LOG_DIR}/cnmf_prep_${RUN_NAME}_%j.out
#SBATCH --error=${LOG_DIR}/cnmf_prep_${RUN_NAME}_%j.err

export CNMF_INPUT_H5AD="${INPUT_H5AD}"
export CNMF_OUTPUT_DIR="${CNMF_OUTPUT_DIR}"
export CNMF_RUN_NAME="${RUN_NAME}"
export CNMF_K_RANGE="${CNMF_K_RANGE}"
export CNMF_N_ITER="${CNMF_N_ITER}"
export CNMF_N_HVG="${CNMF_N_HVG}"
export CNMF_SEED="${CNMF_SEED}"
export CNMF_N_WORKERS="${CNMF_N_WORKERS}"

module purge
source activate ${CONDA_ENV}
cd ${SCRIPT_DIR}

echo "=== Convert TSV to h5ad ==="
python3 -c "
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

df = pd.read_csv('${INPUT_TSV}', sep='\t', index_col='cellID')
adata = sc.AnnData(X=sp.csr_matrix(df.values.astype(float)))
adata.obs_names = df.index.tolist()
adata.var_names = df.columns.tolist()
adata.write('${INPUT_H5AD}')
print(f'Converted: {adata.n_obs} cells x {adata.n_vars} genes')
"

echo ""
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
  echo "  Submitted prep job: ${JOB_PREP}"

  # ===========================================================================
  # Job 2: Factorize (array job)
  # ===========================================================================
  JOB_FACTOR=$(sbatch --parsable --dependency=afterok:${JOB_PREP} --array=0-$((CNMF_N_WORKERS-1)) <<FACTOR_SCRIPT
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -p ${PARTITION}
#SBATCH --job-name=cnPAc_fact
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=2G
#SBATCH -t 1:00:00
#SBATCH --output=${LOG_DIR}/cnmf_factor_${RUN_NAME}_%A_%a.out
#SBATCH --error=${LOG_DIR}/cnmf_factor_${RUN_NAME}_%A_%a.err

export CNMF_OUTPUT_DIR="${CNMF_OUTPUT_DIR}"
export CNMF_RUN_NAME="${RUN_NAME}"
export CNMF_N_WORKERS="${CNMF_N_WORKERS}"

module purge
source activate ${CONDA_ENV}
cd ${SCRIPT_DIR}

echo "Worker \${SLURM_ARRAY_TASK_ID} starting: \$(date)"
python3 cNMF.py factorize
echo "Worker \${SLURM_ARRAY_TASK_ID} complete: \$(date)"
FACTOR_SCRIPT
)
  echo "  Submitted factorize array job: ${JOB_FACTOR}"

  # ===========================================================================
  # Job 3: Combine + K selection plot
  # ===========================================================================
  JOB_KSELECT=$(sbatch --parsable --dependency=afterok:${JOB_FACTOR} <<KSELECT_SCRIPT
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -p ${PARTITION}
#SBATCH --job-name=cnPAc_ksel
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=2G
#SBATCH -t 0:15:00
#SBATCH --output=${LOG_DIR}/cnmf_kselect_${RUN_NAME}_%j.out
#SBATCH --error=${LOG_DIR}/cnmf_kselect_${RUN_NAME}_%j.err

export CNMF_OUTPUT_DIR="${CNMF_OUTPUT_DIR}"
export CNMF_RUN_NAME="${RUN_NAME}"
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
echo "Phase 1 complete for ${RUN_NAME}: \$(date)"
KSELECT_SCRIPT
)
  echo "  Submitted combine + k_select job: ${JOB_KSELECT}"

  LAST_JOB="${JOB_KSELECT}"
done

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=============================================="
echo "Phase 1 Job Chains Submitted"
echo "=============================================="
echo "  Samples: ${SAMPLES[*]}"
echo ""
echo "Monitor with: squeue -u \$USER"
echo ""
echo "When complete, inspect the K selection plots:"
for idx in "${!RUN_NAMES[@]}"; do
  echo "  ${CNMF_OUTPUT_DIR}/${RUN_NAMES[$idx]}/${RUN_NAMES[$idx]}.k_selection.png"
done
echo ""
echo "Then run Phase 2:"
echo "  ./cNMF_postprocess.sh <K> [density]"
echo "=============================================="
