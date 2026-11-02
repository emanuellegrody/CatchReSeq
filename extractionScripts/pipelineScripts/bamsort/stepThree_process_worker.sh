#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=64G
#SBATCH -t 4:00:00
#SBATCH --verbose

# stepThree worker: R processing of the stepTwo alignment tables.
# Submitted by a <dataset>/stepThree_process.sh wrapper; not run directly.
# All configuration is provided via env vars exported by the wrapper.
#
# Runs one job over all samples at once (rather than one job per sample, as
# stepTwo does): process_bamsort takes the full sample list and each sample is
# a cheap read of the alignment CSVs, so the fan-out is not worth the queue
# time. Read10X is the memory-heavy part, which is why the request is 64G.
#
# R_SCRIPT selects the pipeline:
#   bamsort_process.R       viral targets (SIV mac239 / HIV HXB2)
#   bamsort_process_host.R  host targets
# Each documents the env vars it reads.

set -euo pipefail

REQUIRED_VARS=(
    R_SCRIPT FUNCTIONS_R
    RAW_CELLIDS_OUTPUT_DIR RAW_CELLIDS_DATASETS RAW_CELLIDS_SAMPLES
    RAW_CELLIDS_FOLDERS GENOME_PREFIX
    PROCESS_SAMPLES PROCESS_INPUT_DIRS PROCESS_OUTPUT_DIR
)
for v in "${REQUIRED_VARS[@]}"; do
    if [ -z "${!v:-}" ]; then
        echo "Error: required env var $v not set (use the stepThree_process.sh wrapper)"
        exit 1
    fi
done

for f in "$R_SCRIPT" "$FUNCTIONS_R"; do
    if [ ! -f "$f" ]; then
        echo "Error: required file not found: $f"
        exit 1
    fi
done

mkdir -p "$PROCESS_OUTPUT_DIR" "$RAW_CELLIDS_OUTPUT_DIR"

echo "========== stepThree: $(basename "$R_SCRIPT") =========="
echo "Samples:    ${PROCESS_SAMPLES//|/, }"
echo "Input dirs: ${PROCESS_INPUT_DIRS//|/, }"
echo "Output dir: $PROCESS_OUTPUT_DIR"

module purge
module load R/4.4.0
# The cluster's R/4.4.0 module ships no Seurat, so the packages come from a
# user library. Override R_LIBS_USER to point at your own:
#   R_LIBS_USER=/path/to/R/library bash stepThree_process.sh
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/x86_64-pc-linux-gnu-library/4.4/}"

Rscript "$R_SCRIPT"

echo "stepThree finished."
