#!/bin/bash
# stepFive_markers.sh  (Figure1 SIV samples)
# Submits pipelineScripts/bamsort/stepFive_markers.R, which writes the per-cell marker
# tables the CD4 T gate reads. One job covers every sample; the cost is loading
# each Seurat object once, so it is memory-bound rather than CPU-bound.
#
# Run after stepFour_minimums*.R, though it only need the Seurat objects.

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can be
# cloned anywhere
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

# Filesystem root for every input read and output written below. Override to run
# the pipeline against a different copy of the data:
#   DATA_ROOT=/my/data bash <this script>
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"
PLOT_DIR="${SCRIPT_DIR}/../../../plotScripts"

R_SCRIPT="${PIPELINE_DIR}/bamsort/stepFive_markers.R"

# The gate code (cd4t_gate_genes, cd4t_call, marker_counts) lives in the
# plotScripts functions.R source file
export FUNCTIONS_R="${PLOT_DIR}/functions.R"
export RDS_DIR="${DATA_ROOT}/extractedData/Figure1/singleCell/Seurat/RDS/"
# Must match `cd4t_expression_dir` in functions.R, which is where the .Rmd
# figures read these tables from.
export OUTPUT_DIR="${DATA_ROOT}/extractedData/Figure1/singleCell/Seurat/expressionDF/"
export MARKER_SAMPLES="D0|D195|Invitro|Pacute|W0|W2"

[ -f "$R_SCRIPT" ] || { echo "Error: file not found: $R_SCRIPT"; exit 1; }
[ -f "$FUNCTIONS_R" ] || { echo "Error: file not found: $FUNCTIONS_R"; exit 1; }

mkdir -p logs

sbatch \
    --job-name="markers" \
    --account=b1042 \
    --partition=genomics \
    --time=02:00:00 \
    --mem=64G \
    --nodes=1 \
    --ntasks=1 \
    --output="logs/stepFive_markers_%j.txt" \
    --error="logs/stepFive_markers_%j.err" \
    --wrap="Rscript ${R_SCRIPT}"

echo "stepFive markers job submitted."
echo "Output: ${OUTPUT_DIR}"
