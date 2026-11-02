#!/bin/bash
# stepFive_markers.sh  (Figure4 / EGS024 barcoded in vitro)
# Submits pipelineScripts/bamsort/stepFive_markers.R, which writes the per-cell marker
# table the CD4 T gate reads for GFP_invitro, alongside the legacy
# GFP_invitro_CD4T.csv in the Figure4 expressionDF tree. That legacy table is
# log-normalized and carries only 7 of the 17 exclusion markers, which makes it
# a CD3-only gate. The cost is loading the Seurat object once, so the job is
# memory-bound rather than CPU-bound; 64G covers the 1.2 GB object, which is the
# size the Figure1 SIV objects run at.
#
# Depends only on the Seurat object, so it can run at any point after
# GFP_invitro.rds exists. It is numbered last because it is the step the .Rmd
# figures need in place before they can gate.

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can
# be cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"
PLOT_DIR="${SCRIPT_DIR}/../../../plotScripts"

R_SCRIPT="${PIPELINE_DIR}/bamsort/stepFive_markers.R"

# The gate code (cd4t_gate_genes, cd4t_call, marker_counts) lives in the
# plotScripts functions.R -- the analysis library the .Rmd figures source -- not
# in the smaller pipelineScripts one every other step here uses.
export FUNCTIONS_R="${PLOT_DIR}/functions.R"
# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data bash stepFive_markers.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

export RDS_DIR="${DATA_ROOT}/extractedData/Figure4/singleCell/Seurat/RDS/"
# NOT `cd4t_expression_dir`, which points at the Figure1 tree. Figure4 call
# sites must name this directory:
#   apply_cd4t_gate(df, "GFP_invitro", viral_cols = ..., expression.dir = <this>)
#   cd4t_gate_table("GFP_invitro", expression.dir = <this>)
export OUTPUT_DIR="${DATA_ROOT}/extractedData/Figure4/singleCell/Seurat/expressionDF/"
export MARKER_SAMPLES="GFP_invitro"

[ -f "$R_SCRIPT" ] || { echo "Error: file not found: $R_SCRIPT"; exit 1; }
[ -f "$FUNCTIONS_R" ] || { echo "Error: file not found: $FUNCTIONS_R"; exit 1; }

mkdir -p logs

sbatch \
    --job-name="f4markers" \
    --account=b1042 \
    --partition=genomics \
    --time=01:00:00 \
    --mem=64G \
    --nodes=1 \
    --ntasks=1 \
    --output="logs/stepFive_markers_%j.txt" \
    --error="logs/stepFive_markers_%j.err" \
    --wrap="module purge && module load R/4.4.0 && Rscript ${R_SCRIPT}"

echo "stepFive markers job submitted."
echo "Output: ${OUTPUT_DIR}GFP_invitro_CD4T_counts.csv"
