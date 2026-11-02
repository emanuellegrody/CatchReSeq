#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=120G
#SBATCH -t 4:00:00
#SBATCH --job-name=markpan
#SBATCH --output=suppfig3_marker_panels_%j.txt
#SBATCH --error=suppfig3_marker_panels_%j.err
set -uo pipefail

# Script locations are resolved from this script's own path so the repo can be
# cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"
PLOT_DIR="${SCRIPT_DIR}/../../../plotScripts"

# Filesystem roots for the inputs read and outputs written. Override to run
# against a different copy of the data:
#   DATA_ROOT=/my/data PLOT_ROOT=/my/plots sbatch suppfig3_marker_panels.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody/extractedData}"
PLOT_ROOT="${PLOT_ROOT:-/projects/b1042/GoyalLab/egrody/plotData}"

# The gate and marker helpers live in the plotScripts functions.R, the analysis
# library the .Rmd figures source, not the pipelineScripts one.
export FUNCTIONS_R="${PLOT_DIR}/functions.R"
export RDS_DIR="${DATA_ROOT}/Figure1/singleCell/Seurat/RDS/"
export TABLE_DIR="${DATA_ROOT}/Figure1/combined/tiers/"
export OUTPUT_DIR="${PLOT_ROOT}/Supplemental/SuppFig3/"

# msigdbr's curl dependency needs libssl.so.10, absent from the compute nodes
export LD_LIBRARY_PATH=/hpc/software/openssl/1.0.2k/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}

Rscript "${PIPELINE_DIR}/publicData/suppfig3_marker_panels.R"
