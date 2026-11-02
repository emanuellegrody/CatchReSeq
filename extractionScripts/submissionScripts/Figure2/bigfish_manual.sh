#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bigfish
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=4G
#SBATCH -t 1:30:00
#SBATCH --output=logs/bigfish_manual_%j.txt
#SBATCH --error=logs/bigfish_manual_%j.err
#SBATCH --verbose

# Script locations are resolved from this script's own path so the repo can
# be cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

module purge
source activate bigfish_env

# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data bash bigfish_manual.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

BASE_INPUT="${DATA_ROOT}/rawData/EGS023/"
OUTPUT_DIR="${DATA_ROOT}/extractedData/Figure2/bigFISH/manual/"
MASKS_DIR="${DATA_ROOT}/rawData/EGS023/masks/"
SCRIPT="${PIPELINE_DIR}/fish/bigfish_manual_pipeline.py"
CHANNELS="Red Gold70 FarRed"

#update
DATES=("20251117" "20251118" "20251119")

mkdir -p "${OUTPUT_DIR}"

for DATE in "${DATES[@]}"; do
    INPUT_DIR="${BASE_INPUT}/${DATE}/TIFF/"

    echo "========================================"
    echo "Processing: ${INPUT_DIR}"
    echo "========================================"
    
    # Check if input directory exists
    if [ ! -d "${INPUT_DIR}" ]; then
        echo "WARNING: Directory does not exist, skipping: ${INPUT_DIR}"
        continue
    fi
    
    # Run the pipeline
    python "${SCRIPT}" \
        --input_dir "${INPUT_DIR}" \
        --output_dir "${OUTPUT_DIR}" \
        --mask_dir "${MASKS_DIR}" \
        --channels ${CHANNELS}
    
    echo "Finished processing: ${DATE}"
    echo ""
done

INPUT_DIR="${BASE_INPUT}/20251119/TIFF/"
FOV_PATTERN='(\d{8}_[A-Z0-9]+_plate+_\d+(?:_XY\d+)?)'
python "${SCRIPT}" \
        --input_dir "${INPUT_DIR}" \
        --output_dir "${OUTPUT_DIR}" \
        --mask_dir "${MASKS_DIR}" \
        --channels ${CHANNELS} \
        --fov_pattern "${FOV_PATTERN}"

# Capture exit code
EXIT_CODE=$?
exit ${EXIT_CODE}