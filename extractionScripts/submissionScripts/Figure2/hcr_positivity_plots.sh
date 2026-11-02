#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=hcrpos
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH -t 0:30:00
#SBATCH --output=logs/hcr_positivity_plots_%j.txt
#SBATCH --error=logs/hcr_positivity_plots_%j.err
#SBATCH --verbose

# Script locations are resolved from this script's own path so the repo can
# be cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLOT_DIR="${SCRIPT_DIR}/../../../plotScripts"

# hcr_positivity_plots.sh  (Figure2 HCR positivity rule)
# Draws the per-FOV positivity panel and the per-FOV image overlays for the rule
# Figure2.Rmd applies to the HCR reference: n_US + n_SS > 0.
#
# Submit from this directory: the SBATCH log paths are relative to logs/.
# Requires bigfish_manual_replicate1.sh and bigfish_manual_replicate2.sh to have
# written all_cells_stats.csv under both batchthr directories.

module purge
source activate bigfish_env
# source activate sets CONDA_PREFIX but the spack python earlier in PATH still wins.
export PATH="${CONDA_PREFIX}/bin:${PATH}"
if ! python -c "import tifffile, skimage" 2>/dev/null; then
    echo "ERROR: bigfish_env is not active (python=$(which python))"
    exit 1
fi

PLOT_SCRIPTS="${PLOT_DIR}"
# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data bash hcr_positivity_plots.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

OUTPUT_DIR="${DATA_ROOT}/plotData/Figure2/positivityRule"
TABLE_DIR="${DATA_ROOT}/extractedData/Figure2/bigFISH"

for ARM in manual_batchthr manual_replicate2_batchthr; do
    if [ ! -f "${TABLE_DIR}/${ARM}/all_cells_stats.csv" ]; then
        echo "ERROR: missing ${TABLE_DIR}/${ARM}/all_cells_stats.csv"
        exit 1
    fi
done

mkdir -p "${OUTPUT_DIR}"

# FOVs for the image overlays, in three groups: fully masked fields (the least
# biased denominator), the largest old-vs-new disagreements, and well2, which
# carries a much higher spot load and acts as a positive reference.
FOV_LIST="${OUTPUT_DIR}/overlay_fovs.tsv"
cat > "${FOV_LIST}" <<'EOF'
EGS023	20251119_EGS023_plate_006_XY1	A_fullmask
EGS027	20260821_EGS027_well1_004_XY02	A_fullmask
EGS027	20260821_EGS027_well1_003_XY12	A_fullmask
EGS027	20260820_EGS027_well1_004_XY09	B_disagree
EGS027	20260820_EGS027_well1_005_XY09	B_disagree
EGS023	20251119_EGS023_plate_002_XY1	B_disagree
EGS027	20260821_EGS027_well2_001_XY01	C_highsignal
EGS027	20260821_EGS027_well2_001_XY19	C_highsignal
EGS027	20260821_EGS027_well2_001_XY34	C_highsignal
EOF

echo "=== per-FOV positivity panel ==="
python "${PLOT_SCRIPTS}/hcr_positivity_panel.py" \
    --table_dir "${TABLE_DIR}" \
    --output_dir "${OUTPUT_DIR}"
PANEL_CODE=$?

echo "=== per-FOV image overlays ==="
python "${PLOT_SCRIPTS}/hcr_positivity_overlay.py" \
    --fov_list "${FOV_LIST}" \
    --table_dir "${TABLE_DIR}" \
    --output_dir "${OUTPUT_DIR}"
OVERLAY_CODE=$?

EXIT_CODE=0
[ ${PANEL_CODE} -ne 0 ] && EXIT_CODE=${PANEL_CODE}
[ ${OVERLAY_CODE} -ne 0 ] && EXIT_CODE=${OVERLAY_CODE}
echo "Finished HCR positivity plots with exit code ${EXIT_CODE}"
exit ${EXIT_CODE}
