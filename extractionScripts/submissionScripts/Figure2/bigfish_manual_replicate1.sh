#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bigfishR1
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16G
#SBATCH -t 2:00:00
#SBATCH --output=logs/bigfish_manual_replicate1_%j.txt
#SBATCH --error=logs/bigfish_manual_replicate1_%j.err
#SBATCH --verbose

# Script locations are resolved from this script's own path so the repo can
# be cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

# bigfish_manual_replicate1.sh  (Figure2 HCR: EGS023, replicate 1)
# Runs bigfish_manual_pipeline.py over the manually segmented EGS023 FOVs and
# writes <OUTPUT_DIR>/all_cells_stats.csv, the per-cell US/SS/MS table Figure2.Rmd
# reads for the HCR arm. Replicate 2 is EGS027 (bigfish_manual_replicate2.sh).
#
# Submit from this directory: the SBATCH log paths are relative to logs/.
#
# The 58 masks live in one directory but their images are split across three date
# directories, so the pipeline runs once per date against the shared mask
# directory. Each pass discovers only its own date's FOVs; FOVs from the other
# dates have no image there and are simply not discovered.
#
# DATES omits 20251111 and 20251125. Both hold images but no FOV in either
# carries a mask, so a pass over them would process nothing.
#
# FOV_PATTERN captures everything between the date and the "_Camera - " optical
# configuration. EGS023 mixes three FOV namings -- bingo_slide<N>_<NNN>,
# bingo_shorter_slide6_001 and plate_<NNN> -- some with an _XY<N> suffix and some
# without, and this one lookahead covers all of them. It must capture the same
# substring from the image filenames and from the ROIs2Label_ mask filenames,
# since find_mask_for_fov() pairs them on an exact string equality of the capture.
#
# CHANNELS names the Gold channel `Gold70`, which is how this export writes it;
# EGS027 wrote `Gold`. The pipeline maps the three channels positionally onto pol,
# env and nef, so the order here is load-bearing.
#
# The pipeline appends to all_cells_stats.csv when the file already exists, which
# is what accumulates the three dates into one table. The preflight below refuses
# to start if that file is already present, so a resubmission cannot silently
# double the table.

module purge
source activate bigfish_env
# source activate sets CONDA_PREFIX but the spack python earlier in PATH still wins,
# so `python` resolves outside the env. Put the env's bin first, then confirm the
# import: without this check the run reaches the pipeline and fails per date with a
# bare ModuleNotFoundError.
export PATH="${CONDA_PREFIX}/bin:${PATH}"
if ! python -c "import bigfish" 2>/dev/null; then
    echo "ERROR: bigfish_env is not active (python=$(which python))"
    exit 1
fi

# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data bash bigfish_manual_replicate1.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

BASE_INPUT="${DATA_ROOT}/rawData/EGS023"
# Detection-threshold mode, chosen by the first argument:
#   perfov (default) lets every FOV set its own threshold, the bigfish default
#   batch            uses one threshold per channel per date, the mean of that
#                    date's per-FOV automatic thresholds
# Each DATE is a separate invocation, so "batch" is per date, not per run. The
# two modes write to separate directories so they can be compared.
MODE="${1:-perfov}"
case "${MODE}" in
    perfov) OUTPUT_DIR="${DATA_ROOT}/extractedData/Figure2/bigFISH/manual/"; THRESHOLD_ARG="" ;;
    batch)  OUTPUT_DIR="${DATA_ROOT}/extractedData/Figure2/bigFISH/manual_batchthr/"; THRESHOLD_ARG="--batch_threshold" ;;
    *) echo "ERROR: unknown mode '${MODE}', expected perfov or batch"; exit 1 ;;
esac
MASKS_DIR="${DATA_ROOT}/rawData/EGS023/masks/"
SCRIPT="${PIPELINE_DIR}/fish/bigfish_manual_pipeline.py"
CHANNELS="Red Gold70 FarRed"
FOV_PATTERN='(\d{8}_[A-Z0-9]+_.+?)(?=_Camera - )'
DATES=("20251117" "20251118" "20251119")

COMBINED_CSV="${OUTPUT_DIR%/}/all_cells_stats.csv"
if [ -f "${COMBINED_CSV}" ]; then
    echo "ERROR: ${COMBINED_CSV} already exists and the pipeline appends to it."
    echo "       Move or delete it before rerunning, or results will be doubled."
    exit 1
fi

# Preflight every date before running any of them. A channel named wrong, a FOV
# with no mask and a FOV whose analysis channels do not cover every z-slice all
# resolve to a warn-and-skip inside process_fov(), so the pipeline can otherwise
# exit 0 having written nothing.
for DATE in "${DATES[@]}"; do
    INPUT_DIR="${BASE_INPUT}/${DATE}/TIFF/"
    if [ ! -d "${INPUT_DIR}" ]; then
        echo "ERROR: input directory does not exist: ${INPUT_DIR}"
        exit 1
    fi
    echo "--- preflight ${DATE} ---"
    python - "${INPUT_DIR}" "${MASKS_DIR}" "${FOV_PATTERN}" ${CHANNELS} <<'PYEOF'
import os, re, sys, collections
input_dir, masks_dir, fov_pattern = sys.argv[1], sys.argv[2], sys.argv[3]
channels = sys.argv[4:]

fov_re = re.compile(fov_pattern)
z_re = re.compile(r'(Z\d+)')
channel_re = re.compile(r'_Camera - (\S+?)_Z\d+\.tif$')

per_fov = collections.defaultdict(lambda: collections.defaultdict(set))
for filename in os.listdir(input_dir):
    if not filename.endswith('.tif'):
        continue
    fov_match, channel_match, z_match = (fov_re.search(filename),
                                         channel_re.search(filename),
                                         z_re.search(filename))
    if fov_match and channel_match and z_match:
        per_fov[fov_match.group(1)][channel_match.group(1)].add(z_match.group(1))

# find_mask_for_fov() only considers files with the ROIs2Label_ prefix.
mask_fovs = set()
for filename in os.listdir(masks_dir):
    if filename.startswith('ROIs2Label_') and filename.endswith('.tif'):
        match = fov_re.search(filename)
        if match:
            mask_fovs.add(match.group(1))

image_fovs = set(per_fov)
if not image_fovs:
    print("ERROR: FOV pattern matched no image filenames in {}".format(input_dir))
    sys.exit(1)
if not mask_fovs:
    print("ERROR: no ROIs2Label_*.tif masks in {}".format(masks_dir))
    sys.exit(1)

to_process = image_fovs & mask_fovs
if not to_process:
    print("ERROR: no image FOV in this date matches a mask FOV")
    sys.exit(1)

# discover_z_slices() unions the z-slices of every channel present for a FOV,
# including channels this run does not analyse, and load_and_stack_fov() then
# raises FileNotFoundError if an analysis channel is missing any of them.
missing_channel, ragged_z = [], []
for fov in sorted(to_process):
    z_union = set().union(*per_fov[fov].values())
    for channel in channels:
        if channel not in per_fov[fov]:
            missing_channel.append((fov, channel))
        elif z_union - per_fov[fov][channel]:
            ragged_z.append((fov, channel, sorted(z_union - per_fov[fov][channel])))

# discover_z_slices() and load_and_stack_fov() select files by `fov in filename`,
# so one FOV id being a prefix of another would pull in the wrong files.
collisions = [(a, b) for a in to_process for b in image_fovs if a != b and a in b]

print("Preflight: {} image FOVs, {} of them masked".format(
    len(image_fovs), len(to_process)))
print("  channels present: {}".format(
    ", ".join(sorted({c for fov in per_fov for c in per_fov[fov]}))))

for fov, channel in missing_channel:
    print("  ERROR: {} has no {} channel".format(fov, channel))
for fov, channel, missing in ragged_z:
    print("  ERROR: {} channel {} missing z-slices {}".format(fov, channel, missing))
for a, b in collisions:
    print("  ERROR: FOV id {} is a substring of {}".format(a, b))

if missing_channel or ragged_z or collisions:
    sys.exit(1)
PYEOF
    if [ $? -ne 0 ]; then
        echo "ERROR: preflight failed for ${DATE}, nothing has been run"
        exit 1
    fi
done

mkdir -p "${OUTPUT_DIR}"

EXIT_CODE=0
for DATE in "${DATES[@]}"; do
    INPUT_DIR="${BASE_INPUT}/${DATE}/TIFF/"
    echo "========================================"
    echo "EGS023 (replicate 1) manual big-FISH: ${DATE}"
    echo "  input:  ${INPUT_DIR}"
    echo "  masks:  ${MASKS_DIR}"
    echo "  output: ${OUTPUT_DIR}"
    echo "  mode:   ${MODE}"
    echo "========================================"

    python "${SCRIPT}" \
        --input_dir "${INPUT_DIR}" \
        --output_dir "${OUTPUT_DIR}" \
        --mask_dir "${MASKS_DIR}" \
        --channels ${CHANNELS} \
        --fov_pattern "${FOV_PATTERN}" \
        ${THRESHOLD_ARG}

    DATE_CODE=$?
    echo "Finished ${DATE} with exit code ${DATE_CODE}"
    [ ${DATE_CODE} -ne 0 ] && EXIT_CODE=${DATE_CODE}
done

echo "Finished EGS023 replicate 1 with exit code ${EXIT_CODE}"
exit ${EXIT_CODE}
