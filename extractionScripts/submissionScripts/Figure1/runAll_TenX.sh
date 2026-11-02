#!/bin/bash
# runAll_TenX.sh  (Figure1 TenX)
# Chains stepTwo -> stepThree -> stepFour for the TenX arm, waiting for each
# stage's SLURM jobs to drain before submitting the next. Each stage reads the
# inputs the previous one wrote.
#
# There is no stepOne here: the TenX cellranger runs are the original cellranger counts.
# stepTwo checks that every BAM is present and aborts if one is missing.
#
#   bash runAll_TenX.sh            # from stepTwo (default)
#   bash runAll_TenX.sh stepThree  # from stepThree, reusing the alignment tables
#
# Run detached, e.g.:
#   screen -dmS tenx_figure1 bash runAll_TenX.sh
#   screen -r tenx_figure1
#
# Progress goes to runAll_TenX.log next to this script; per-job logs stay in logs/.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Everything this script prints goes to runAll_TenX.log beside it as well as to the
# terminal. A detached `screen -dmS` run has no terminal at all, which is the
# usual way this is launched, so without the tee the progress log is lost. The
# log is appended to and each run is headed by its start time, so a rerun does
# not discard the record of the previous one.
LOG_FILE="${SCRIPT_DIR}/runAll_TenX.log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "===== $(date '+%Y-%m-%d %H:%M:%S')  $(basename "$0") $*  (pid $$) ====="

# Filesystem root for every input read and output written below. Override to run
# the pipeline against a different copy of the data:
#   DATA_ROOT=/my/data bash <this script>
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

START_STEP="${1:-stepTwo}"
case "$START_STEP" in
    stepTwo) START_N=2 ;;
    stepThree) START_N=3 ;;
    stepFour) START_N=4 ;;
    *) echo "Usage: $0 [stepTwo|stepThree|stepFour]"; exit 1 ;;
esac

TenX_ROOT="${DATA_ROOT}/extractedData/Figure1/singleCell"
SAMPLES_CSV="${SCRIPT_DIR}/stepTwo_bamsort_TenX.csv"

say() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# Block until no queued or running job of this user has a name starting with $1.
wait_for_jobs() {
    local prefix="$1"
    local n empty=0
    while true; do
        n=$(squeue -u "$USER" -h -o "%j" | grep -c "^${prefix}")
        if [ "$n" -eq 0 ]; then
            empty=$((empty + 1))
            [ "$empty" -ge 2 ] && break
        else
            empty=0
            say "  waiting on $n ${prefix}* job(s)"
        fi
        sleep 120
    done
}

# Fail loudly rather than feeding a half-finished stage into the next one.
check_alignment() {
    local missing=0 sample n
    while IFS=, read -r sample _rest; do
        sample="$(echo "$sample" | tr -d '[:space:]\r')"
        [ -z "$sample" ] && continue
        n=$(ls "${TenX_ROOT}/bamsort/alignment/${sample}_bamsort_alignment_"*.csv 2>/dev/null | wc -l)
        say "  ${sample}: $n alignment tables"
        [ "$n" -eq 0 ] && missing=$((missing + 1))
    done < <(tail -n +2 "$SAMPLES_CSV")
    return $missing
}

say "===== starting from ${START_STEP} ====="

if [ "$START_N" -le 2 ]; then
    say "===== stepTwo: bamsort alignment + splice ====="
    bash stepTwo_bamsort_TenX.sh || { say "stepTwo submission failed"; exit 1; }
    sleep 30
    wait_for_jobs "tenxbam-"
    say "stepTwo jobs finished; checking alignment tables"
    if ! check_alignment; then
        say "ABORTING: a sample has no alignment tables (see logs/)"
        exit 1
    fi
    say "splice/combined: $(ls "${TenX_ROOT}/bamsort/splice/combined/" 2>/dev/null | tr '\n' ' ')"
fi

if [ "$START_N" -le 3 ]; then
    say "===== stepThree: per-UMI category tables ====="
    bash stepThree_process_TenX.sh || { say "stepThree submission failed"; exit 1; }
    sleep 30
    wait_for_jobs "tenxThree"
    say "stepThree jobs finished"
    say "reads/: $(ls "${TenX_ROOT}/bamsort/reads/" 2>/dev/null | tr '\n' ' ')"
fi

say "===== stepFour: detection thresholds ====="
srun -A b1042 -p genomics -N 1 -n 1 --mem=32G -t 00:30:00 \
     bash -c "module purge; module load R/4.4.0; Rscript ${SCRIPT_DIR}/../../pipelineScripts/bamsort/stepFour_minimums_TenX.R" \
     || { say "stepFour failed"; exit 1; }
say "stepFour finished"

say "reads/minimum/: $(ls "${TenX_ROOT}/bamsort/reads/minimum/" 2>/dev/null | tr '\n' ' ')"
say "===== TenX pipeline complete ====="
