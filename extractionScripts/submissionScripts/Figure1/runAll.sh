#!/bin/bash
# runAll.sh  (Figure1 Ca*Re)
# Chains stepOne -> stepTwo -> stepThree -> stepFour, waiting for each stage's
# SLURM jobs to drain before submitting the next. Each stage reads the inputs the
# previous one wrote.
#
# Can rerun starting at a certain step if desired, for example:
#   bash runAll.sh            # from stepOne (default)
#   bash runAll.sh stepTwo    # from stepTwo, reusing the existing count BAMs
#
# Run detached, e.g.:
#   screen -dmS care_figure1 bash runAll.sh
#   screen -r care_figure1
#
# Progress goes to runAll.log next to this script; per-job logs stay in logs/.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Everything this script prints goes to runAll.log beside it as well as to the
# terminal. A detached `screen -dmS` run has no terminal at all, which is the
# usual way this is launched, so without the tee the progress log is lost. The
# log is appended to and each run is headed by its start time, so a rerun does
# not discard the record of the previous one.
LOG_FILE="${SCRIPT_DIR}/runAll.log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "===== $(date '+%Y-%m-%d %H:%M:%S')  $(basename "$0") $*  (pid $$) ====="

# Filesystem root for every input read and output written below. Override to run
# the pipeline against a different copy of the data:
#   DATA_ROOT=/my/data bash <this script>
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

START_STEP="${1:-stepOne}"
case "$START_STEP" in
    stepOne) START_N=1 ;;
    stepTwo) START_N=2 ;;
    stepThree) START_N=3 ;;
    stepFour) START_N=4 ;;
    *) echo "Usage: $0 [stepOne|stepTwo|stepThree|stepFour]"; exit 1 ;;
esac

CARE_ROOT="${DATA_ROOT}/extractedData/Figure1/CaRe"
GROUP_NAMES="viremia ART invitro uninfected"

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
check_counts() {
    local missing=0 bam
    while IFS=, read -r _id name group _fq; do
        name="$(echo "$name" | tr -d '[:blank:]\r')"
        group="$(echo "$group" | tr -d '[:blank:]\r')"
        [ -z "$name" ] && continue
        bam="${CARE_ROOT}/${group}/counts/Mmul_10_mac239_${name}/outs/possorted_genome_bam.bam"
        if [ ! -f "$bam" ]; then
            say "  MISSING: $bam"
            missing=$((missing + 1))
        fi
    done < <(tail -n +2 stepOne_mkcounts.csv)
    return $missing
}

check_alignment() {
    local missing=0 group n
    for group in $GROUP_NAMES; do
        n=$(ls "${CARE_ROOT}/${group}/bamsort/alignment/"*_bamsort_alignment_*.csv 2>/dev/null | wc -l)
        say "  ${group}: $n alignment tables"
        [ "$n" -eq 0 ] && missing=$((missing + 1))
    done
    return $missing
}

say "===== starting from ${START_STEP} ====="

if [ "$START_N" -le 1 ]; then
    say "===== stepOne: cellranger count ====="
    bash stepOne_mkcounts.sh || { say "stepOne submission failed"; exit 1; }
    sleep 30
    wait_for_jobs "mkcounts-"
    say "stepOne jobs finished; checking BAMs"
    if ! check_counts; then
        say "ABORTING: some cellranger runs did not produce a BAM (see logs/)"
        exit 1
    fi
    say "all 30 count BAMs present"
else
    say "skipping stepOne; verifying the existing count BAMs"
    if ! check_counts; then
        say "ABORTING: cannot start at ${START_STEP} without every count BAM"
        exit 1
    fi
    say "all 30 count BAMs present"
fi

if [ "$START_N" -le 2 ]; then
    say "===== stepTwo: bamsort alignment + splice ====="
    bash stepTwo_bamsort.sh || { say "stepTwo submission failed"; exit 1; }
    sleep 30
    wait_for_jobs "bamsort-"
    say "stepTwo jobs finished; checking alignment tables"
    if ! check_alignment; then
        say "ABORTING: a group has no alignment tables (see logs/)"
        exit 1
    fi
fi

if [ "$START_N" -le 3 ]; then
    say "===== stepThree: per-UMI category tables ====="
    bash stepThree_process.sh || { say "stepThree submission failed"; exit 1; }
    sleep 30
    wait_for_jobs "stepThree-"
    say "stepThree jobs finished"

    for group in $GROUP_NAMES; do
        say "${group} reads/: $(ls "${CARE_ROOT}/${group}/bamsort/reads/" 2>/dev/null | tr '\n' ' ')"
    done
fi

say "===== stepFour: detection thresholds ====="
srun -A b1042 -p genomics -N 1 -n 1 --mem=32G -t 00:30:00 \
     bash -c "module purge; module load R/4.4.0; Rscript ${SCRIPT_DIR}/../../pipelineScripts/bamsort/stepFour_minimums_CaRe.R" \
     || { say "stepFour failed"; exit 1; }
say "stepFour finished"

for group in $GROUP_NAMES; do
    say "${group} reads/minimum/: $(ls "${CARE_ROOT}/${group}/bamsort/reads/minimum/" 2>/dev/null | tr '\n' ' ')"
done
say "===== pipeline complete ====="
