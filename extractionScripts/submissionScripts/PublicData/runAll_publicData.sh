#!/bin/bash
# runAll_publicData.sh
# One entry point for the published human datasets of Figure 3: Collora et al.
# 2022 (PRJNA777601), HIV-seq / Frouard et al. (GSE266455), and Giron et al.
# 2024 (GSE254483). Takes raw GEO/SRA/ENA data through to the inputs the
# Figure3.Rmd `# Public data` and `# Giron et al` sections read.
#
# Note this covers the PUBLIC data only. The in-house EGS026 re-sequencing of
# the Giron pools is a separate pipeline in submissionScripts/Figure3, and the
# SIV arm is Figure1/Figure2.
#
# Stages, in order:
#   reference   cellranger mkref of the joint GRCh38 + HXB2 genome
#   download    SRA/ENA fastqs for all three datasets, plus the published
#               count matrices the protein gates are scored on
#   cellranger  stepOne per dataset (5' v2 for Collora/HIV-seq, 3' v3
#               feature-barcode for Giron)
#   bamsort     stepTwo per dataset: HXB2 region tables + splice sites
#   adtgate     the two protein CD4 T gates, writing the CD4T_cellIDs.csv
#               files functions.R reads
#
# Each stage waits on the one before it through sbatch job dependencies, so the
# whole thing can be launched in one go and left alone. Stages are individually
# idempotent: anything already on disk is skipped, so a partial run can be
# resumed by re-running the same command.
#
# Usage:
#   bash runAll_publicData.sh                 # every stage
#   bash runAll_publicData.sh cellranger bamsort
#   bash runAll_publicData.sh --from bamsort  # that stage and all later ones
#   bash runAll_publicData.sh --dry-run       # print the plan, submit nothing
#
# Run from a login node: the two published-matrix downloads use wget directly
# rather than through SLURM, matching how download_GSE266455.sh was written.
#
# NOT a stage: the GEO metadata corrections. samples_manifest.csv for HIV-seq
# already carries the corrected sample names in SAMPLE, with GEO's own name
# kept in ORIGINAL_GEO_SAMPLE, so a run started from this script produces
# correctly labelled output directly. fix_PID8027_labels.py and
# fix_PID1052_labels.py in pipelineScripts repair output that was already built
# under GEO's names; they are not needed here.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUB="${SCRIPT_DIR}"          # the per-dataset dirs nest under this one
PIPE="${SCRIPT_DIR}/../../pipelineScripts"

# Everything this script prints goes to runAll_publicData.log beside it as well
# as to the terminal. A detached `screen -dmS` run has no terminal, which is the
# usual way this is launched, so without the tee the record of what was
# submitted is lost. The log is appended to and each run is headed by its start
# time and arguments, so a --dry-run and a real run are told apart and a rerun
# does not discard the previous record.
LOG_FILE="${SCRIPT_DIR}/runAll_publicData.log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "===== $(date '+%Y-%m-%d %H:%M:%S')  $(basename "$0") $*  (pid $$) ====="
# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data bash runAll_publicData.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

DATA="${DATA_ROOT}"
GENOMES="${DATA}/genomes"
PUBLIC="${DATA}/publicData"

ALL_STAGES=(reference download cellranger bamsort adtgate)
DRY_RUN=false

# --- argument parsing --------------------------------------------------------
STAGES=()
while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --from)
            [ -n "${2:-}" ] || { echo "--from needs a stage name"; exit 1; }
            found=false
            for s in "${ALL_STAGES[@]}"; do
                [ "$s" = "$2" ] && found=true
                $found && STAGES+=("$s")
            done
            $found || { echo "Unknown stage: $2"; exit 1; }
            shift 2 ;;
        -h|--help) sed -n '2,45p' "${BASH_SOURCE[0]}"; exit 0 ;;
        -*) echo "Unknown option: $1"; exit 1 ;;
        *)  STAGES+=("$1"); shift ;;
    esac
done
[ ${#STAGES[@]} -eq 0 ] && STAGES=("${ALL_STAGES[@]}")

wants() { for s in "${STAGES[@]}"; do [ "$s" = "$1" ] && return 0; done; return 1; }

# logs/ is gitignored, so a fresh clone has none, and the stepTwo scripts carry
# a relative #SBATCH --output=logs/... that sbatch will not create for them.
for d in HIV HIVseq Giron; do mkdir -p "${SUB}/${d}/logs"; done
mkdir -p "${SCRIPT_DIR}/logs"

# Array sizes come from the manifests, so a manifest that reads as empty would
# produce --array=1-0 and an sbatch error that does not say which file is at
# fault. Fail here instead, naming it.
require_count() {
    local n="$1" what="$2" src="$3"
    [ "${n:-0}" -gt 0 ] || { echo "Error: no ${what} found in ${src}"; exit 1; }
}

# --- submission helpers ------------------------------------------------------
# DEPS accumulates the job ids a stage produced; the next stage waits on them.
DEPS=""

dep_arg() { [ -n "$DEPS" ] && echo "--dependency=afterok:${DEPS}" || echo ""; }

# sbatch one job and return its id. Honours --dry-run.
submit() {
    if $DRY_RUN; then echo "DRY: sbatch $*" >&2; echo ""; return 0; fi
    local out; out="$(sbatch --parsable "$@")"
    echo "${out%%;*}"
}

add_dep() { [ -n "${1:-}" ] && DEPS="${DEPS:+${DEPS}:}$1"; }

# Run a stepOne script, which sbatches one job per sample itself. SBATCH_DEPENDENCY
# is read by sbatch as the default for --dependency, so exporting it here makes
# every job the script submits inherit this stage's dependency without the
# script needing to know about it.
run_stepone() {
    local dir="$1" label="$2"
    echo "--- cellranger: ${label} ---"
    if $DRY_RUN; then echo "DRY: (cd ${dir} && bash stepOne_cellranger.sh)"; return 0; fi
    local ids
    ids="$(cd "$dir" && SBATCH_DEPENDENCY="${DEPS:+afterok:${DEPS}}" \
           bash stepOne_cellranger.sh | tee /dev/stderr \
           | grep -oE 'Submitted batch job [0-9]+' | grep -oE '[0-9]+' || true)"
    for id in $ids; do NEW_DEPS="${NEW_DEPS:+${NEW_DEPS}:}$id"; done
}

banner() { echo; echo "=============== $1 ==============="; }

# =============================================================================
# reference
# =============================================================================
if wants reference; then
    banner "reference: joint GRCh38 + HXB2"
    if [ -d "${GENOMES}/GRCh38_HXB2" ]; then
        echo "already built: ${GENOMES}/GRCh38_HXB2 -- skipping"
    else
        id="$(cd "${SUB}/HIV" && submit --job-name=mkref mkref_GRCh38_HXB2.sh)"
        echo "mkref job: ${id:-<dry>}"
        add_dep "$id"
    fi
fi

# =============================================================================
# download
# =============================================================================
if wants download; then
    banner "download"

    # Published count matrices. These carry the CITE-seq protein the gates are
    # scored on, which our own cellranger runs do not reproduce, so they are a
    # pipeline input rather than a convenience. wget on the login node.
    if [ -d "${PUBLIC}/HIVseq/matrices" ] && \
       [ "$(ls -1 "${PUBLIC}/HIVseq/matrices" 2>/dev/null | wc -l)" -gt 0 ]; then
        echo "HIV-seq published matrices present -- skipping"
    elif $DRY_RUN; then
        echo "DRY: bash ${PIPE}/publicData/download_GSE266455.sh"
    else
        echo "--- fetching HIV-seq published matrices (GSE266455) ---"
        bash "${PIPE}/publicData/download_GSE266455.sh"
    fi

    COLLORA_TAR="${PUBLIC}/Colloraetal2022/processed/GSE187515_processed_data_files.tar.gz"
    if [ -d "${PUBLIC}/Colloraetal2022/processed/matrices" ]; then
        echo "Collora published matrices present -- skipping"
    elif $DRY_RUN; then
        echo "DRY: fetch + untar GSE187515_processed_data_files.tar.gz"
    else
        echo "--- fetching Collora published matrices (GSE187515) ---"
        mkdir -p "${PUBLIC}/Colloraetal2022/processed"
        [ -f "$COLLORA_TAR" ] || wget -q --tries=5 --waitretry=10 -O "$COLLORA_TAR" \
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE187nnn/GSE187515/suppl/GSE187515_processed_data_files.tar.gz"
        tar -xzf "$COLLORA_TAR" -C "${PUBLIC}/Colloraetal2022/processed"
    fi

    # FASTQs. Array sizes are derived rather than taken from the #SBATCH
    # directives, so adding a sample to a manifest does not silently leave its
    # runs undownloaded.
    NEW_DEPS=""

    N_COLLORA=$(tail -n +2 "${SUB}/HIV/samples_manifest.csv" | tr -d '\r' \
                | cut -d, -f3 | tr '|' '\n' | grep -c . || true)
    require_count "$N_COLLORA" "Collora SRA runs" "${SUB}/HIV/samples_manifest.csv"
    echo "--- Collora: ${N_COLLORA} SRA runs ---"
    id="$(submit $(dep_arg) --array="1-${N_COLLORA}%8" "${PIPE}/publicData/download_collora_sra.sh")"
    echo "  job ${id:-<dry>}"; [ -n "$id" ] && NEW_DEPS="${NEW_DEPS:+${NEW_DEPS}:}$id"

    N_HIVSEQ=$(grep -c . "${PUBLIC}/HIVseq/scripts/HIVseq_SRR.txt")
    require_count "$N_HIVSEQ" "HIV-seq SRA runs" "${PUBLIC}/HIVseq/scripts/HIVseq_SRR.txt"
    echo "--- HIV-seq: ${N_HIVSEQ} SRA runs ---"
    id_hs="$(submit $(dep_arg) --array="1-${N_HIVSEQ}%8" "${PIPE}/publicData/download_sra.sh")"
    echo "  job ${id_hs:-<dry>}"

    # The 16 runs whose barcode read SRA flags as technical: re-dumped with
    # --include-technical once the first pass has been attempted.
    N_REDL=$(grep -c . "${PUBLIC}/HIVseq/scripts/redownload_SRR.txt")
    require_count "$N_REDL" "HIV-seq re-dump runs" "${PUBLIC}/HIVseq/scripts/redownload_SRR.txt"
    echo "--- HIV-seq: ${N_REDL} runs re-dumped with --include-technical ---"
    redl_dep=""; [ -n "$id_hs" ] && redl_dep="--dependency=afterany:${id_hs}"
    id="$(submit $redl_dep --array="1-${N_REDL}%6" "${PIPE}/publicData/redownload_sra.sh")"
    echo "  job ${id:-<dry>}"; [ -n "$id" ] && NEW_DEPS="${NEW_DEPS:+${NEW_DEPS}:}$id"

    N_GIRON=$(cut -f1 "${PUBLIC}/Gironetal2024/scripts/download_manifest.tsv" | sort -u | grep -c . || true)
    require_count "$N_GIRON" "Giron ENA runs" "${PUBLIC}/Gironetal2024/scripts/download_manifest.tsv"
    echo "--- Giron: ${N_GIRON} ENA runs ---"
    id="$(submit $(dep_arg) --array="1-${N_GIRON}%4" "${PIPE}/publicData/download_ena.sh")"
    echo "  job ${id:-<dry>}"; [ -n "$id" ] && NEW_DEPS="${NEW_DEPS:+${NEW_DEPS}:}$id"

    DEPS="$NEW_DEPS"
fi

# =============================================================================
# cellranger
# =============================================================================
if wants cellranger; then
    banner "cellranger (stepOne)"
    NEW_DEPS=""
    run_stepone "${SUB}/HIV"    "Collora"
    run_stepone "${SUB}/HIVseq" "HIV-seq"
    run_stepone "${SUB}/Giron"  "Giron"
    DEPS="$NEW_DEPS"
fi

# =============================================================================
# bamsort
# =============================================================================
if wants bamsort; then
    banner "bamsort (stepTwo)"
    NEW_DEPS=""
    for d in HIV HIVseq Giron; do
        echo "--- bamsort: ${d} ---"
        id="$(cd "${SUB}/${d}" && submit $(dep_arg) stepTwo_bamsort.sh)"
        echo "  job ${id:-<dry>}"
        [ -n "$id" ] && NEW_DEPS="${NEW_DEPS:+${NEW_DEPS}:}$id"
    done
    DEPS="$NEW_DEPS"
fi

# =============================================================================
# adtgate
# =============================================================================
# Two gates rather than one because the panels differ: Collora stained isotype
# controls and gets the percentile rule that paper defines, HIV-seq did not and
# gets an Otsu cut with a depth guard. Each writes its dataset's
# CD4T_cellIDs.csv, which functions.R reads through adt_cd4t_paths.
if wants adtgate; then
    banner "adtgate"
    mkdir -p "${SCRIPT_DIR}/logs"
    for g in collora_adtgate adtgate_final; do
        echo "--- ${g}.R ---"
        id="$(submit $(dep_arg) \
                --job-name="${g}" \
                --account=b1042 --partition=genomics \
                --nodes=1 --ntasks=4 --mem=64G --time=4:00:00 \
                --output="${SCRIPT_DIR}/logs/${g}_%j.out" \
                --error="${SCRIPT_DIR}/logs/${g}_%j.err" \
                --wrap="module purge; module load R/4.4.0; Rscript ${PIPE}/${g}.R")"
        echo "  job ${id:-<dry>}"
    done
fi

banner "submitted"
cat <<EOF
Next: knit the Figure3.Rmd sections that read this output.
  # Public data   -> Load, Extraction, Integrate, then the Plotting chunks
  # Giron et al   -> Extraction, then Plotting
Each chunk's prerequisites are stated in the prose above it. The Seurat steps
are memory-heavy; run them under srun rather than on a login node.
EOF
