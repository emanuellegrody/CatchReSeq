#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --time=08:00:00
#SBATCH --job-name=Giron_dl
#SBATCH --array=1-8%4
#SBATCH --output=/projects/b1042/GoyalLab/egrody/publicData/Gironetal2024/logs/dl_%A_%a.out
#SBATCH --error=/projects/b1042/GoyalLab/egrody/publicData/Gironetal2024/logs/dl_%A_%a.err

# Download Giron et al. 2024 (SRA SRP486539 / BioProject PRJNA1070656): eight
# paired-end Illumina runs (HIV-A/B/C/D). The already-gzipped fastqs are pulled
# straight from ENA and each file is checked against ENA's published md5, so a
# truncated transfer fails loudly instead of silently leaving a short fastq.
# One array task per run; download_manifest.tsv holds one row per file:
#   run  sample  url  md5  bytes
set -euo pipefail

DEST="${DEST:-/projects/b1042/GoyalLab/egrody/publicData/Gironetal2024/fastq}"
MANIFEST="${MANIFEST:-/projects/b1042/GoyalLab/egrody/publicData/Gironetal2024/scripts/download_manifest.tsv}"
mkdir -p "${DEST}"
cd "${DEST}"

# Runs in stable manifest order; this task handles the index-th one.
mapfile -t RUNS < <(cut -f1 "${MANIFEST}" | awk '!seen[$0]++')
RUN=${RUNS[$((SLURM_ARRAY_TASK_ID-1))]:-}
[ -n "${RUN}" ] || { echo "no run at index ${SLURM_ARRAY_TASK_ID}"; exit 1; }
echo "==================== ${RUN} (task ${SLURM_ARRAY_TASK_ID}) ===================="

verify() {  # $1=file $2=expected md5 ; true if present and md5 matches
  [ -s "$1" ] || return 1
  [ "$(md5sum "$1" | cut -d' ' -f1)" = "$2" ]
}

while IFS=$'\t' read -r run sample url md5 bytes; do
  [ "$run" = "$RUN" ] || continue
  fname=$(basename "$url")
  for attempt in 1 2 3; do
    if verify "${DEST}/${fname}" "$md5"; then
      echo "  ${fname}: md5 OK (${md5})"; break
    fi
    echo "[$(date)] ${fname}: download attempt ${attempt} (expect ${bytes} bytes)"
    curl -fSL --retry 5 --retry-delay 10 -o "${DEST}/${fname}" "https://${url}" || true
    if verify "${DEST}/${fname}" "$md5"; then
      echo "  ${fname}: md5 OK (${md5})"; break
    fi
    echo "  ${fname}: md5 mismatch (got $(md5sum "${DEST}/${fname}" 2>/dev/null | cut -d' ' -f1 || echo none), want ${md5})"
    if [ "$attempt" -eq 3 ]; then echo "ERROR: ${fname} failed md5 after 3 attempts"; exit 1; fi
    rm -f "${DEST}/${fname}"
  done
done < "${MANIFEST}"

echo "[$(date)] ${RUN} done:"
awk -F'\t' -v r="$RUN" '$1==r{print $3}' "${MANIFEST}" | while read -r u; do ls -lh "${DEST}/$(basename "$u")"; done
