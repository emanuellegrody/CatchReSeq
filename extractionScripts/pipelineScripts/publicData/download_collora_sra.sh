#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=16GB
#SBATCH --time=8:00:00
#SBATCH --job-name=Coll_sra
#SBATCH --array=1-56%8
#SBATCH --output=/projects/b1042/GoyalLab/egrody/publicData/Colloraetal2022/scripts/logs/sra_%A_%a.out
#SBATCH --error=/projects/b1042/GoyalLab/egrody/publicData/Colloraetal2022/scripts/logs/sra_%A_%a.err

# Download the Collora et al. 2022 (PRJNA777601) TenX SRA runs:
# prefetch -> fasterq-dump -> gzip. One array task per run. Idempotent: a run
# whose gzipped FASTQs are already present is skipped, so failed tasks can be
# re-run safely.
#
# Accessions are read from the RUNS column of the HIV submission manifest
# rather than a separate list, so the runs downloaded are exactly the runs
# stepOne_cellranger.sh will align and the two cannot drift apart. The manifest
# holds the 18 unstimulated TenX samples (56 runs); antigen-stimulated, ADT and
# VDJ libraries are deliberately absent.
#
# Unlike the HIV-seq runs, these store the 10X barcode read as a biological
# read (26 bp or 150 bp), so plain --split-files keeps it and the
# --include-technical workaround redownload_sra.sh needs does not apply here.
#
# The array size must match the accession count. Submit through
# submissionScripts/PublicData/runAll_publicData.sh, which sizes it, or by hand:
#   N=$(tail -n +2 ../../submissionScripts/PublicData/HIV/samples_manifest.csv | tr -d '\r' \
#       | cut -d, -f3 | tr '|' '\n' | grep -c .)
#   sbatch --array=1-${N}%8 download_collora_sra.sh

set -euo pipefail

module purge
module load sratoolkit/3.2.1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${DEST:-/projects/b1042/GoyalLab/egrody/publicData/Colloraetal2022/fastq}"
MANIFEST="${MANIFEST:-${SCRIPT_DIR}/../../submissionScripts/PublicData/HIV/samples_manifest.csv}"
[ -f "${MANIFEST}" ] || { echo "Missing manifest: ${MANIFEST}"; exit 1; }
mkdir -p "${DEST}"
cd "${DEST}"

# One accession per line, in manifest order, so a task id maps to a fixed run.
ACC=$(tail -n +2 "${MANIFEST}" | tr -d '\r' | cut -d, -f3 | tr '|' '\n' \
      | grep . | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ -z "${ACC}" ]; then echo "no accession at index ${SLURM_ARRAY_TASK_ID}"; exit 1; fi
echo "==================== ${ACC} (task ${SLURM_ARRAY_TASK_ID}) ===================="

if ls "${DEST}/${ACC}"*.fastq.gz >/dev/null 2>&1; then
  echo "${ACC}: gzipped fastq already present, skipping"
  exit 0
fi

# Node-local scratch for prefetch/extraction temp, then deliver to DEST.
TMP=${SLURM_TMPDIR:-/tmp}/${ACC}_$$
mkdir -p "${TMP}"
trap 'rm -rf "${TMP}"' EXIT

echo "[$(date)] prefetch ${ACC}"
prefetch --max-size 100G --output-directory "${TMP}" "${ACC}"

echo "[$(date)] fasterq-dump ${ACC}"
fasterq-dump --split-files --threads "${SLURM_NTASKS}" \
    --temp "${TMP}" --outdir "${DEST}" "${TMP}/${ACC}/${ACC}.sra"

echo "[$(date)] gzip ${ACC} fastqs"
for fq in "${DEST}/${ACC}"*.fastq; do
  [ -f "$fq" ] && gzip -f "$fq"
done

echo "[$(date)] ${ACC} done:"
ls -lh "${DEST}/${ACC}"*.fastq.gz
