#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=16GB
#SBATCH --time=8:00:00
#SBATCH --job-name=HIVseq_sra
#SBATCH --array=1-35%8
# Relative to the directory sbatch is run from; create it first with
# `mkdir -p logs`, as sbatch will not create it and the job fails if it is absent.
#SBATCH --output=logs/download_sra_%A_%a.out
#SBATCH --error=logs/download_sra_%A_%a.err

# Download HIVseq public SRA accessions: prefetch -> fasterq-dump -> gzip.
# One array task per accession (8 concurrent). Idempotent: skips an accession
# whose gzipped FASTQ(s) already exist, so failed tasks can be re-run safely.
set -euo pipefail

module purge
module load sratoolkit/3.2.1

# Data root for this dataset. CARE_PUBLIC_ROOT names the publicData directory
# directly; DATA_ROOT names its parent, and is what runAll_publicData.sh sets.
PUBLIC="${CARE_PUBLIC_ROOT:-${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}/publicData}"
DEST="${DEST:-${PUBLIC}/HIVseq/fastq}"
ACC_LIST="${ACC_LIST:-${PUBLIC}/HIVseq/scripts/HIVseq_SRR.txt}"
mkdir -p "${DEST}"
cd "${DEST}"

ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ACC_LIST}")
if [ -z "${ACC}" ]; then echo "no accession at line ${SLURM_ARRAY_TASK_ID}"; exit 1; fi
echo "==================== ${ACC} (task ${SLURM_ARRAY_TASK_ID}) ===================="

# Idempotency: if any gzipped fastq for this accession is already present, skip.
if ls "${DEST}/${ACC}"*.fastq.gz >/dev/null 2>&1; then
  echo "${ACC}: gzipped fastq already present, skipping"
  exit 0
fi

# Use node-local scratch for prefetch/extraction temp, then deliver to DEST.
TMP=${SLURM_TMPDIR:-/tmp}/${ACC}_$$
mkdir -p "${TMP}"
trap 'rm -rf "${TMP}"' EXIT

echo "[$(date)] prefetch ${ACC}"
prefetch --max-size 100G --output-directory "${TMP}" "${ACC}"

echo "[$(date)] fasterq-dump ${ACC}"
fasterq-dump --split-files --threads ${SLURM_NTASKS} \
    --temp "${TMP}" --outdir "${DEST}" "${TMP}/${ACC}/${ACC}.sra"

echo "[$(date)] gzip ${ACC} fastqs"
for fq in "${DEST}/${ACC}"*.fastq; do
  [ -f "$fq" ] && gzip -f "$fq"
done

echo "[$(date)] ${ACC} done:"
ls -lh "${DEST}/${ACC}"*.fastq.gz
