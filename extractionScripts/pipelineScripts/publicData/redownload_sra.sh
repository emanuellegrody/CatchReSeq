#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=24GB
#SBATCH --time=20:00:00
#SBATCH --job-name=HIVseq_redl
#SBATCH --array=1-16%6
#SBATCH --output=/projects/b1042/GoyalLab/egrody/publicData/HIVseq/logs/redl_%A_%a.out
#SBATCH --error=/projects/b1042/GoyalLab/egrody/publicData/HIVseq/logs/redl_%A_%a.err

# Re-download 16 HIVseq TenX runs whose barcode read was dropped by the first
# attempt. Root cause: in these runs the 10X cell-barcode read (~30 bp) is
# flagged TECHNICAL in SRA, and plain `fasterq-dump --split-files` discards
# technical reads -> only the cDNA read (_4) survived. Fix: --include-technical,
# then drop the two short sample-index reads (i7/i5, ~10 bp), leaving the
# barcode read and the cDNA read (their _N numbering, e.g. _3/_4, is what the
# cellranger worker auto-detects).
#
# The download and the post-processing are separate resumable stages: if the
# uncompressed .fastq are already on disk (download finished), the run skips
# straight to dropping index reads + gzip without pulling the data again.
set -euo pipefail

module purge
module load sratoolkit/3.2.1
module load pigz/2.7-gcc-12.3.0

DEST="${DEST:-/projects/b1042/GoyalLab/egrody/publicData/HIVseq/fastq}"
ACC_LIST="${ACC_LIST:-/projects/b1042/GoyalLab/egrody/publicData/HIVseq/scripts/redownload_SRR.txt}"
mkdir -p "${DEST}"
cd "${DEST}"

ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ACC_LIST}")
if [ -z "${ACC}" ]; then echo "no accession at line ${SLURM_ARRAY_TASK_ID}"; exit 1; fi
echo "==================== ${ACC} (task ${SLURM_ARRAY_TASK_ID}) ===================="

# Length of the first read in a fastq(.gz). Run inside a pipefail-free subshell:
# `head -2` closes the pipe early, so the decompressor gets SIGPIPE (exit 141),
# which under `set -o pipefail` + `set -e` would otherwise abort the whole job.
first_read_len() { ( set +o pipefail; zcat -f "$1" 2>/dev/null | head -2 | tail -1 | tr -d '\n' | wc -c ); }

# Idempotency: skip if a barcode-length read (20-40 bp) is already present gzipped.
for fq in "${DEST}/${ACC}"_*.fastq.gz; do
  [ -e "$fq" ] || continue
  L=$(first_read_len "$fq")
  if [ "$L" -ge 20 ] && [ "$L" -le 40 ]; then
    echo "${ACC}: barcode read already present (${fq}, ${L}bp), skipping"; exit 0
  fi
done

# Reuse an already-completed download: if the uncompressed .fastq are present we
# go straight to post-processing rather than re-downloading. Only fetch from SRA
# when they are absent.
shopt -s nullglob
existing_fastq=( "${DEST}/${ACC}"_*.fastq )
shopt -u nullglob
if [ "${#existing_fastq[@]}" -gt 0 ]; then
  echo "[$(date)] ${ACC}: ${#existing_fastq[@]} uncompressed .fastq already present, skipping download"
else
  rm -f "${DEST}/${ACC}"_*.fastq.gz
  TMP=${DEST}/_tmp_download/${ACC}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
  mkdir -p "${TMP}"
  trap 'rm -rf "${TMP}"' EXIT

  echo "[$(date)] prefetch ${ACC}"
  prefetch --max-size 200G --output-directory "${TMP}" "${ACC}"

  echo "[$(date)] fasterq-dump --include-technical ${ACC}"
  fasterq-dump --split-files --include-technical --threads ${SLURM_NTASKS} \
      --temp "${TMP}" --outdir "${DEST}" "${TMP}/${ACC}/${ACC}.sra"
fi

echo "[$(date)] drop short index reads (<20 bp), keep barcode + cDNA"
have_bc=0; have_cdna=0
for fq in "${DEST}/${ACC}"_*.fastq; do
  [ -f "$fq" ] || continue
  L=$(first_read_len "$fq")
  if   [ "$L" -lt 20 ]; then echo "  drop index read $(basename "$fq") (${L}bp)"; rm -f "$fq"
  elif [ "$L" -le 40 ]; then echo "  keep barcode read $(basename "$fq") (${L}bp)"; have_bc=1
  else echo "  keep cDNA read $(basename "$fq") (${L}bp)"; have_cdna=1
  fi
done
if [ "$have_bc" -ne 1 ] || [ "$have_cdna" -ne 1 ]; then
  echo "ERROR: ${ACC} did not yield both a barcode read and a cDNA read"; exit 1
fi

echo "[$(date)] gzip ${ACC} fastqs"
for fq in "${DEST}/${ACC}"_*.fastq; do [ -f "$fq" ] && pigz -f -p "${SLURM_NTASKS}" "$fq"; done

echo "[$(date)] ${ACC} done:"
ls -lh "${DEST}/${ACC}"_*.fastq.gz
