#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=brie2
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=16G
#SBATCH -t 10:00:00
#SBATCH --output=logs/brie2_%j.txt
#SBATCH --error=logs/brie2_%j.err
#SBATCH --verbose

# Per-sample BRIE2 worker: brie-count -> brie-quant -> per-cell CSV export.
# Submitted by Figure4/brie2.sh; not run directly. One worker serves both arms
# (the macaque host targets and the mouse liver targets).
#
# Args: $1=sample_name $2=counts_dir $3=output_dir $4=gff_file $5=out_name
#       $6=export_csv (filename inside out_name) [$7..=event ids to export]
#
# PIPELINE_DIR is exported by the wrapper so the export script is found without
# an absolute path.

module purge
source activate brie2

SAMPLE=$1
SHORT_NAME=$5
BAM=${2}/${SAMPLE}/outs/possorted_genome_bam.bam
OUTDIR=${3}/${SHORT_NAME}
GFF=$4
EXPORT_CSV=$6
shift 6
EVENTS=("$@")

mkdir -p ${OUTDIR}

if [ -f "${OUTDIR}/filtered_barcodes.tsv" ]; then
  BARCODES=${OUTDIR}/filtered_barcodes.tsv
  echo "Using pre-filtered barcodes for ${SAMPLE}: $(wc -l < ${BARCODES}) cells"
else
  BARCODES=${2}/${SAMPLE}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
  echo "Using Cell Ranger filtered barcodes for ${SAMPLE}"
fi

echo "Running brie-count for ${SAMPLE}..."
brie-count -a ${GFF} \
  -s ${BAM} \
  -b ${BARCODES} \
  -o ${OUTDIR} \
  -p 12

if [ ! -f "${OUTDIR}/brie_count.h5ad" ]; then
  echo "brie-count failed for ${SAMPLE}, skipping brie-quant."
  exit 1
fi

echo "Running brie-quant for ${SAMPLE}..."
brie-quant -i ${OUTDIR}/brie_count.h5ad \
  -o ${OUTDIR}/brie_quant.h5ad \
  -p 12

echo "Exporting per-cell BRIE2 isoform counts to ${EXPORT_CSV} for ${SAMPLE}..."
python "${PIPELINE_DIR}/isoform/export_brie_per_cell.py" \
  ${OUTDIR}/brie_count.h5ad \
  ${OUTDIR}/${EXPORT_CSV} \
  "${EVENTS[@]}"

if [ ! -f "${OUTDIR}/${EXPORT_CSV}" ]; then
  echo "WARNING: ${EXPORT_CSV} was not written for ${SAMPLE}."
fi

echo "Finished ${SAMPLE}."
