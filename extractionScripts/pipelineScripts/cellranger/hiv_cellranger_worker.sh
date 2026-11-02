#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH --verbose

# Per-SAMPLE cellranger worker for the public Collora et al. 2022 HIV data
# (10X 5' v2 TenX), aligning to the joint human GRCh38 + HXB2 reference so host
# and HIV are mapped in one pass. The resulting BAM (with CB/UB/xf tags, incl.
# on HIV reads mapping the HXB2 contig) and valid-cell list feed downstream.
#
# One GEO sample (GSM) = several SRA runs (lanes). All of a sample's runs are
# combined into ONE cellranger run, so cell calling sees the full library.
#
# SRA fasterq-dump names files SRRxxxx_1/_2.fastq.gz, which cellranger will not
# parse, so each run is symlinked to the 10X convention as its own lane:
#   run i, _1 (barcode read) -> ${SAMPLE}_S1_L00i_R1_001.fastq.gz
#   run i, _2 (cDNA read)    -> ${SAMPLE}_S1_L00i_R2_001.fastq.gz
# Both 26x98 and 150x150 TenX geometries work with --chemistry=SC5P-R2 (barcode
# is read from R1's first 26 bp; R2 is the cDNA).
#
# Required env vars (from wrapper):
#   SAMPLE             clean sample name (cellranger run-dir suffix, FASTQ prefix)
#   RUNS               '|'-separated SRA run accessions for this sample
#   FASTQ_SRC_DIR      dir containing ${RUN}_1.fastq.gz and _2.fastq.gz
#   TRANSCRIPTOME      cellranger reference (joint human+HIV GRCh38_HXB2)
#   OUTPUT_PATH        dir to cd into before cellranger count
# Optional:
#   CELLRANGER_BIN     default /home/egy2296/packages/cellranger-7.2.0/cellranger
#   OUTS_ONLY          'TRUE' to keep only outs/ after the run (default FALSE)

set -euo pipefail

for v in SAMPLE RUNS FASTQ_SRC_DIR TRANSCRIPTOME OUTPUT_PATH; do
    if [ -z "${!v:-}" ]; then echo "Error: required env var $v not set"; exit 1; fi
done

CELLRANGER_BIN="${CELLRANGER_BIN:-/home/egy2296/packages/cellranger-7.2.0/cellranger}"
genome_name=$(basename "$TRANSCRIPTOME")
RUN_DIR="${OUTPUT_PATH}/${genome_name}_${SAMPLE}"

# Idempotency: skip if this sample's cellranger run already finished.
if [ -f "${RUN_DIR}/outs/possorted_genome_bam.bam" ]; then
    echo "cellranger already done for $SAMPLE, skipping: ${RUN_DIR}/outs"
    exit 0
fi

source activate cellranger
mkdir -p "$OUTPUT_PATH"
cd "$OUTPUT_PATH"

# Stage 10X-named symlinks for cellranger, one lane per SRA run.
RUNS="${RUNS//[$'\r\n']/}"          # tolerate CRLF manifests
IFS='|' read -r -a RUN_ARR <<< "$RUNS"
FQ_DIR="${OUTPUT_PATH}/fastqs_${SAMPLE}"
rm -rf "$FQ_DIR"; mkdir -p "$FQ_DIR"
lane=0
for run in "${RUN_ARR[@]}"; do
    lane=$((lane+1))
    L=$(printf "L%03d" "$lane")
    # SRA fasterq-dump numbers the barcode/cDNA reads differently per run
    # depending on whether index reads were stored first: a run may split to
    # _1(barcode)/_2(cDNA) or _3(barcode)/_4(cDNA). Auto-detect the pair; the
    # barcode read is always the lower-numbered of the two.
    if   [ -f "${FASTQ_SRC_DIR}/${run}_1.fastq.gz" ]; then bc=1; cd=2
    elif [ -f "${FASTQ_SRC_DIR}/${run}_3.fastq.gz" ]; then bc=3; cd=4
    else echo "Error: no barcode read (_1 or _3) found for $run"; exit 1; fi
    r1="${FASTQ_SRC_DIR}/${run}_${bc}.fastq.gz"
    r2="${FASTQ_SRC_DIR}/${run}_${cd}.fastq.gz"
    for f in "$r1" "$r2"; do
        if [ ! -f "$f" ]; then echo "Error: FASTQ not found: $f"; exit 1; fi
    done
    ln -sf "$r1" "${FQ_DIR}/${SAMPLE}_S1_${L}_R1_001.fastq.gz"
    ln -sf "$r2" "${FQ_DIR}/${SAMPLE}_S1_${L}_R2_001.fastq.gz"
    echo "  lane $L  <- $run (R1=_${bc} R2=_${cd})"
done
echo "Staged ${lane} lane(s) for $SAMPLE in $FQ_DIR"

# Clear any stale partial run dir so cellranger does not refuse to start.
[ -d "$RUN_DIR" ] && rm -rf "$RUN_DIR"

"$CELLRANGER_BIN" count \
    --id="${genome_name}_${SAMPLE}" \
    --fastqs="$FQ_DIR" \
    --sample="$SAMPLE" \
    --transcriptome="$TRANSCRIPTOME" \
    --chemistry=SC5P-R2 \
    --nosecondary

if [ "${OUTS_ONLY:-FALSE}" = "TRUE" ]; then
    echo "OUTS_ONLY=TRUE: removing non-outs entries under $RUN_DIR"
    find "$RUN_DIR" -mindepth 1 -maxdepth 1 ! -name outs -exec rm -rf {} +
fi

echo "cellranger done: ${RUN_DIR}/outs/"
