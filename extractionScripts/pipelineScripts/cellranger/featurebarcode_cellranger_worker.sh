#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH --verbose

# Per-POOL cellranger worker for the 10X 3' v3 pools with TotalSeq-B feature
# barcoding: Giron et al. 2024 (GSE254483) and the EGS026 re-sequenced pools.
# Each pool has a Gene Expression library and an Antibody Capture library
# (12 CITE-seq antibodies + 4 donor hashtags). cellranger count runs in
# feature-barcode mode against the joint human GRCh38 + HXB2 reference, so the
# matrix carries genes + HXB2 + ADT + HTO and the possorted BAM has CB/UB/xf
# tags on HIV reads, feeding the bamsort process (like HIVseq/Collora).
#
# The two pool sets differ only in where the Gene Expression reads come from.
# Giron's are an SRA run, named _1/_2 by ENA and staged here as 10X-named
# symlinks. EGS026's are bcl-convert output already named to the 10X convention
# (${TenX_FASTQ_SAMPLE}_S#_L001_R#_001.fastq.gz) and are read in place. EGS026
# has no Antibody Capture library of its own, so it borrows the Giron run for
# the same pool; both libraries come from one GEM well, so their cell barcodes
# agree.
#
# ENA names each run's reads _1 (R1: 16bp cell barcode + 12bp UMI) and _2
# (R2: cDNA for the TenX library, feature barcode for the Antibody Capture library).
#
# Required env vars (from wrapper):
#   SAMPLE            pool name (HIV-A..D); cellranger run-dir suffix
#   FB_RUN            SRA run accession of the Antibody Capture library
#   FB_FASTQ_SRC_DIR  dir containing ${FB_RUN}_1.fastq.gz and ${FB_RUN}_2.fastq.gz
#                     (falls back to FASTQ_SRC_DIR)
#   TRANSCRIPTOME     cellranger reference (joint human+HIV GRCh38_HXB2)
#   FEATURE_REF       feature_reference.csv (ADT + hashtag barcodes)
#   OUTPUT_PATH       dir to cd into / write the run
#
# Gene Expression source -- set exactly one of:
#   TenX_RUN                             SRA run accession; staged as symlinks,
#                                        with FASTQ_SRC_DIR as its source dir
#   TenX_FASTQ_DIR + TenX_FASTQ_SAMPLE   already 10X-named; read in place
#
# Optional:
#   CELLRANGER_BIN    default /home/egy2296/packages/cellranger-7.2.0/cellranger
#   OUTS_ONLY         'TRUE' to keep only outs/ after the run (default FALSE)

set -euo pipefail

# The Antibody Capture source dir is FB_FASTQ_SRC_DIR; Giron's wrapper sets the
# one dir it uses for both libraries as FASTQ_SRC_DIR, so accept that too.
FB_FASTQ_SRC_DIR="${FB_FASTQ_SRC_DIR:-${FASTQ_SRC_DIR:-}}"

for v in SAMPLE FB_RUN FB_FASTQ_SRC_DIR TRANSCRIPTOME FEATURE_REF OUTPUT_PATH; do
    if [ -z "${!v:-}" ]; then echo "Error: required env var $v not set"; exit 1; fi
done

# Exactly one Gene Expression source. Erroring on both-set rather than silently
# preferring one keeps a half-edited wrapper from running against the wrong reads.
if [ -n "${TenX_RUN:-}" ] && [ -n "${TenX_FASTQ_DIR:-}" ]; then
    echo "Error: set either TenX_RUN or TenX_FASTQ_DIR, not both"; exit 1
fi
if [ -z "${TenX_RUN:-}" ] && [ -z "${TenX_FASTQ_DIR:-}" ]; then
    echo "Error: no Gene Expression source: set TenX_RUN, or TenX_FASTQ_DIR + TenX_FASTQ_SAMPLE"
    exit 1
fi
if [ -n "${TenX_RUN:-}" ] && [ -z "${FASTQ_SRC_DIR:-}" ]; then
    echo "Error: TenX_RUN needs FASTQ_SRC_DIR to resolve ${TenX_RUN}_1/_2.fastq.gz"; exit 1
fi
if [ -n "${TenX_FASTQ_DIR:-}" ] && [ -z "${TenX_FASTQ_SAMPLE:-}" ]; then
    echo "Error: TenX_FASTQ_DIR needs TenX_FASTQ_SAMPLE (the FASTQ sample prefix)"; exit 1
fi

CELLRANGER_BIN="${CELLRANGER_BIN:-/home/egy2296/packages/cellranger-7.2.0/cellranger}"
genome_name=$(basename "$TRANSCRIPTOME")
RUN_DIR="${OUTPUT_PATH}/${genome_name}_${SAMPLE}"

# Idempotency: skip if this pool's cellranger run already finished.
if [ -f "${RUN_DIR}/outs/possorted_genome_bam.bam" ]; then
    echo "cellranger already done for $SAMPLE, skipping: ${RUN_DIR}/outs"
    exit 0
fi

[ -f "$FEATURE_REF" ] || { echo "Error: feature reference not found: $FEATURE_REF"; exit 1; }

source activate cellranger
mkdir -p "$OUTPUT_PATH"
cd "$OUTPUT_PATH"

# Stage 10X-named symlinks, one lane per library (_1 -> R1, _2 -> R2).
stage_lib() {           # $1=run  $2=fastq-subdir  $3=sample-prefix  $4=src-dir
    local run="$1" dir="$2" pref="$3" src="$4"
    local r1="${src}/${run}_1.fastq.gz" r2="${src}/${run}_2.fastq.gz"
    [ -f "$r1" ] && [ -f "$r2" ] || { echo "Error: FASTQ not found for $run ($r1 / $r2)"; exit 1; }
    rm -rf "$dir"; mkdir -p "$dir"
    ln -sf "$r1" "${dir}/${pref}_S1_L001_R1_001.fastq.gz"
    ln -sf "$r2" "${dir}/${pref}_S1_L001_R2_001.fastq.gz"
}

# Gene Expression library: staged from an SRA run, or read in place.
if [ -n "${TenX_RUN:-}" ]; then
    TenX_DIR="${OUTPUT_PATH}/fastqs_${SAMPLE}_GEX"
    TenX_SAMPLE="${SAMPLE}tenx"
    stage_lib "$TenX_RUN" "$TenX_DIR" "$TenX_SAMPLE" "$FASTQ_SRC_DIR"
    echo "TenX library: $TenX_RUN staged in $TenX_DIR"
else
    TenX_DIR="$TenX_FASTQ_DIR"
    TenX_SAMPLE="$TenX_FASTQ_SAMPLE"
    [ -d "$TenX_DIR" ] || { echo "Error: TenX FASTQ dir not found: $TenX_DIR"; exit 1; }
    # The TenX FASTQs are read in place; confirm the prefix actually matches something.
    tenx_found=$(find "$TenX_DIR" -maxdepth 1 -name "${TenX_SAMPLE}_S*_R1_001.fastq.gz" | wc -l)
    [ "$tenx_found" -gt 0 ] || { echo "Error: no TenX FASTQs matching ${TenX_SAMPLE} in $TenX_DIR"; exit 1; }
    echo "TenX library: ${TenX_SAMPLE} (${tenx_found} lane(s)) in $TenX_DIR"
fi

# Antibody Capture library: always staged from an SRA run.
FB_DIR="${OUTPUT_PATH}/fastqs_${SAMPLE}_FB"
stage_lib "$FB_RUN" "$FB_DIR" "${SAMPLE}fb" "$FB_FASTQ_SRC_DIR"
echo "FB library:  $FB_RUN staged in $FB_DIR"

# libraries.csv declaring the two 10X libraries for this pool.
LIB_CSV="${OUTPUT_PATH}/libraries_${SAMPLE}.csv"
{
    echo "fastqs,sample,library_type"
    echo "${TenX_DIR},${TenX_SAMPLE},Gene Expression"
    echo "${FB_DIR},${SAMPLE}fb,Antibody Capture"
} > "$LIB_CSV"

# Clear any stale partial run dir so cellranger does not refuse to start.
[ -d "$RUN_DIR" ] && rm -rf "$RUN_DIR"

"$CELLRANGER_BIN" count \
    --id="${genome_name}_${SAMPLE}" \
    --libraries="$LIB_CSV" \
    --feature-ref="$FEATURE_REF" \
    --transcriptome="$TRANSCRIPTOME" \
    --chemistry=SC3Pv3 \
    --nosecondary

if [ "${OUTS_ONLY:-FALSE}" = "TRUE" ]; then
    echo "OUTS_ONLY=TRUE: removing non-outs entries under $RUN_DIR"
    find "$RUN_DIR" -mindepth 1 -maxdepth 1 ! -name outs -exec rm -rf {} +
fi

echo "cellranger done: ${RUN_DIR}/outs/"
