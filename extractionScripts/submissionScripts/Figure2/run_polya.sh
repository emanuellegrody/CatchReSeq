#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=run_polya
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH -t 4:00:00
#SBATCH --output=logs/run_polya.%j.out
#SBATCH --error=logs/run_polya.%j.err
#SBATCH --verbose

set -euo pipefail

mkdir -p logs

# -----------------------------------------------------------------------------
# Which steps to run: pass a digit string (e.g. 1, 12, 13, 123). Default = all.
#   1 = polya_comparative_stats.py   (SuppFig 2A)
#   2 = bam_polya_virus.py           (Figure 2D)
#   3 = bam_polya_host.py            (SuppFig 2C)
#   4 = polya_splice_proximity.py    (HXB2/SIVmac239 polyA vs splice sites)
#   5 = bam_polya_longread.py        (SuppFig 2D)
# Submit:  sbatch run_polya.sh 13     Run locally:  ./run_polya.sh 2
# -----------------------------------------------------------------------------
STEPS="${1:-12345}"
if ! [[ "$STEPS" =~ ^[12345]+$ ]]; then
    echo "ERROR: STEPS must contain only digits 1-5 (got '$STEPS')" >&2
    exit 1
fi
run_step() { [[ "$STEPS" == *"$1"* ]]; }
echo "Running steps: $STEPS"

module purge
eval "$(conda shell.bash hook)"
source activate cellranger
module load python
module load samtools

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
# Script locations are resolved from this script's own path so the repo can be
# cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

# Filesystem roots for every input read and output written below. Override to
# run against a different copy of the data:
#   DATA_ROOT=/my/data bash run_polya.sh 2
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

HOST_GENOME="${DATA_ROOT}/genomes/Mmul_10_mac239/fasta/genome.fa"
HOST_GTF="${DATA_ROOT}/genomes/Mmul_10_mac239/genes/protein_coding_genes.gtf"
GENE_LIST="${DATA_ROOT}/genomes/Mmul_10_mac239/genes/pc_gene_names.txt"
VIRAL_GENOME="${DATA_ROOT}/genomes/mac239/fasta/genome.fa"

HOST_BAM="${DATA_ROOT}/extractedData/Figure2/300cy/counts/polyA/Mmul_10_mac239/JK85_D13_STAR_Aligned.sortedByCoord.out.bam"
SIV_BAM="${DATA_ROOT}/extractedData/Figure2/300cy/counts/polyA/STAR/JK85_D13_STAR_Aligned.sortedByCoord.out.bam"
# Cell Ranger possorted BAM carries xf/CB/UB tags. The polyA scripts use this
# to drop reads whose xf is not in {17, 25}
CELLRANGER_BAM="${DATA_ROOT}/extractedData/Figure2/300cy/counts/Mmul_10_mac239annot_JK85_D13/outs/possorted_genome_bam.bam"
# Seurat-filtered cellID whitelist for the host polyA sweep (step 3). Must
# include a 'cellID' column with the Cell Ranger CB format (XXX-1).
SEURAT_CELLIDS="${DATA_ROOT}/extractedData/Figure1/singleCell/Seurat/UMAPcoords/D13UMAP_coords.csv"
# Long-read (ONT cDNA) inputs for step 5. These barcodes belong to the Invitro
# Seurat object, not D13: 234 of 339 long-read cellIDs are in InvitroUMAP_coords,
# and none are in D13UMAP_coords. The BAM is aligned to the 10279bp mac239, which
# bam_polya_longread.py rebuilds from the 9602bp FASTA below.
LONGREAD_BAM="${DATA_ROOT}/extractedData/Figure2/longRead/nf-core/cDNA/bam/dedup/cDNA.dedup.sorted.bam"
LONGREAD_FASTA="${DATA_ROOT}/genomes/Mmul_10_mac239/inputs/mac239.fa"
LONGREAD_SEURAT_CELLIDS="${DATA_ROOT}/extractedData/Figure1/singleCell/Seurat/UMAPcoords/InvitroUMAP_coords.csv"

# -----------------------------------------------------------------------------
# Outputs
# -----------------------------------------------------------------------------
# Results are data, so they go under DATA_ROOT rather than into the repo. This is
# the directory Figure2.Rmd and SuppFigs.Rmd read the polyA tables from.
OUT_DIR="${POLYA_OUT_DIR:-${DATA_ROOT}/extractedData/Figure2/polyA}"
COMP_DIR="${OUT_DIR}/polya_comparative"
BAM_OUT_DIR="${OUT_DIR}/bam_polya_virus"
HOST_OUT_DIR="${OUT_DIR}/bam_polya_host"
SPLICE_OUT_DIR="${OUT_DIR}/polya_splice_proximity"
LONGREAD_OUT_DIR="${OUT_DIR}/bam_polya_longread"

mkdir -p "$COMP_DIR" "$BAM_OUT_DIR" "$HOST_OUT_DIR" "$SPLICE_OUT_DIR" \
         "$LONGREAD_OUT_DIR"

# -----------------------------------------------------------------------------
# 1) Genome-level polyA density: host vs viral (input for SuppFig 2A)
# -----------------------------------------------------------------------------
if run_step 1; then
    echo "========================================"
    echo "[1] polya_comparative_stats.py"
    echo "========================================"
    python3 "${PIPELINE_DIR}/polya/polya_comparative_stats.py" \
        -host  "$HOST_GENOME" \
        -viral "$VIRAL_GENOME" \
        -gtf   "$HOST_GTF" \
        -genes "$GENE_LIST" \
        -out   "$COMP_DIR"
fi

# -----------------------------------------------------------------------------
# 2) SIV BAM -> per-read polyA structure CSV (input for Figure 2D)
# -----------------------------------------------------------------------------
if run_step 2; then
    echo ""
    echo "========================================"
    echo "[2] bam_polya_virus.py"
    echo "========================================"
    python3 "${PIPELINE_DIR}/polya/bam_polya_virus.py" \
        "$SIV_BAM" \
        "${BAM_OUT_DIR}/d400_class" \
        --cellranger-bam "$CELLRANGER_BAM" \
        --xf-allow 0,17,25
fi

# -----------------------------------------------------------------------------
# 3) Host BAM -> per-gene explained reads sweep (input for SuppFig 2C)
# -----------------------------------------------------------------------------
if run_step 3; then
    echo ""
    echo "========================================"
    echo "[3] bam_polya_host.py"
    echo "========================================"
    python3 "${PIPELINE_DIR}/polya/bam_polya_host.py" \
        "$HOST_BAM" \
        "$GENE_LIST" \
        "$HOST_GENOME" \
        "$HOST_GTF" \
        --out "${HOST_OUT_DIR}/host_polya_sweep.csv" \
        --cellranger-bam "$CELLRANGER_BAM" \
        --seurat-cellids "$SEURAT_CELLIDS" \
        --min-umi-reads 2
fi

# -----------------------------------------------------------------------------
# 4) Internal polyA tracts upstream of donors vs downstream of acceptors,
#    for HXB2 and SIVmac239. Genome FASTAs and splice-site tables are defaults
#    inside the script; nothing but an output directory is required.
# -----------------------------------------------------------------------------
if run_step 4; then
    echo ""
    echo "========================================"
    echo "[4] polya_splice_proximity.py"
    echo "========================================"
    python3 "${PIPELINE_DIR}/polya/polya_splice_proximity.py" \
        -out "$SPLICE_OUT_DIR"
fi

# -----------------------------------------------------------------------------
# 5) Long-read BAM -> per-molecule priming-site classification plus the two
#    control end sets (input for SuppFig 2D)
# -----------------------------------------------------------------------------
if run_step 5; then
    echo ""
    echo "========================================"
    echo "[5] bam_polya_longread.py"
    echo "========================================"
    python3 "${PIPELINE_DIR}/polya/bam_polya_longread.py" \
        "$LONGREAD_BAM" \
        "$LONGREAD_FASTA" \
        --out-dir "$LONGREAD_OUT_DIR" \
        --seurat-cellids "$LONGREAD_SEURAT_CELLIDS"
fi

echo ""
echo "All done."
