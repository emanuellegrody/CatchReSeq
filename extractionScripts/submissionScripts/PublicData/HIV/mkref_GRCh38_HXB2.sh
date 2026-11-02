#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=64G
#SBATCH -t 8:00:00
#SBATCH --job-name=mkref_GRCh38_HXB2
#SBATCH --output=logs/mkref_GRCh38_HXB2_%j.txt
#SBATCH --error=logs/mkref_GRCh38_HXB2_%j.err
#SBATCH --verbose

# Build a joint human + HXB2 cellranger reference (mirrors the Mmul_10_mac239
# combined reference): HXB2 (K03455.1, 9719 bp) appended to GRCh38 as a single
# contig with one whole-genome single-exon gene, so cellranger aligns host +
# HIV in one pass and stamps CB/UB/xf tags on HIV reads.
#
# Run: cd submissionScripts/HIV && sbatch mkref_GRCh38_HXB2.sh

set -euo pipefail

# cellranger binary. Override to point at your own install:
#   CELLRANGER=/path/to/cellranger bash mkref_GRCh38_HXB2.sh
CELLRANGER="${CELLRANGER:-/home/egy2296/packages/cellranger-7.2.0/cellranger}"
# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   GENOME_ROOT=/my/genomes bash mkref_GRCh38_HXB2.sh
GENOME_ROOT="${GENOME_ROOT:-/projects/b1042/GoyalLab/egrody/genomes}"

GENOMES="${GENOME_ROOT}"
GRCH38="${GENOMES}/refdata-gex-GRCh38-2020-A"
BUILD="${GENOMES}/GRCh38_HXB2_inputs"
COMBINED_FA="${BUILD}/combined_genome.fa"
COMBINED_GTF="${BUILD}/combined_genes.gtf"

source activate cellranger

# --- combined genome fasta: GRCh38 contigs + HXB2 appended --------------------
echo "$(date) building combined fasta"
cat "${GRCH38}/fasta/genome.fa" > "$COMBINED_FA"
[ -n "$(tail -c1 "$COMBINED_FA")" ] && echo >> "$COMBINED_FA"   # ensure trailing newline
cat "${BUILD}/HXB2.fa" >> "$COMBINED_FA"

# --- combined GTF: GRCh38 genes + HXB2 single-exon gene -----------------------
echo "$(date) building combined gtf"
cat "${GRCH38}/genes/genes.gtf" > "$COMBINED_GTF"
[ -n "$(tail -c1 "$COMBINED_GTF")" ] && echo >> "$COMBINED_GTF"
cat "${BUILD}/HXB2.gtf" >> "$COMBINED_GTF"

# sanity: HXB2 must appear once in each combined input
echo "HXB2 fasta headers: $(grep -c '^>HXB2$' "$COMBINED_FA")"
echo "HXB2 gtf lines:     $(grep -c $'^HXB2\t' "$COMBINED_GTF")"

# --- cellranger mkref (creates ${GENOMES}/GRCh38_HXB2) ------------------------
echo "$(date) running cellranger mkref"
cd "$GENOMES"
"$CELLRANGER" mkref \
    --genome=GRCh38_HXB2 \
    --fasta="$COMBINED_FA" \
    --genes="$COMBINED_GTF" \
    --nthreads=12 \
    --memgb=60

echo "$(date) mkref done -> ${GENOMES}/GRCh38_HXB2"
