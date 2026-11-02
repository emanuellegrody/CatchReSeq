#!/bin/bash
# memoryseq_nfcore_rnaseq.sh  (Figure4 / LD002 MemorySeq clone panel)
# Runs pipelineScripts/memoryseq/nfcore_rnaseq_worker.sh over the 15 LD002
# clone libraries to produce the gene-level TPM matrix the heritable-gene call
# reads. This is the first of the two LD002 steps:
#
#   memoryseq_nfcore_rnaseq.sh   -> star_salmon/salmon.merged.gene_tpm.tsv
#   memoryseq_heritable_genes.sh -> potential_heritable_genes_filtered.csv
#
# Panel J of Figure4.Rmd reads the second file.
#
# Unlike every other wrapper here this one does not sbatch: nf-core/rnaseq
# submits its own jobs from the Nextflow process, which has to stay alive for
# the duration. Run it inside a screen session on a login node, or inside an
# allocation long enough to cover the whole pipeline.
#
# LD002 is human (GRCh38), not macaque. The clones are a human CD4 T line;
# the EGS024 in vitro arm its genes are compared against in Panel J is aligned
# to Mmul_10. Panel J joins the two on gene symbol -- see the note in that
# chunk.

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can
# be cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

WORKER="${PIPELINE_DIR}/memoryseq/nfcore_rnaseq_worker.sh"

# Filesystem roots for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data bash memoryseq_nfcore_rnaseq.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"
GENOME_ROOT="${GENOME_ROOT:-${DATA_ROOT}/genomes}"

# The committed sheet holds bare FASTQ names, which the worker joins to
# FASTQ_ROOT. The reads were sequenced as 20240925_LD02/Emmie/ and are filed
# under rawData/LD002/ here, matching every other accession in the repository.
export SAMPLE_SHEET="${SCRIPT_DIR}/memoryseq_nfcore_rnaseq.csv"
export FASTQ_ROOT="${FASTQ_ROOT:-${DATA_ROOT}/rawData/LD002/Emmie}"

export OUTPUT_DIR="${DATA_ROOT}/extractedData/LD002/nf-core"
export GTF="${GENOME_ROOT}/GRCh38/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
export FASTA="${GENOME_ROOT}/GRCh38/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
export WORK_DIR="${DATA_ROOT}/extractedData/LD002/nfcore_work"

export NFCORE_VERSION="3.17.0"
export SLURM_ACCOUNT="b1042"
export SLURM_PARTITION="genomics"

[ -f "$WORKER" ] || { echo "Error: file not found: $WORKER"; exit 1; }
[ -f "$SAMPLE_SHEET" ] || { echo "Error: file not found: $SAMPLE_SHEET"; exit 1; }

# The reference LD002 was aligned against is the public 10x GRCh38 build, which
# is not one of the references under GENOME_ROOT -- those are the two custom
# cellranger builds the viral arms need. It has to be fetched before this runs:
#
#   curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
#
# GRCh38_HXB2/ is on disk and is GRCh38 with the HXB2 contig appended, so it
# aligns LD002 equivalently apart from reads that would map to HXB2 -- of which
# an uninfected human clone panel has none. Point GTF and FASTA at it to run
# without the download, but note in any methods text which build was used,
# because TPM is computed over whichever annotation is supplied:
#
#   GTF=${GENOME_ROOT}/GRCh38_HXB2/genes/genes.gtf.gz \
#   FASTA=${GENOME_ROOT}/GRCh38_HXB2/fasta/genome.fa \
#     bash memoryseq_nfcore_rnaseq.sh
for reference in "$GTF" "$FASTA"; do
    [ -f "$reference" ] || {
        echo "Error: reference not found: ${reference}"
        echo "       See the reference note in $(basename "${BASH_SOURCE[0]}")."
        exit 1
    }
done

bash "$WORKER"
