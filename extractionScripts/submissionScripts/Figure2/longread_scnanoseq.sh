#!/bin/bash
# longread_scnanoseq.sh  (Figure2 / EGS007 long-read arm)
# Runs pipelineScripts/longread/scnanoseq_worker.sh over the five EGS007
# Nanopore libraries to produce the barcode-tagged, deduplicated BAMs and the
# IsoQuant counts the Figure2 long-read panels read:
#
#   longread_scnanoseq.sh      -> <sample>/bam/dedup/<sample>.dedup.sorted.bam
#   bamsort_isoform_longread.sh -> longRead/bamsort/isoform/cDNA_isoforms.csv
#   run_polya.sh                -> polyA/bam_polya_longread/
#
# Figure2.Rmd reads the dedup BAM directly for the coverage panel, and the two
# tables above for the isoform and polyA panels.
#
# Unlike every other wrapper here this one does not sbatch: nf-core/scnanoseq
# submits its own jobs from the Nextflow process, which has to stay alive for
# the duration. Run it inside a screen session on a login node, or inside an
# allocation long enough to cover the whole pipeline.
#
# EGS007 is the in vitro SIV sample. All five libraries are 10x 3' v3 chemistry
# and are quantified against the viral contig alone, not against Mmul_10.

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can
# be cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

WORKER="${PIPELINE_DIR}/longread/scnanoseq_worker.sh"

# Filesystem roots for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data bash longread_scnanoseq.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"
GENOME_ROOT="${GENOME_ROOT:-${DATA_ROOT}/genomes}"

# The committed sheet holds bare FASTQ names, which the worker joins to
# FASTQ_ROOT.
export SAMPLE_SHEET="${SCRIPT_DIR}/longread_scnanoseq.csv"
export FASTQ_ROOT="${FASTQ_ROOT:-${DATA_ROOT}/rawData/EGS007/Sequencing/longRead}"

export OUTPUT_DIR="${DATA_ROOT}/extractedData/Figure2/longRead/nf-core"
export WORK_DIR="${DATA_ROOT}/extractedData/Figure2/longRead/nfcore_work"

# The viral-only reference. The 9602 bp cellranger builds under GENOME_ROOT
# drop the 5' LTR U3 and the 3' U5 so that short reads cannot multimap between
# the two LTRs, a problem long reads do not have, so this arm aligns against
# the full-length 10279 bp contig instead. The two are therefore not
# interchangeable: every coordinate in bamsort_isoform_longread.sh is offset
# against the 10279 bp contig, and the dedup BAM header records that length.
# See the reference note below.
export FASTA="${FASTA:-${GENOME_ROOT}/mac239_longread/mac239.fa}"
export GTF="${GTF:-${GENOME_ROOT}/mac239_longread/mac239.gtf}"

# Pipeline revision and parameters, as published.
export NFCORE_VERSION="1.0.0"
export BARCODE_FORMAT="10X_3v3"
export KMER_SIZE="14"
export MIN_Q_SCORE="10"
export MIN_LENGTH="1"
export RETAIN_INTRONS="true"

# SLURM account and partition. These are the Northwestern Quest values; override
# them for your own site:
#   SLURM_ACCOUNT=myacct SLURM_PARTITION=mypart bash <this script>
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-b1042}"
export SLURM_PARTITION="${SLURM_PARTITION:-genomics}"

[ -f "$WORKER" ] || { echo "Error: file not found: $WORKER"; exit 1; }
[ -f "$SAMPLE_SHEET" ] || { echo "Error: file not found: $SAMPLE_SHEET"; exit 1; }

# The 10279 bp mac239 FASTA this arm was aligned against is not on disk. What is
# under GENOME_ROOT/mac239/ is the 9602 bp cellranger build, which the viral
# arms of Figures 1 and 3 use and which will not reproduce these coordinates.
# Rebuild the full-length contig from its accession, place it as the two files
# named above, and the pipeline runs unchanged. The published dedup BAM header
# records the length to check a rebuild against:
#
#   samtools view -H <sample>.dedup.sorted.bam | grep '^@SQ'
#   @SQ  SN:mac239  LN:10279
for reference in "$GTF" "$FASTA"; do
    [ -f "$reference" ] || {
        echo "Error: reference not found: ${reference}"
        echo "       See the reference note in $(basename "${BASH_SOURCE[0]}")."
        exit 1
    }
done

bash "$WORKER"
