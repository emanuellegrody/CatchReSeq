#!/bin/bash
# memoryseq_variability.sh  (Figure4 / LD002 MemorySeq clone panel)
# Submits pipelineScripts/memoryseq/memoryseq_variability.R, which turns the
# nf-core gene TPM matrix into the over-dispersed-gene list behind Figure4.Rmd
# Panel J.
#
# Depends only on memoryseq_nfcore_rnaseq.sh having produced the matrices. It
# reads two tables and fits two GLMs, so it is neither memory- nor CPU-bound.
#
# The settings below reproduce the published list exactly -- 472 rows, 245 from
# the ExpSample arm and 227 from the Controls arm, the same gene_id/arm pairs as
# potential_heritable_genes.csv. Three of them are defects, documented at their
# definitions in the R script and summarised here. To run the call as the method
# intends:
#
#   CONTROL_SAMPLES= COUNTS_MATRIX= HITS_FROM=ExpSample bash memoryseq_variability.sh
#
# That changes Panel J's numbers; the Panel J prose records by how much.

set -euo pipefail

# Script locations are resolved from this script's own path so the repo can
# be cloned anywhere; nothing below reaches outside it.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${SCRIPT_DIR}/../../pipelineScripts"

R_SCRIPT="${PIPELINE_DIR}/memoryseq/memoryseq_variability.R"

# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data bash memoryseq_variability.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"
SALMON_DIR="${DATA_ROOT}/extractedData/LD002/nf-core/star_salmon"

export TPM_MATRIX="${SALMON_DIR}/salmon.merged.gene_tpm.tsv"

# DEFECT: the published run applied the 1-unit expression floor to raw counts,
# not RPM. Unset COUNTS_MATRIX to apply it to TPM instead.
export COUNTS_MATRIX="${COUNTS_MATRIX-${SALMON_DIR}/salmon.merged.gene_counts.tsv}"

export OUTPUT_DIR="${OUTPUT_DIR:-${DATA_ROOT}/extractedData/LD002/memoryseq/}"
export SAMPLE_PATTERN="^Clone_"

# DEFECT: these four are clones, relabelled as controls by the original script
# (`controls[c(2,9,11,13),]$controls <- 'Controls'`, indexing the clone list in
# the order the samtools_stats files sort). LD002 has no control libraries.
# Unset CONTROL_SAMPLES to fit the 15-clone panel whole.
export CONTROL_SAMPLES="${CONTROL_SAMPLES-Clone_02|Clone_17|Clone_23|Clone_25}"

# DEFECT: the published list keeps hits from both arms, so some of the rows
# Panel J reads were called from the 4-library fit. Set HITS_FROM=ExpSample to
# keep one arm, as the paper's own script does.
export HITS_FROM="${HITS_FROM:-both}"

# Hits at rearranged TCR/Ig loci are clone-specific by genomic rearrangement,
# not by gene regulation, so they are excluded. This pattern is the one that
# reproduces the published file; it is slightly over-broad, also catching TROAP
# and TRIP13. A precise version:
#   EXCLUDE_PATTERN='^(IG[HKL]|TR[ABGD]|IGLL)' bash memoryseq_variability.sh
export EXCLUDE_PATTERN="${EXCLUDE_PATTERN-^(IG|TR)}"

export MEAN_FLOOR="1"
export MEAN_TPM_CUTOFF="1.5"
export PERCENTILE_CUTOFF="0.98"

[ -f "$R_SCRIPT" ] || { echo "Error: file not found: $R_SCRIPT"; exit 1; }
[ -f "$TPM_MATRIX" ] || {
    echo "Error: TPM matrix not found: $TPM_MATRIX"
    echo "       Run memoryseq_nfcore_rnaseq.sh first."
    exit 1
}

mkdir -p logs

sbatch \
    --job-name="ld002herit" \
    --account=b1042 \
    --partition=genomics \
    --time=00:20:00 \
    --mem=16G \
    --nodes=1 \
    --ntasks=1 \
    --output="logs/memoryseq_variability_%j.txt" \
    --error="logs/memoryseq_variability_%j.err" \
    --wrap="module purge && module load R/4.4.0 && Rscript ${R_SCRIPT}"

echo "memoryseq variability job submitted."
echo "Output: ${OUTPUT_DIR}potential_heritable_genes_filtered.csv (Panel J reads this)"
