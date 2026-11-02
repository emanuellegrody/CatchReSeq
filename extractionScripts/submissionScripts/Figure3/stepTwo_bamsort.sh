#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH -t 2:00:00
#SBATCH --job-name=egs026-bamsort
#SBATCH --output=logs/stepTwo_bamsort_%j.txt
#SBATCH --error=logs/stepTwo_bamsort_%j.err
#SBATCH --verbose

# HIV bamsort step (mirrors the Mmul_mac239 TenX bamsort). For each EGS026
# pool, from the joint-genome (GRCh38_HXB2) cellranger BAM (CB/UB/xf tags on
# HIV reads):
#   1. alignment: HXB2 reads/UMIs per HIV splice region ->
#        <sample>_bamsort_alignment_<target>.csv   (UMI,cellID,count)
#   2. splice sites: junctions (>=20M + >=1000N + M, xf:i:25) per molecule ->
#        splice/<sample>_splicesites.csv   (CB,UB,...,Donor,Acceptor,...)
# The alignment tables feed stepThree_process.sh (process_bamsort) and then
# set_minimums in the analysis scripts; the splice sites
# confirm MS (D4->A7 donor) like the SIV SPLICE_FILTER_DONOR step.
#
# Run AFTER stepOne_cellranger.sh finishes. Idempotent and re-runnable: samples
# whose BAM is not yet present are skipped, so it can be re-run to pick up
# newly-finished samples.  Run: cd submissionScripts/Figure3 && sbatch stepTwo_bamsort.sh

set -uo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
MANIFEST="${SCRIPT_DIR}/samples_manifest.csv"
COORDS_CSV="${SCRIPT_DIR}/../../pipelineScripts/bamsort/coordinates/bamsort_coordinates_HXB2.csv"
ALIGN_SCRIPT="${SCRIPT_DIR}/../../pipelineScripts/bamsort/bamsort_alignment_reads.py"
SPLICE_SCRIPT="${SCRIPT_DIR}/../../pipelineScripts/bamsort/bamsort_splice.py"

### UPDATE PATHS
# Filesystem root for every input read and output written below. Override to
# run this against a different copy of the data:
#   DATA_ROOT=/my/data bash stepTwo_bamsort.sh
DATA_ROOT="${DATA_ROOT:-/projects/b1042/GoyalLab/egrody}"

COUNTS_ROOT="${DATA_ROOT}/extractedData/Figure3/singleCell/counts"
GENOME_PREFIX="GRCh38_HXB2"          # cellranger run-dir prefix (basename of the reference)
VIRAL_CONTIG="HXB2"                  # contig name in the joint reference / BAM
OUTPUT_ROOT="${DATA_ROOT}/extractedData/Figure3/singleCell/bamsort"

for f in "$MANIFEST" "$COORDS_CSV" "$ALIGN_SCRIPT" "$SPLICE_SCRIPT"; do
    if [ ! -f "$f" ]; then
        # bamsort_alignment_reads{,_host}.py are withheld from the public
        # release, so say that rather than failing as a generic missing file.
        case "$(basename "$f")" in
            bamsort_alignment_reads.py|bamsort_alignment_reads_host.py)
                echo "Error: $(basename "$f") is not present."
                echo
                echo "It is withheld from the public release of this repository and is"
                echo "available on request. See PRIVATE_SCRIPTS.md for what it does, its"
                echo "interface, and how to ask for it. Drop the file in"
                echo "extractionScripts/pipelineScripts/ and re-run; nothing else changes."
                exit 2
                ;;
        esac
        echo "Error: missing $f"
        exit 1
    fi
done

source activate cellranger
module load python
module load samtools

mkdir -p "$OUTPUT_ROOT"
SPLICE_DIR="${OUTPUT_ROOT}/splice"; mkdir -p "$SPLICE_DIR"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
SUMMARY="${OUTPUT_ROOT}/all_samples_alignment_summary.csv"
echo "sample,target,start,end,barcodes,total_reads" > "$SUMMARY"

# splice regions (skip header)
mapfile -t COORD_LINES < <(tail -n +2 "$COORDS_CSV")

# optional sample filter: pass sample names as args; else all in manifest
WANTED=" $* "

processed=0; skipped=0
while IFS=, read -r SAMPLE _rest; do
    SAMPLE="$(echo "$SAMPLE" | tr -d '[:space:]')"
    [ -z "$SAMPLE" ] && continue
    if [ "$#" -gt 0 ] && [[ "$WANTED" != *" $SAMPLE "* ]]; then continue; fi

    BAM="${COUNTS_ROOT}/${GENOME_PREFIX}_${SAMPLE}/outs/possorted_genome_bam.bam"
    if [ ! -f "$BAM" ]; then
        echo "SKIP $SAMPLE: cellranger BAM not found yet ($BAM)"
        skipped=$((skipped+1)); continue
    fi

    echo "========== $SAMPLE =========="
    for line in "${COORD_LINES[@]}"; do
        IFS=, read -r TARGET START END <<< "$line"
        [ -z "${TARGET:-}" ] && continue
        OUT="${OUTPUT_ROOT}/${SAMPLE}_bamsort_alignment_${TARGET}.csv"
        # --windows charges a read spanning two adjacent targets to one of them
        if python3 "$ALIGN_SCRIPT" "$BAM" "$VIRAL_CONTIG" "$START" "$END" "$OUT" \
               --windows "$COORDS_CSV"; then
            rows=$(($(wc -l < "$OUT") - 1))
            reads=$(awk -F, 'NR>1{s+=$3} END{print s+0}' "$OUT")
            echo "$SAMPLE,$TARGET,$START,$END,$rows,$reads" >> "$SUMMARY"
            echo "  $TARGET ($START-$END): $rows barcodes, $reads reads"
        else
            echo "  ERROR on $TARGET ($START-$END)"
            echo "$SAMPLE,$TARGET,$START,$END,ERROR,ERROR" >> "$SUMMARY"
        fi
    done

    # --- splice sites: subset to the viral contig, then bamsort_splice --------
    SPLICE_BAM="${TMP}/${SAMPLE}.${VIRAL_CONTIG}.bam"
    samtools view -b "$BAM" "$VIRAL_CONTIG" > "$SPLICE_BAM" 2>/dev/null
    samtools index "$SPLICE_BAM"
    SPLICE_OUT="${SPLICE_DIR}/${SAMPLE}_splicesites.csv"
    if python3 "$SPLICE_SCRIPT" "$SPLICE_BAM" "$SPLICE_OUT" --reference "$VIRAL_CONTIG" >/dev/null 2>&1; then
        echo "  splice junctions: $(($(wc -l < "$SPLICE_OUT") - 1))"
    else
        echo "  ERROR extracting splice sites"
    fi
    rm -f "$SPLICE_BAM" "${SPLICE_BAM}.bai"

    processed=$((processed+1))
done < <(tail -n +2 "$MANIFEST")

echo
echo "Done. processed=$processed skipped=$skipped"
echo "Alignment summary: $SUMMARY"
echo "Splice sites:      $SPLICE_DIR"
