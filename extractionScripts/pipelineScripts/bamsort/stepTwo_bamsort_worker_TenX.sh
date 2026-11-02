#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16G
#SBATCH -t 8:00:00
#SBATCH --verbose

# Per-sample worker for single-cell TenX libraries: alignment -> splice
# extract+filter. Submitted by a stepTwo_bamsort_TenX.sh wrapper; not run
# directly. All configuration is provided via env vars exported by the wrapper.
#
# The amplicon worker (stepTwo_bamsort_worker.sh) builds every BAM path as
# <root><sample><target><end> and loops over a non-empty target list, because a
# Ca*Re sample is several amplicon libraries. A TenX sample is one library with
# no target suffix, so this worker takes one BAM path per sample and runs the
# same three stages over it:
#   1. bamsort_alignment_reads.py over every region in COORDS_CSV
#   2. bamsort_splice.py           junctions on the viral contig
#   3. bamsort_splice_sites.py     canonical-donor filter -> the correction file
#        stepThree reads as <sample>_splicesites_filtered.csv, plus every
#        junction as <sample>_splicesites_all.csv
#
# Naming matters downstream: the alignment tables must be written as
# <sample>_bamsort_alignment_<region>.csv with no target token in the sample
# part. extract_category reads that token to decide whether a region is
# isoform-resolved, and its absence is what makes a region fall to the TenX
# reading (A1_D4 -> any-D4, D4_A7 -> US-SS, A7_LTR -> any-A7).
#
# Time is set at 8h because bamsort_splice.py streams the whole BAM rather than
# fetching the viral contig, and the TenX BAMs run 3-23 GB.

set -euo pipefail

REQUIRED_VARS=(
    SAMPLE ALIGNMENT_BAM
    COORDS_CSV OUTPUT_ROOT ALIGNMENT_SCRIPT
    SPLICE_SCRIPT SPLICE_FILTER_SCRIPT SPLICE_OUTPUT_DIR SPLICE_FILTER_DONOR
)
for v in "${REQUIRED_VARS[@]}"; do
    if [ -z "${!v:-}" ]; then
        echo "Error: required env var $v not set (use the stepTwo_bamsort_TenX.sh wrapper)"
        exit 1
    fi
done

# Viral contig to pull alignment reads and junctions from.
VIRAL_CONTIG="${VIRAL_CONTIG:-mac239}"

# Shortest N-CIGAR gap called a splice junction. Must match SPLICE_MIN_INTRON in
# the Ca*Re stepTwo or the two arms are not comparable: the script default of
# 1000 is above the D3 acceptors (D3-A3 528, D3-A4 600, D3-A5 741).
SPLICE_MIN_INTRON="${SPLICE_MIN_INTRON:-500}"

module purge
source activate cellranger
module load python
module load samtools

for f in "$COORDS_CSV" "$ALIGNMENT_SCRIPT" "$SPLICE_SCRIPT" "$SPLICE_FILTER_SCRIPT"; do
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
        echo "Error: required file not found: $f"
        exit 1
    fi
done

if [ ! -f "$ALIGNMENT_BAM" ]; then
    echo "Error: alignment BAM not found: $ALIGNMENT_BAM"
    exit 1
fi

mkdir -p "$OUTPUT_ROOT" "$SPLICE_OUTPUT_DIR"

# =============================================================================
# Step 1: alignment (all coords)
# =============================================================================
echo "========== STEP 1: bamsort_alignment_reads ($SAMPLE) =========="

SUMMARY_FILE="${OUTPUT_ROOT}${SAMPLE}_alignment_summary.txt"
echo "Sample,Target,Start,End,Rows,Total_Reads" > "$SUMMARY_FILE"

python3 <<EOF
import csv, subprocess, os, sys

bam_file = "${ALIGNMENT_BAM}"
sample = "${SAMPLE}"

coords = []
with open("${COORDS_CSV}", "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        keys = list(row.keys())
        coords.append({
            "target": row.get("Target", row.get(keys[0], "")).strip(),
            "start":  row.get("Start",  row.get(keys[1], "")).strip(),
            "end":    row.get("End",    row.get(keys[2], "")).strip(),
        })

print(f"Found {len(coords)} coordinate ranges")

for coord in coords:
    target, start, end = coord["target"], coord["start"], coord["end"]
    output_file = f"${OUTPUT_ROOT}{sample}_bamsort_alignment_{target}.csv"
    print(f"  Region: {target} ({start}-{end})")

    # --windows passes the whole coordinate table so the script can charge a
    # boundary-spanning read to one region only (primary_count).
    cmd = ["python3", "${ALIGNMENT_SCRIPT}", bam_file, "${VIRAL_CONTIG}", start, end, output_file,
           "--windows", "${COORDS_CSV}"]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        if os.path.exists(output_file):
            with open(output_file, "r") as f:
                rows = sum(1 for _ in f) - 1
            total_reads = 0
            with open(output_file, "r") as f:
                reader = csv.reader(f); next(reader)
                for row in reader:
                    if len(row) > 2:
                        try: total_reads += int(row[2])
                        except (ValueError, IndexError): pass
            with open("${SUMMARY_FILE}", "a") as sf:
                sf.write(f"{sample},{target},{start},{end},{rows},{total_reads}\n")
            print(f"    Success: {rows} unique barcodes, {total_reads} reads")
        else:
            print("    ERROR: output not created")
            with open("${SUMMARY_FILE}", "a") as sf:
                sf.write(f"{sample},{target},{start},{end},ERROR,ERROR\n")
            sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"    ERROR: {e.stderr}")
        with open("${SUMMARY_FILE}", "a") as sf:
            sf.write(f"{sample},{target},{start},{end},ERROR,ERROR\n")
        sys.exit(1)
EOF

echo "Alignment complete. Summary: $SUMMARY_FILE"

# =============================================================================
# Step 2a: splice extraction
# =============================================================================
echo "========== STEP 2a: bamsort_splice ($SAMPLE) =========="

SPLICE_CSV="${SPLICE_OUTPUT_DIR}/${SAMPLE}_splicesites.csv"

python3 "$SPLICE_SCRIPT" "$ALIGNMENT_BAM" "$SPLICE_CSV" \
    --reference "$VIRAL_CONTIG" \
    --min-intron "$SPLICE_MIN_INTRON"

# =============================================================================
# Step 2b: splice site filter -> correction file
# =============================================================================
if [ -n "${PROCESS_CORRECTION_DIR:-}" ]; then
    echo "========== STEP 2b: bamsort_splice_sites filter ($SAMPLE) =========="
    mkdir -p "$PROCESS_CORRECTION_DIR"

    # An empty --targets entry makes the filter read <sample>_splicesites.csv and
    # write <sample>_splicesites_filtered.csv, which is the name process_bamsort
    # looks for. The amplicon run passes real suffixes here and concatenates the
    # per-target files afterwards; a TenX sample has only the one file, so it
    # lands on the right name directly.
    python3 "$SPLICE_FILTER_SCRIPT" \
        "$SPLICE_OUTPUT_DIR" \
        "$PROCESS_CORRECTION_DIR" \
        --samples "$SAMPLE" \
        --targets "" \
        --filter-donor "$SPLICE_FILTER_DONOR" \
        --sites "${SPLICE_SITES:-mac239}"

    # _filtered is the correction file; _all keeps every junction.
    for KIND in filtered all; do
        COMBINED_FILE="${PROCESS_CORRECTION_DIR%/}/${SAMPLE}_splicesites_${KIND}.csv"
        if [ -f "$COMBINED_FILE" ]; then
            echo "Splice junctions (${KIND}): $COMBINED_FILE ($(($(wc -l < "$COMBINED_FILE") - 1)) records)"
        else
            echo "Warning: no ${KIND} splice file written for $SAMPLE"
        fi
    done
fi

echo "stepTwo finished for sample: $SAMPLE"
echo "Run stepThree_process_TenX.sh once every sample's stepTwo job has finished."
