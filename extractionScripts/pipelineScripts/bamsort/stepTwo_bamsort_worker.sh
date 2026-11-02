#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH --verbose

# Per-sample worker: alignment -> splice extract+filter.
# Submitted by Figure1/stepTwo_bamsort.sh; not run directly.
# All configuration is provided via env vars exported by the wrapper.
#
# Combining the per-region alignment tables into per-UMI categories is
# stepThree_process.sh, which gathers across every sample and so has to wait
# for all of these jobs to finish.

set -euo pipefail

REQUIRED_VARS=(
    SAMPLE ALIGNMENT_BAM_ID
    COORDS_CSV INPUT_ROOT INPUT_END OUTPUT_ROOT ALIGNMENT_SCRIPT
    SPLICE_SCRIPT SPLICE_FILTER_SCRIPT SPLICE_INPUT_ROOT SPLICE_INPUT_END
    SPLICE_OUTPUT_DIR SPLICE_TARGETS SPLICE_FILTER_DONOR
)
for v in "${REQUIRED_VARS[@]}"; do
    if [ -z "${!v:-}" ]; then
        echo "Error: required env var $v not set (use stepTwo_bamsort.sh wrapper)"
        exit 1
    fi
done

# Reconstruct splice targets array from |-separated string
IFS='|' read -r -a SPLICE_TARGETS_ARR <<< "$SPLICE_TARGETS"

# Viral contig to pull alignment reads from, i.e. the viral chromosome name in
# the reference the BAM was aligned against.
VIRAL_CONTIG="${VIRAL_CONTIG:-mac239}"

# Shortest N-CIGAR gap called a splice junction. bamsort_splice.py defaults to
# 1000, which is above the D3 acceptors (D3-A3 528, D3-A4 600, D3-A5 741) and so
# hides them entirely.
SPLICE_MIN_INTRON="${SPLICE_MIN_INTRON:-500}"

# One sample may have several amplicon libraries, each with its own BAM. This
# step runs over every ID in ALIGNMENT_BAM_IDS; stepThree's process_bamsort then
# gathers all of their alignment CSVs back together under the one sample.
ALIGNMENT_BAM_IDS="${ALIGNMENT_BAM_IDS:-$ALIGNMENT_BAM_ID}"
IFS='|' read -r -a ALIGNMENT_BAM_IDS_ARR <<< "$ALIGNMENT_BAM_IDS"

module purge
source activate cellranger
module load python
module load samtools

for f in "$COORDS_CSV" "$ALIGNMENT_SCRIPT" \
         "$SPLICE_SCRIPT" "$SPLICE_FILTER_SCRIPT"; do
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

mkdir -p "$OUTPUT_ROOT" "$SPLICE_OUTPUT_DIR"

# =============================================================================
# Step 1: alignment (each amplicon library of this sample, all coords)
# =============================================================================
for ALIGNMENT_BAM_ID in "${ALIGNMENT_BAM_IDS_ARR[@]}"; do
echo "========== STEP 1: bamsort_alignment_reads ($SAMPLE / $ALIGNMENT_BAM_ID) =========="

BAM_FILE="${INPUT_ROOT}${ALIGNMENT_BAM_ID}${INPUT_END}"
if [ ! -f "$BAM_FILE" ]; then
    echo "Error: alignment BAM not found: $BAM_FILE"
    exit 1
fi

SUMMARY_FILE="${OUTPUT_ROOT}${ALIGNMENT_BAM_ID}_alignment_summary.txt"
echo "Sample,Target,Start,End,Rows,Total_Reads" > "$SUMMARY_FILE"

python3 <<EOF
import csv, subprocess, os, sys

bam_file = "${BAM_FILE}"
sample = "${ALIGNMENT_BAM_ID}"

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
done

# =============================================================================
# Step 2: splice extraction + filtering
# =============================================================================
echo "========== STEP 2a: bamsort_splice ($SAMPLE x targets) =========="

for TARGET_NAME in "${SPLICE_TARGETS_ARR[@]}"; do
    SPLICE_BAM="${SPLICE_INPUT_ROOT}${SAMPLE}${TARGET_NAME}${SPLICE_INPUT_END}"
    SPLICE_CSV="${SPLICE_OUTPUT_DIR}/${SAMPLE}${TARGET_NAME}_splicesites.csv"

    if [ ! -f "$SPLICE_BAM" ]; then
        echo "Warning: splice BAM not found, skipping: $SPLICE_BAM"
        continue
    fi

    echo "Processing: ${SAMPLE}${TARGET_NAME}"
    # Restrict junctions to the viral contig of this run; bamsort_splice.py
    # defaults to mac239, which is correct only for the SIV sample sets.
    python3 "$SPLICE_SCRIPT" "$SPLICE_BAM" "$SPLICE_CSV" \
        --reference "${VIRAL_CONTIG:-mac239}" \
        --min-intron "$SPLICE_MIN_INTRON"
done

if [ -n "${PROCESS_CORRECTION_DIR:-}" ]; then
    echo "========== STEP 2b: bamsort_splice_sites filter ($SAMPLE) =========="
    mkdir -p "$PROCESS_CORRECTION_DIR"
    # SPLICE_SITES picks the reference splice-site coordinates the junctions
    # are labelled against (mac239 or HXB2).
    python3 "$SPLICE_FILTER_SCRIPT" \
        "$SPLICE_OUTPUT_DIR" \
        "$PROCESS_CORRECTION_DIR" \
        --samples "$SAMPLE" \
        --targets "${SPLICE_TARGETS_ARR[@]}" \
        --filter-donor "$SPLICE_FILTER_DONOR" \
        --sites "${SPLICE_SITES:-mac239}"

    # bamsort_splice_sites.py writes one CSV per target for each of its two
    # outputs; this gathers each into one file per sample (header once).
    #   _filtered: the --filter-donor junctions, which process_bamsort in
    #     stepThree reads as <sample>_splicesites_filtered.csv
    #   _all: every junction, kept as a record of the splicing seen
    for KIND in filtered all; do
        COMBINED_FILE="${PROCESS_CORRECTION_DIR%/}/${SAMPLE}_splicesites_${KIND}.csv"
        : > "$COMBINED_FILE"
        wrote_header=0
        for TARGET_NAME in "${SPLICE_TARGETS_ARR[@]}"; do
            TARGET_FILE="${PROCESS_CORRECTION_DIR%/}/${SAMPLE}${TARGET_NAME}_splicesites_${KIND}.csv"
            [ -f "$TARGET_FILE" ] || continue
            if [ "$wrote_header" -eq 0 ]; then
                cat "$TARGET_FILE" >> "$COMBINED_FILE"
                wrote_header=1
            else
                tail -n +2 "$TARGET_FILE" >> "$COMBINED_FILE"
            fi
        done
        if [ "$wrote_header" -eq 0 ]; then
            echo "Warning: no per-target ${KIND} splice files for $SAMPLE; removing empty $COMBINED_FILE"
            rm -f "$COMBINED_FILE"
        else
            echo "Combined ${KIND} splice junctions: $COMBINED_FILE ($(($(wc -l < "$COMBINED_FILE") - 1)) records)"
        fi
    done

    echo "Splice correction complete. Filtered: $PROCESS_CORRECTION_DIR"
fi

echo "stepTwo finished for sample: $SAMPLE"
echo "Run stepThree_process.sh once every sample's stepTwo job has finished."
