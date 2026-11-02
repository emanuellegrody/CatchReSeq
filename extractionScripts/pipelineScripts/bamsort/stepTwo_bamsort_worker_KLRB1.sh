#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=4G
#SBATCH -t 0:30:00
#SBATCH --verbose

# Per-sample worker (KLRB1): whole-gene bamsort alignment.
# Submitted by Figure1/stepTwo_bamsort_KLRB1.sh; not run directly.
#
# KLRB1 is a host gene captured across its whole locus, so this runs
# bamsort_alignment_reads.py and emits one cellID/UMI/count row per barcode.
# That is a different question from the host worker
# (stepTwo_bamsort_worker_host.sh), which classifies reads at a single exon as
# inclusion or exclusion for the Figure5 splice targets, and from the viral
# worker (stepTwo_bamsort_worker.sh), which follows the alignment step with
# splice-junction extraction. Neither applies here.
#
# Each row of the samples CSV carries its own BAM and output location: the TenX
# and Ca*Re arms of the KLRB1 comparison live in separate trees, so there is no
# shared INPUT_ROOT to derive them from.
#
# There is no stepThree for this arm. Figure1.Rmd collapses the UMI table and
# applies the detection threshold inline, in the chunks under "# KLRB1".

set -euo pipefail

REQUIRED_VARS=(
    SAMPLE BAM_FILE OUTPUT_PREFIX OUTPUT_ROOT COORDS_CSV ALIGNMENT_SCRIPT
)
for v in "${REQUIRED_VARS[@]}"; do
    if [ -z "${!v:-}" ]; then
        echo "Error: required env var $v not set (use stepTwo_bamsort_KLRB1.sh wrapper)"
        exit 1
    fi
done

module purge
source activate cellranger
module load python
module load samtools

for f in "$COORDS_CSV" "$ALIGNMENT_SCRIPT" "$BAM_FILE"; do
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

mkdir -p "$OUTPUT_ROOT"

echo "========== bamsort_alignment_reads (KLRB1: $SAMPLE) =========="
echo "BAM:    $BAM_FILE"
echo "Output: ${OUTPUT_ROOT}${OUTPUT_PREFIX}_bamsort_alignment_<target>.csv"

SUMMARY_FILE="${OUTPUT_ROOT}${OUTPUT_PREFIX}_alignment_summary.txt"
echo "Sample,Target,Chr,Start,End,Rows,Total_Reads" > "$SUMMARY_FILE"

python3 <<PYEOF
import csv, subprocess, os, sys

bam_file = "${BAM_FILE}"
prefix   = "${OUTPUT_PREFIX}"

coords = []
with open("${COORDS_CSV}", "r") as f:
    rows = (line for line in f if line.strip() and not line.lstrip().startswith("#"))
    reader = csv.DictReader(rows)
    for row in reader:
        coords.append({
            "target": row["Target"].strip(),
            "chr":    row["Chr"].strip(),
            "start":  row["Start"].strip(),
            "end":    row["End"].strip(),
        })

print(f"Found {len(coords)} coordinate range(s)")

for coord in coords:
    target, chrom, start, end = coord["target"], coord["chr"], coord["start"], coord["end"]
    output_file = f"${OUTPUT_ROOT}{prefix}_bamsort_alignment_{target}.csv"
    print(f"  Region: {target} ({chrom}:{start}-{end})")

    # No --windows: a single whole-gene window means primary_count would equal
    # count, so there is no boundary read to charge to one region over another.
    cmd = ["python3", "${ALIGNMENT_SCRIPT}", bam_file, chrom, start, end, output_file]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"    ERROR: {e.stderr}")
        with open("${SUMMARY_FILE}", "a") as sf:
            sf.write(f"{prefix},{target},{chrom},{start},{end},ERROR,ERROR\n")
        sys.exit(1)

    if not os.path.exists(output_file):
        print("    ERROR: output not created")
        with open("${SUMMARY_FILE}", "a") as sf:
            sf.write(f"{prefix},{target},{chrom},{start},{end},ERROR,ERROR\n")
        sys.exit(1)

    with open(output_file, "r") as f:
        rows_n = sum(1 for _ in f) - 1
    total_reads = 0
    with open(output_file, "r") as f:
        reader = csv.reader(f); next(reader)
        for row in reader:
            if len(row) > 2:
                try: total_reads += int(row[2])
                except (ValueError, IndexError): pass

    with open("${SUMMARY_FILE}", "a") as sf:
        sf.write(f"{prefix},{target},{chrom},{start},{end},{rows_n},{total_reads}\n")
    print(f"    Success: {rows_n} unique barcode-UMI rows, {total_reads} reads")
PYEOF

echo "stepTwo finished for KLRB1 sample: $SAMPLE"
