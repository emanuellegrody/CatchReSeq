#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=4G
#SBATCH -t 1:00:00
#SBATCH --verbose

# Per-sample worker (host targets): bamsort alignment only.
# Submitted by Figure5/stepTwo_bamsort.sh; not run directly.
# All configuration is provided via env vars exported by the wrapper.
#
# Combining the per-target alignment tables into per-UMI counts is
# stepThree_process.sh, which gathers across every sample and so has to wait
# for all of these jobs to finish.

set -euo pipefail

REQUIRED_VARS=(
    SAMPLE
    COORDS_CSV INPUT_ROOT INPUT_END OUTPUT_ROOT ALIGNMENT_SCRIPT
)
for v in "${REQUIRED_VARS[@]}"; do
    if [ -z "${!v:-}" ]; then
        echo "Error: required env var $v not set (use stepTwo_bamsort.sh wrapper)"
        exit 1
    fi
done

module purge
source activate cellranger
module load python
module load samtools

for f in "$COORDS_CSV" "$ALIGNMENT_SCRIPT"; do
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

# =============================================================================
# Alignment (single sample, all coords). Coords CSV must have columns:
#   Target,Chr,Start,End
# =============================================================================
echo "========== bamsort_alignment_reads (host: $SAMPLE) =========="

BAM_FILE="${INPUT_ROOT}${SAMPLE}${INPUT_END}"
if [ ! -f "$BAM_FILE" ]; then
    echo "Error: alignment BAM not found: $BAM_FILE"
    exit 1
fi

SUMMARY_FILE="${OUTPUT_ROOT}${SAMPLE}_alignment_summary.txt"
echo "Sample,Target,Chr,Start,End,Inc_Rows,Inc_Reads,Exc_Rows,Exc_Reads" > "$SUMMARY_FILE"

python3 <<EOF
import csv, subprocess, os, sys

bam_file = "${BAM_FILE}"
sample   = "${SAMPLE}"

coords = []
with open("${COORDS_CSV}", "r") as f:
    # Skip blank lines and #-prefixed comments before/after the header
    rows = (line for line in f if line.strip() and not line.lstrip().startswith("#"))
    reader = csv.DictReader(rows)
    for row in reader:
        coords.append({
            "target": row["Target"].strip(),
            "chr":    row["Chr"].strip(),
            "start":  row["Start"].strip(),
            "end":    row["End"].strip(),
        })

print(f"Found {len(coords)} coordinate ranges")

def csv_rows_and_reads(path):
    if not os.path.exists(path):
        return 0, 0
    with open(path, "r") as f:
        rows = sum(1 for _ in f) - 1
    total = 0
    with open(path, "r") as f:
        reader = csv.reader(f); next(reader)
        for r in reader:
            if len(r) > 2:
                try: total += int(r[2])
                except (ValueError, IndexError): pass
    return rows, total

def strip_probe_version(name):
    """Drop a trailing "v<digits>" probe-version suffix.

    A library may carry a probe-version suffix the coordinate-table target
    name does not: sample "WT_Slkv4" targets the region named "Slk". Names
    that merely end in a digit (BCL2L1) are left alone.
    """
    head = name.rstrip("0123456789")
    if head != name and head.endswith("v"):
        return head[:-1]
    return name


for coord in coords:
    target, chrom, start, end = coord["target"], coord["chr"], coord["start"], coord["end"]
    output_file = f"${OUTPUT_ROOT}{sample}_bamsort_alignment_{target}.csv"
    exclusion_file = f"${OUTPUT_ROOT}{sample}_bamsort_alignment_{target}_exclusion.csv"
    # Strip the optional genotype prefix (e.g. WT_/KO_) and any probe-version
    # suffix (e.g. Slkv4 -> Slk) before matching. Accept both exact matches
    # (target == sample_target) and multi-row gene targets whose name extends
    # the sample (e.g. sample "PTPRC" matches targets "PTPRC_constant",
    # "PTPRC_exon4", ...).
    sample_target = sample.split("_", 1)[1] if "_" in sample else sample
    sample_target = strip_probe_version(sample_target)
    is_match = (target == sample_target) or target.startswith(sample_target + "_")
    label = "Match" if is_match else "Off-target"
    print(f"  Region: {target} ({chrom}:{start}-{end}) [{label}]")

    cmd = ["python3", "${ALIGNMENT_SCRIPT}",
           bam_file, chrom, start, end, output_file,
           "--exclusion-output", exclusion_file]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        if os.path.exists(output_file):
            inc_rows, inc_reads = csv_rows_and_reads(output_file)
            exc_rows, exc_reads = csv_rows_and_reads(exclusion_file)
            with open("${SUMMARY_FILE}", "a") as sf:
                sf.write(f"{sample},{target},{chrom},{start},{end},"
                         f"{inc_rows},{inc_reads},{exc_rows},{exc_reads}\n")

            if is_match:
                print(f"    Success: inclusion={inc_rows} barcodes/{inc_reads} reads, "
                      f"exclusion={exc_rows} barcodes/{exc_reads} reads")
            else:
                print(f"    Off-target: inclusion={inc_rows}/{inc_reads}, "
                      f"exclusion={exc_rows}/{exc_reads} (outputs not saved)")
                # Off-target alignment CSVs are removed by default to avoid
                # cluttering the output directory with non-matching results.
                # Comment out the os.remove() calls below to keep them.
                for p in (output_file, exclusion_file):
                    if os.path.exists(p):
                        os.remove(p)
                stats_log = os.path.join(
                    os.path.dirname(output_file),
                    f"stats_{os.path.splitext(os.path.basename(output_file))[0]}.txt",
                )
                if os.path.exists(stats_log):
                    os.remove(stats_log)
        else:
            print("    ERROR: output not created")
            with open("${SUMMARY_FILE}", "a") as sf:
                sf.write(f"{sample},{target},{chrom},{start},{end},ERROR,ERROR,ERROR,ERROR\n")
            sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"    ERROR: {e.stderr}")
        with open("${SUMMARY_FILE}", "a") as sf:
            sf.write(f"{sample},{target},{chrom},{start},{end},ERROR,ERROR\n")
        sys.exit(1)
EOF

echo "stepTwo finished for sample: $SAMPLE"
echo "Run stepThree_process.sh once every sample's stepTwo job has finished."
