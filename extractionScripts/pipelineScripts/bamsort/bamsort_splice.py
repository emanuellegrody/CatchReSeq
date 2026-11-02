#!/usr/bin/env python3

"""
BAM Splice Site Analysis Tool

Description:
    Extracts splice junction information from BAM files by identifying reads with
    large introns (N ≥1000bp) preceded by substantial alignments (M ≥20bp).

Input:
    - BAM file: Aligned reads
    - Output CSV path: File to write splice junction data

Output CSV Columns:
    - CB: Cell barcode (10x Genomics format)
    - UB: UMI barcode (Unique Molecular Identifier)
    - CIGAR: CIGAR string of the read alignment
    - Read sequence: Full query sequence of the read
    - Sequence after the N: Sequence portion after the splice junction
    - Before_N_M_coords: Genomic coordinates of alignment before splice (chr:start-end)
    - After_N_M_coords: Genomic coordinates of alignment after splice (chr:start-end)
    - Donor: Last nucleotide position of the splice donor site (0-based)
    - Acceptor: First nucleotide position of the splice acceptor site (0-based)
    - Reference: Reference chromosome/contig name

Filtering Criteria:
    - Default: Only processes reads mapped to the --reference contig (default 'mac239')
    - With --include-all-splicing: Processes reads from all references
    - Only includes reads with xf tag = 25 (specific alignment flag)
    - Requires alignment pattern: ≥20M + ≥1000N + M (match-intron-match)

Usage:
    python3 bamsort_splice.py input.bam output.csv [--include-all-splicing]

"""

import pysam
import csv
import sys
import argparse
import os


def find_splice_junctions(read, min_intron_size=1000, min_match_size=20):
    """
    Find splice junctions by looking for N≥min_intron_size, then checking
    for immediately preceding M≥min_match_size.

    Returns list of (seq_after_N, beforeM_coords, afterM_coords) tuples.

    Fixed version with correct position tracking and sequence extraction.
    """
    seq = read.query_sequence
    if seq is None or not read.cigartuples:
        return []

    junctions = []

    # Pre-calculate all positions to avoid tracking bugs
    ref_positions = []
    query_positions = []

    ref_pos = read.reference_start
    query_pos = 0

    # Build position maps for each CIGAR operation
    for op, length in read.cigartuples:
        ref_positions.append(ref_pos)
        query_positions.append(query_pos)

        if op == 0:  # M - match/mismatch
            query_pos += length
            ref_pos += length
        elif op == 1:  # I - insertion (consumes query)
            query_pos += length
        elif op == 2:  # D - deletion (consumes reference)
            ref_pos += length
        elif op == 3:  # N - intron (consumes reference)
            ref_pos += length
        elif op == 4:  # S - soft clip (consumes query)
            query_pos += length
        elif op in (5, 6):  # H, P - no consumption
            pass

    # Now scan for splice junctions with accurate position tracking
    for i, (op, length) in enumerate(read.cigartuples):

        # Found a potential intron
        if op == 3 and length >= min_intron_size:  # N operation ≥ threshold

            # Look immediately before this N for an M operation
            if i > 0:
                prev_op, prev_length = read.cigartuples[i - 1]

                if prev_op == 0 and prev_length >= min_match_size:  # M operation ≥ threshold

                    # Calculate coordinates for the M before N
                    before_m_start = ref_positions[i] - prev_length
                    before_m_end = ref_positions[i]
                    beforeM_coords = (read.reference_name, before_m_start, before_m_end)

                    # Calculate coordinates for M after N (if it exists)
                    after_intron_ref_pos = ref_positions[i] + length
                    afterM_coords = None
                    seq_after_N = ""

                    if i + 1 < len(read.cigartuples):
                        next_op, next_length = read.cigartuples[i + 1]
                        if next_op == 0:  # Next operation is M
                            afterM_coords = (
                            read.reference_name, after_intron_ref_pos, after_intron_ref_pos + next_length)

                            # Extract sequence after the intron (fixed)
                            query_start_after_intron = query_positions[i + 1]
                            seq_after_N = seq[query_start_after_intron:query_start_after_intron + next_length]

                    junctions.append((seq_after_N, beforeM_coords, afterM_coords))

    return junctions


def validate_inputs(bam_file, output_csv):
    """Validate input files and paths."""
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM file not found: {bam_file}")

    if not os.access(bam_file, os.R_OK):
        raise PermissionError(f"Cannot read BAM file: {bam_file}")

    # Check if output directory is writable
    output_dir = os.path.dirname(output_csv) or '.'
    if not os.access(output_dir, os.W_OK):
        raise PermissionError(f"Cannot write to output directory: {output_dir}")


def main(bam_file, output_csv, include_all_splicing=False,
         min_intron_size=1000, min_match_size=20, reference="mac239"):
    # Validate inputs
    validate_inputs(bam_file, output_csv)

    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        raise RuntimeError(f"Failed to open BAM file: {e}")

    # Track statistics
    total_reads = 0
    reads_processed = 0
    reads_wrong_ref = 0
    reads_wrong_xf = 0
    reads_no_splice = 0
    junctions_found = 0

    try:
        with open(output_csv, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([
                "CB", "UB", "CIGAR", "Read sequence",
                "Sequence after the N",
                "Before_N_M_coords", "After_N_M_coords",
                "Donor", "Acceptor", "Reference"
            ])

            for read in bam:
                total_reads += 1

                if read.is_unmapped or read.cigartuples is None:
                    continue

                reads_processed += 1

                # Reference filtering
                if not include_all_splicing:
                    if read.reference_name != reference:
                        reads_wrong_ref += 1
                        continue

                # Check xf tag
                if not (read.has_tag("xf") and read.get_tag("xf") == 25):
                    reads_wrong_xf += 1
                    continue

                # Find all splice junctions in this read
                junctions = find_splice_junctions(read, min_intron_size, min_match_size)

                if junctions:
                    # Extract read-level information
                    cb = read.get_tag("CB") if read.has_tag("CB") else ""
                    ub = read.get_tag("UB") if read.has_tag("UB") else ""
                    cigarstring = read.cigarstring
                    sequence = read.query_sequence or ""

                    # Process each junction found in this read
                    for seq_after_N, beforeM_coords, afterM_coords in junctions:
                        before_str = f"{beforeM_coords[0]}:{beforeM_coords[1]}-{beforeM_coords[2]}" if beforeM_coords else ""
                        after_str = f"{afterM_coords[0]}:{afterM_coords[1]}-{afterM_coords[2]}" if afterM_coords else ""

                        # Extract donor and acceptor positions
                        donor = beforeM_coords[2] - 1 if beforeM_coords else ""  # Last base of before_M (0-based)
                        acceptor = afterM_coords[1] if afterM_coords else ""  # First base of after_M

                        writer.writerow([cb, ub, cigarstring, sequence, seq_after_N,
                                         before_str, after_str, donor, acceptor, read.reference_name])
                        junctions_found += 1
                else:
                    reads_no_splice += 1

    finally:
        bam.close()

    # Print statistics
    print(f"Processing complete!")
    print(f"Total reads examined: {total_reads}")
    print(f"Reads processed (mapped with CIGAR): {reads_processed}")
    if not include_all_splicing:
        print(f"Reads from non-{reference} references: {reads_wrong_ref}")
    print(f"Reads without xf:i:25 tag: {reads_wrong_xf}")
    print(f"Reads without qualifying splice pattern: {reads_no_splice}")
    print(f"Splice junctions found: {junctions_found}")
    print(f"Parameters: min_intron={min_intron_size}bp, min_match={min_match_size}bp")
    print(f"Output file: {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract splice junctions using reverse lookup (N-first) method",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default: N≥1000bp, M≥20bp, mac239 only
  python3 bamsort_splice.py input.bam output.csv

  # Restrict to a different viral contig (e.g. the HIV reference)
  python3 bamsort_splice.py input.bam output.csv --reference HXB2

  # Custom thresholds
  python3 bamsort_splice.py input.bam output.csv --min-intron 500 --min-match 15

  # Include all references
  python3 bamsort_splice.py input.bam output.csv --include-all-splicing
        """
    )

    parser.add_argument("bam_file", help="Path to input BAM file")
    parser.add_argument("output_csv", help="Path to output CSV file")
    parser.add_argument("--include-all-splicing", action="store_true",
                        help="Include splice junctions from all references (default: --reference only)")
    parser.add_argument("--reference", default="mac239",
                        help="Contig to restrict splice junctions to (default: mac239)")
    parser.add_argument("--min-intron", type=int, default=1000,
                        help="Minimum intron size in bp (default: 1000)")
    parser.add_argument("--min-match", type=int, default=20,
                        help="Minimum match size before intron in bp (default: 20)")

    args = parser.parse_args()

    print(f"Input BAM: {args.bam_file}")
    print(f"Output CSV: {args.output_csv}")
    print(f"Reference filtering: {'All references' if args.include_all_splicing else args.reference + ' only'}")
    print(f"Minimum intron size: {args.min_intron}bp")
    print(f"Minimum match size: {args.min_match}bp")
    print("-" * 60)

    try:
        main(args.bam_file, args.output_csv, args.include_all_splicing,
             args.min_intron, args.min_match, args.reference)
    except (FileNotFoundError, PermissionError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)