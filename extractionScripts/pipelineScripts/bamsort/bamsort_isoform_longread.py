#!/usr/bin/env python3

"""
BAM Isoform Classifier for Long Read Single-Cell RNA-seq

Description:
    Classifies reads into isoforms (US, MS, SS, or any) based on alignment.

Isoform Classification Logic:
    US (Unspliced): Alignment of ≥50bp within region D1/431-A1/4658, with no splice
                    junction (N) spanning that region (+/- 5bp tolerance)
    MS (Multiply Spliced): D4/6043-A7/8252 spliced out, AND either
                           (a) alignment before D4/6043 + alignment after A7/8252, OR
                           (b) alignment before D1/431 + alignment at/after A7/8252
                               with no alignment in the full D1-A7 span
                               (doubly-spliced: both D1-A1 and D4-A7 removed)
    SS (Singly Spliced): Alignment before D1/431, splice gap D1/431-A1/4658,
                         alignment in D4/6043-A7/8252
    any: All other alignment patterns

Input:
    - BAM file: mac239-aligned long reads with cell barcode (CB) and UMI tags (optional)
    - Output CSV path: File to write isoform classifications

Output CSV Columns:
    - cellID: Cell barcode from CB tag
    - UMI: Unique molecular identifier from UB/UR tag (optional)
    - isoform: Classification (US, MS, SS, or any)

Technical Details:
    - Uses CIGAR string parsing to identify aligned (M) and splice gap (N) regions
    - Tolerance zone: +/- 5bp for all coordinate boundaries
    - Only processes primary mapped reads aligned to mac239
    - Deduplicates output by cellID + UMI (one classification per pair)
    - Prioritizes isoform classification in order: US → MS → SS → any

Usage:
    python3 bamsort_isoform_longread.py input.bam output.csv [--tolerance 10]
"""

import pysam
import csv
import sys
import argparse
import os
from collections import defaultdict


# Define genomic regions with biological meaning
REGION_US_START = 431
REGION_US_END = 4658
REGION_MS_GAP_START = 6043
REGION_MS_GAP_END = 8252
REGION_SS_GAP_START = 431
REGION_SS_GAP_END = 4658
REGION_SS_AFTER_START = 6043
REGION_SS_AFTER_END = 8252

# Technical parameters
MIN_ALIGNMENT_LENGTH = 50  # Minimum bp for US classification


def parse_alignment_regions(read):
    """
    Parse CIGAR string to extract aligned regions (M) and splice gaps (N).

    Returns:
        tuple: (aligned_regions, gap_regions) each as [(start, end), ...]
    """
    if not read.cigartuples or read.reference_start is None or read.is_unmapped:
        return [], []

    aligned_regions = []
    gap_regions = []
    ref_pos = read.reference_start

    for op, length in read.cigartuples:
        if op == 0:  # M - match/mismatch (aligned region)
            aligned_regions.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op == 2:  # D - deletion (consumes reference)
            ref_pos += length
        elif op == 3:  # N - intron/skip (splice junction)
            gap_regions.append((ref_pos, ref_pos + length))
            ref_pos += length
        # I (insertion), S (soft clip), H (hard clip), P (padding) don't consume reference

    return aligned_regions, gap_regions


def has_alignment_in_region(aligned_regions, region_start, region_end,
                            min_length=1, tolerance=10):
    """
    Check if any aligned region overlaps with target region within tolerance.

    Args:
        aligned_regions: List of (start, end) tuples
        region_start: Start of target region
        region_end: End of target region
        min_length: Minimum overlap length required (bp)
        tolerance: +/- bp tolerance for boundaries

    Returns:
        bool: True if sufficient overlap exists
    """
    # Apply tolerance to expand the target region
    region_start_with_tolerance = region_start - tolerance
    region_end_with_tolerance = region_end + tolerance

    for align_start, align_end in aligned_regions:
        # Calculate overlap
        overlap_start = max(align_start, region_start_with_tolerance)
        overlap_end = min(align_end, region_end_with_tolerance)
        overlap_length = max(0, overlap_end - overlap_start)

        if overlap_length >= min_length:
            return True

    return False



def lacks_alignment_in_region(aligned_regions, region_start, region_end,
                              tolerance=10):
    """
    Check that no aligned region overlaps with the core of the target region.

    Uses tolerance to SHRINK the query region (more permissive — small
    alignment fragments near boundaries are ignored).

    Args:
        aligned_regions: List of (start, end) tuples
        region_start: Start of expected gap region
        region_end: End of expected gap region
        tolerance: +/- bp tolerance for boundaries

    Returns:
        bool: True if no alignment exists in the core region
    """
    core_start = region_start + tolerance
    core_end = region_end - tolerance

    if core_start >= core_end:
        return True

    for align_start, align_end in aligned_regions:
        overlap_start = max(align_start, core_start)
        overlap_end = min(align_end, core_end)
        if overlap_end > overlap_start:
            return False

    return True


def classify_isoform(aligned_regions, tolerance=10, offset=0):
    """
    Classify read into isoform based on alignment and splice pattern.

    Classification is mutually exclusive and prioritized as:
    1. US (Unspliced): any alignment within D1/431-A1/4658 — overrides all others
    2. MS (Multiply Spliced): D4/6043-A7/8252 intron removed AND either
       (a) alignment before D4 + alignment after A7, or
       (b) alignment before D1 + alignment at/after A7 only, with the full
           D1-A7 span empty (doubly-spliced — both introns removed)
    3. SS (Singly Spliced): no alignment in D1/431-A1/4658 (intron removed),
       alignment before D1 and in D4/6043-A7/8252, D4-A7 region IS aligned
    4. any: Everything else

    Args:
        aligned_regions: List of (start, end) tuples from M operations
        tolerance: +/- bp tolerance for boundaries
        offset: bp by which the BAM's coordinate system exceeds that of the
            9602bp mac239.fa the REGION_* constants are written in. Read
            coordinates are shifted down by this amount before classification,
            which is equivalent to shifting every donor/acceptor site up.

    Returns:
        str: Isoform classification
    """
    if offset:
        aligned_regions = [(s - offset, e - offset) for s, e in aligned_regions]
    if not aligned_regions:
        return "any"

    # Check US first: any alignment in D1-A1 overrides everything
    if has_alignment_in_region(aligned_regions, REGION_US_START,
                                REGION_US_END,
                                min_length=MIN_ALIGNMENT_LENGTH,
                                tolerance=tolerance):
        return "US"

    d4a7_spliced = lacks_alignment_in_region(aligned_regions,
                                             REGION_MS_GAP_START,
                                             REGION_MS_GAP_END, tolerance)
    d1a1_spliced = lacks_alignment_in_region(aligned_regions,
                                             REGION_SS_GAP_START,
                                             REGION_SS_GAP_END, tolerance)

    # Check MS: D4-A7 intron removed, alignment before D4 and after A7
    if d4a7_spliced:
        has_before = has_alignment_in_region(aligned_regions, 0,
                                             REGION_MS_GAP_START,
                                             min_length=1, tolerance=tolerance)
        has_after = has_alignment_in_region(aligned_regions, REGION_MS_GAP_END,
                                            float('inf'),
                                            min_length=1, tolerance=tolerance)
        if has_before and has_after:
            return "MS"

    # Doubly-spliced MS: alignment before D1 and at/after A7 only, with no
    # alignment anywhere in the D1-A7 span (both D1-A1 and D4-A7 removed).
    # The above d4a7-only check also catches this (since "before D1" ⊂
    # "before D4"), but this branch makes the doubly-spliced pattern
    # explicit.
    d1a7_spliced = lacks_alignment_in_region(aligned_regions,
                                             REGION_US_START,
                                             REGION_MS_GAP_END,
                                             tolerance)
    if d1a7_spliced:
        has_before_d1 = has_alignment_in_region(aligned_regions, 0,
                                                REGION_US_START,
                                                min_length=1,
                                                tolerance=tolerance)
        has_after_a7 = has_alignment_in_region(aligned_regions,
                                               REGION_MS_GAP_END,
                                               float('inf'),
                                               min_length=1,
                                               tolerance=tolerance)
        if has_before_d1 and has_after_a7:
            return "MS"

    # Check SS: D1-A1 intron removed, alignment before D1 and in D4-A7,
    # but D4-A7 must be aligned (otherwise it's MS or any)
    if d1a1_spliced and not d4a7_spliced:
        has_before = has_alignment_in_region(aligned_regions, 0,
                                             REGION_SS_GAP_START,
                                             min_length=1, tolerance=tolerance)
        has_after = has_alignment_in_region(aligned_regions, REGION_SS_AFTER_START,
                                            REGION_SS_AFTER_END,
                                            min_length=1, tolerance=tolerance)
        if has_before and has_after:
            return "SS"

    # Default: any other pattern
    return "any"


def extract_tag_safe(read, *tag_names):
    """
    Safely extract tags from read, trying multiple tag names.

    Returns:
        tuple: Values for requested tags (None if not found)
    """
    values = []
    for tag_name in tag_names:
        try:
            values.append(read.get_tag(tag_name))
        except KeyError:
            values.append(None)
    return tuple(values)


def process_bam_file(bam_file, output_csv, tolerance=5, reference="mac239",
                     split_bam=False, offset=0, merge_supplementary=False):
    """
    Process BAM file and classify reads into isoforms.

    Only processes primary alignments to the specified reference. Output is
    deduplicated by (cellID, UMI) — one classification per pair.

    Args:
        bam_file: Path to input BAM file
        output_csv: Path to output CSV file
        tolerance: +/- bp tolerance for region boundaries
        reference: Reference name to filter for (default: mac239)
        split_bam: If True, write separate BAM files for each isoform class
    """
    # Validate inputs
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM file not found: {bam_file}")

    if not os.access(bam_file, os.R_OK):
        raise PermissionError(f"Cannot read BAM file: {bam_file}")

    # Check output directory is writable
    output_dir = os.path.dirname(output_csv) or '.'
    if not os.access(output_dir, os.W_OK):
        raise PermissionError(f"Cannot write to output directory: {output_dir}")

    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        raise RuntimeError(f"Failed to open BAM file: {e}")

    # Open split BAM writers if requested
    bam_writers = {}
    if split_bam:
        output_base = os.path.splitext(output_csv)[0]
        for iso in ['US', 'SS', 'MS', 'any']:
            bam_path = f"{output_base}_{iso}.bam"
            bam_writers[iso] = pysam.AlignmentFile(bam_path, "wb", header=bam.header)

    # Statistics tracking
    stats = {
        'total_reads': 0,
        'unmapped_reads': 0,
        'secondary_supplementary': 0,
        'wrong_reference': 0,
        'no_cell_barcode': 0,
        'no_umi': 0,
        'duplicate_cb_umi': 0,
        'classified_reads': 0,
        'isoform_counts': defaultdict(int)
    }

    seen_cb_umi = {}  # (cb, umi) -> dict of per-molecule fields

    try:
        # Group alignments by read name. minimap2 represents some spliced long
        # reads as a primary plus one or more supplementary alignments; with
        # merge_supplementary the blocks of all of a read's alignments are
        # pooled, so a split-mapped junction is not lost.
        by_name = defaultdict(list)
        for read in bam:
            stats['total_reads'] += 1

            if read.is_unmapped:
                stats['unmapped_reads'] += 1
                continue

            if read.is_secondary or (read.is_supplementary
                                     and not merge_supplementary):
                stats['secondary_supplementary'] += 1
                continue

            if read.reference_name != reference:
                stats['wrong_reference'] += 1
                continue

            by_name[read.query_name].append(read)

        for query_name, reads in by_name.items():
            # CB/UMI and the representative record come from the primary
            # alignment; supplementary records carry the same tags but a
            # primary is always present for a mapped read.
            primary = [r for r in reads if not r.is_supplementary] or reads
            rep = primary[0]

            cb, ub, ur = extract_tag_safe(rep, "CB", "UB", "UR")
            umi = ub if ub is not None else (ur if ur is not None else "")

            if cb is None:
                stats['no_cell_barcode'] += 1
                continue

            if not umi:
                stats['no_umi'] += 1

            # Deduplicate by (CB, UMI) — keep first classification
            key = (cb, umi)
            if key in seen_cb_umi:
                stats['duplicate_cb_umi'] += 1
                continue

            aligned_regions = []
            gap_regions = []
            for r in reads:
                ar, gr = parse_alignment_regions(r)
                aligned_regions.extend(ar)
                gap_regions.extend(gr)
            aligned_regions.sort()

            isoform = classify_isoform(aligned_regions, tolerance=tolerance,
                                       offset=offset)

            # ref_start lets downstream analysis group molecules by alignment
            # start position; n_junctions counts splice gaps, the number a
            # molecule can report on. Both are in BAM coordinates.
            seen_cb_umi[key] = {
                'isoform': isoform,
                'ref_start': min(s for s, e in aligned_regions),
                'ref_end': max(e for s, e in aligned_regions),
                'aligned_bases': sum(e - s for s, e in aligned_regions),
                'n_junctions': sum(1 for s, e in gap_regions if e - s >= 50),
            }
            stats['classified_reads'] += 1
            stats['isoform_counts'][isoform] += 1

            if split_bam:
                bam_writers[isoform].write(rep)

        bam.close()

        # Close and index split BAMs
        if split_bam:
            output_base = os.path.splitext(output_csv)[0]
            for iso, writer in bam_writers.items():
                writer.close()
                bam_path = f"{output_base}_{iso}.bam"
                tmp_path = bam_path + ".unsorted.tmp"
                os.rename(bam_path, tmp_path)
                pysam.sort("-o", bam_path, tmp_path)
                os.remove(tmp_path)
                pysam.index(bam_path)

        # Write deduplicated results
        with open(output_csv, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['cellID', 'UMI', 'isoform', 'ref_start',
                             'ref_end', 'aligned_bases', 'n_junctions'])
            for (cb, umi), m in seen_cb_umi.items():
                writer.writerow([cb, str(umi), m['isoform'], m['ref_start'],
                                 m['ref_end'], m['aligned_bases'],
                                 m['n_junctions']])

    except Exception as e:
        if 'bam' in locals():
            bam.close()
        for w in bam_writers.values():
            w.close()
        raise RuntimeError(f"Error during processing: {e}")

    return stats


def print_statistics(stats, tolerance, offset=0, merge_supplementary=False):
    """
    Print processing statistics to console.
    """
    print("Processing complete!")
    print("-" * 60)
    print(f"Total reads examined: {stats['total_reads']}")
    print(f"Unmapped reads: {stats['unmapped_reads']}")
    print(f"Secondary/supplementary (skipped): {stats['secondary_supplementary']}")
    print(f"Wrong reference (skipped): {stats['wrong_reference']}")
    print(f"Reads without cell barcode (skipped): {stats['no_cell_barcode']}")
    print(f"Reads without UMI (still classified): {stats['no_umi']}")
    print(f"Duplicate CB+UMI (skipped): {stats['duplicate_cb_umi']}")
    print(f"Unique CB+UMI pairs classified: {stats['classified_reads']}")
    print("-" * 60)
    print("Isoform distribution:")
    for isoform in ['US', 'MS', 'SS', 'any']:
        count = stats['isoform_counts'][isoform]
        pct = 100 * count / max(1, stats['classified_reads'])
        print(f"  {isoform}: {count} ({pct:.1f}%)")

    # Paper values: a molecule is "resolved" when it carries a single-isoform
    # call (US, SS or MS); "any" means its alignment pattern is compatible with
    # more than one isoform. Percentages are over all classified molecules and
    # over resolved molecules respectively.
    n_total = stats['classified_reads']
    n_US = stats['isoform_counts']['US']
    n_SS = stats['isoform_counts']['SS']
    n_MS = stats['isoform_counts']['MS']
    n_resolved = n_US + n_SS + n_MS
    print("-" * 60)
    print(f"Resolved to a single isoform: {n_resolved}/{n_total} = "
          f"{100 * n_resolved / max(1, n_total):.1f}%")
    print(f"Of resolved, unspliced (US):  {n_US}/{n_resolved} = "
          f"{100 * n_US / max(1, n_resolved):.1f}%")
    print(f"Of resolved, spliced (SS+MS): {n_SS + n_MS}/{n_resolved} = "
          f"{100 * (n_SS + n_MS) / max(1, n_resolved):.1f}%")
    print("-" * 60)
    print(f"Tolerance used: +/- {tolerance}bp")
    print(f"Coordinate offset applied: {offset}bp")
    print(f"Supplementary alignments merged: {merge_supplementary}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Classify long-read RNA-seq isoforms based on alignment patterns",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Isoform Definitions:
  US (Unspliced):        ≥50bp alignment within 431-4658, no splice gap there
  MS (Multiply Spliced): Splice gap 6043-8252 with alignment before 6043 and
                         after 8252 — OR — doubly-spliced (alignment in 0-431
                         and 8252+ only, full 431-8252 span empty)
  SS (Singly Spliced):   Splice gap 431-4658, alignment before 431 and in 6043-8252
  any:                   All other patterns

Tolerance:
  All coordinate boundaries use +/- 5bp tolerance by default to account for
  biological variation and alignment uncertainty.

Examples:
  # Default tolerance (10bp)
  python3 bamsort_isoform_classifier.py input.bam output.csv

  # Custom tolerance
  python3 bamsort_isoform_classifier.py input.bam output.csv --tolerance 10

  # Strict boundaries (no tolerance)
  python3 bamsort_isoform_classifier.py input.bam output.csv --tolerance 0
        """
    )

    parser.add_argument("bam_file", help="Path to input BAM file")
    parser.add_argument("output_csv", help="Path to output CSV file")
    parser.add_argument("--tolerance", type=int, default=10,
                        help="Coordinate tolerance in bp (default: 10)")
    parser.add_argument("--split-bam", action="store_true",
                        help="Write separate BAM files for each isoform class "
                             "(US, SS, MS, any)")
    parser.add_argument("--coordinate-offset", type=int, default=0,
                        help="bp by which the BAM reference exceeds the 9602bp "
                             "mac239.fa these donor/acceptor coordinates are "
                             "written in (default: 0). The EGS007 long-read "
                             "BAM is aligned to a 10279bp mac239 and needs 553")
    parser.add_argument("--merge-supplementary", action="store_true",
                        help="Pool the aligned blocks of each read's primary "
                             "and supplementary alignments before classifying, "
                             "so split-mapped splice junctions are retained")

    args = parser.parse_args()

    try:
        stats = process_bam_file(args.bam_file, args.output_csv, args.tolerance,
                                 split_bam=args.split_bam,
                                 offset=args.coordinate_offset,
                                 merge_supplementary=args.merge_supplementary)
        print_statistics(stats, args.tolerance,
                         offset=args.coordinate_offset,
                         merge_supplementary=args.merge_supplementary)

    except (FileNotFoundError, PermissionError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        sys.exit(1)