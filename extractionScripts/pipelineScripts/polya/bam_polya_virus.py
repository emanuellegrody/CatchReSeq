#!/usr/bin/env python3
"""
SIV BAM PolyA structure analyzer (processing only).

PolyA-tract evidence here comes from the read sequences themselves: for each
Read1 with a quality alignment we extract the cell barcode (last 16 bp) and
UMI (preceding 12 bp), find the longest polyA stretch immediately upstream,
classify the read as 'terminal' (Read1 maps within +/-400 bp of the SIVmac239
polyA tail at coord 9601) or 'undetermined', and record the paired Read2
position and Read1<->Read2 distance.

When --cellranger-bam is supplied, the Cell Ranger possorted BAM's mac239
chromosome is pre-scanned for {read_name -> (xf, CB, UB)}; pairs are then
restricted to xf in {17, 25} and the CSV is annotated with corrected CB/UB
columns matching Seurat cellID format (16 bp + '-1' suffix).

Output: <prefix>_polya_structure_analysis.csv consumed by Figure2.Rmd Panel D.

Usage:
    python bam_polya_virus.py input.bam output_prefix [--sample-size N]
                                                      [--min-polya-length N]
                                                      [--min-a-fraction F]
                                                      [--cellranger-bam PATH]
                                                      [--cellranger-chrom CHR]
"""

import argparse
import os
import re
import sys
import warnings
from collections import defaultdict

import pandas as pd
import pysam

warnings.filterwarnings('ignore')

CIGAR_OP_RE = re.compile(r'(\d+)([MIDNSHP=X])')
# Terminal-polyA priming-site coordinate in the mac239 reference at
# /projects/b1042/GoyalLab/egrody/genomes/mac239/fasta/mac239.fa (length 9602).
# 3'LTR PAS (AATAAAG) at 9579-9585; cleavage ~15-30 nt downstream coincides
# with the genome end. Used only by the legacy Python summary line; Panel D
# uses the R-side classifier in Figure2.Rmd which sets terminal_site itself.
SIV_TERMINAL_COORD = 9602
TERMINAL_DIST_BP = 400
DEFAULT_XF_KEEP = (17, 25)


def load_cellranger_lookup(cr_bam_path, chrom='mac239'):
    """Build {query_name -> (xf, CB, UB)} from one chromosome of the Cell Ranger
    possorted BAM. Reads without CB are skipped to save memory; xf=None is
    stored for reads with no xf tag, so the caller can decide whether to keep
    them."""
    print(f"Loading Cell Ranger tags from {cr_bam_path} (chrom={chrom})...")
    lookup = {}
    with pysam.AlignmentFile(cr_bam_path, 'rb') as bam:
        if chrom not in bam.references:
            sys.exit(f"ERROR: chromosome '{chrom}' not in Cell Ranger BAM "
                     f"references")
        for r in bam.fetch(chrom):
            if r.is_secondary or r.is_supplementary:
                continue
            if not r.has_tag('CB'):
                continue
            lookup[r.query_name] = (r.get_tag('xf') if r.has_tag('xf') else None,
                                    r.get_tag('CB'),
                                    r.get_tag('UB') if r.has_tag('UB') else None)
    print(f"  Loaded {len(lookup):,} reads with CB tags from {chrom}")
    return lookup


def extract_10x_barcode_umi(sequence):
    """Extract 16 bp cell barcode and 12 bp UMI from the 3' end of Read 1."""
    if not sequence or len(sequence) < 28:
        return None, None
    return sequence[-16:], sequence[-28:-16]


def find_polya_stretch_from_end(sequence, min_length=6, min_a_fraction=2/3):
    """Longest A-rich stretch starting from the end. Returns
    (start, end, length, a_fraction) or None."""
    if not sequence:
        return None
    sequence = sequence.upper()
    seq_len = len(sequence)
    best = None
    best_len = 0
    for plen in range(min_length, seq_len + 1):
        start = seq_len - plen
        stretch = sequence[start:seq_len]
        a_fraction = stretch.count('A') / plen
        if a_fraction >= min_a_fraction and plen > best_len:
            best = (start, seq_len, plen, a_fraction)
            best_len = plen
    return best


def is_quality_read(cigarstring):
    """First two M-ops in the CIGAR sum to >=90 bp."""
    m_lengths = [int(L) for L, op in CIGAR_OP_RE.findall(cigarstring) if op == 'M']
    return sum(m_lengths[:2]) >= 90


def classify_by_pos(read1_pos, terminal_coord=SIV_TERMINAL_COORD,
                    dist_from_coord=TERMINAL_DIST_BP):
    """Terminal if Read1 leftmost pos is within `dist_from_coord` of the
    SIV polyA tail coordinate."""
    if abs(int(read1_pos) - terminal_coord) <= dist_from_coord:
        return 'terminal'
    return 'undetermined'


def analyze_read_structure(bam_file, sample_size=None,
                           min_polya_length=6, min_a_fraction=2/3,
                           cr_lookup=None, xf_keep=None):
    """Stream BAM, pair Read1+Read2 by name, classify, return per-pair DataFrame.
    If cr_lookup is provided ({query_name -> (xf, CB, UB)}), pairs are kept only
    when the read has a valid CB tag AND its xf is in xf_keep (set of ints).
    `None` is allowed as an xf value if it's in xf_keep (i.e., reads with no xf
    tag at all can be retained by including None in xf_keep)."""
    print(f"Analyzing polyA structure in: {bam_file}")
    print("Read1 layout: [genomic][polyA][UMI:12bp][cellID:16bp]")
    print(f"Terminal classification: |Read1 pos - {SIV_TERMINAL_COORD}| <= {TERMINAL_DIST_BP} bp")
    if cr_lookup is not None:
        print(f"xf filter ACTIVE: keeping pairs with xf in "
              f"{sorted(v for v in xf_keep if v is not None)}"
              f"{' (incl. reads with no xf tag)' if None in xf_keep else ''} "
              f"(lookup size = {len(cr_lookup):,})")
    print()

    read_pairs = defaultdict(dict)
    total_reads = 0

    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for read in bam:
            total_reads += 1
            if total_reads % 500_000 == 0:
                complete = sum(1 for p in read_pairs.values() if 'r1' in p and 'r2' in p)
                print(f"  Processed {total_reads:,} reads, {complete:,} complete pairs")
            if read.is_unmapped or not read.is_paired or not read.query_sequence:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
            info = {
                'chr': read.reference_name,
                'pos': read.reference_start,
                'end': read.reference_end,
                'is_reverse': read.is_reverse,
                'cigar': read.cigarstring,
                'mapq': read.mapping_quality,
                'sequence': read.query_sequence,
            }
            if read.is_read1:
                read_pairs[read.query_name]['r1'] = info
            elif read.is_read2:
                read_pairs[read.query_name]['r2'] = info
            if sample_size and len(read_pairs) >= sample_size * 2:
                break

    print(f"Done streaming BAM. Total reads: {total_reads:,}")

    pairs = list(read_pairs.items())
    if sample_size:
        pairs = pairs[:sample_size]

    rows = []
    skipped_no_cr = skipped_xf = 0
    for name, pair in pairs:
        if 'r1' not in pair or 'r2' not in pair:
            continue
        r1, r2 = pair['r1'], pair['r2']
        cb, umi = extract_10x_barcode_umi(r1['sequence'])
        if not cb or not umi:
            continue
        if not is_quality_read(r1['cigar']):
            continue

        cb_corrected = ub_corrected = xf_val = None
        if cr_lookup is not None:
            entry = cr_lookup.get(name)
            if entry is None:
                skipped_no_cr += 1
                continue
            xf_val, cb_corrected, ub_corrected = entry
            if xf_val not in xf_keep:
                skipped_xf += 1
                continue

        pre_barcode = r1['sequence'][:-28]
        polya = find_polya_stretch_from_end(pre_barcode, min_polya_length, min_a_fraction)
        if polya:
            polya_start, polya_end, polya_length, a_fraction = polya
            genomic = pre_barcode[:polya_start]
        else:
            polya_start = polya_end = polya_length = a_fraction = None
            genomic = pre_barcode

        classification = classify_by_pos(r1['pos'])

        read_distance = None
        if r1['chr'] == r2['chr'] and r1['pos'] is not None and r2['pos'] is not None:
            read_distance = abs(r2['pos'] - r1['pos'])

        rows.append({
            'read_name': name,
            'cell_barcode': cb,
            'umi': umi,
            'cb_corrected': cb_corrected,
            'ub_corrected': ub_corrected,
            'xf': xf_val,
            'read1_pos': r1['pos'],
            'read1_end': r1['end'],
            'read1_is_reverse': r1['is_reverse'],
            'read1_mapq': r1['mapq'],
            'read2_pos': r2['pos'],
            'read2_end': r2['end'],
            'read2_is_reverse': r2['is_reverse'],
            'read2_mapq': r2['mapq'],
            'polya_found': polya is not None,
            'polya_start_in_sequence': polya_start,
            'polya_end_in_sequence': polya_end,
            'polya_length': polya_length,
            'polya_a_fraction': a_fraction,
            'genomic_sequence_length': len(genomic),
            'genomic_sequence': genomic,
            'genomic_classification': classification,
            'read1_to_read2_distance': read_distance,
            'r1_pre_barcode_sequence': pre_barcode,
            'r2_sequence': r2['sequence'],
        })

    if cr_lookup is not None:
        xf_display = sorted(v for v in xf_keep if v is not None)
        if None in xf_keep:
            xf_display = xf_display + ['no-tag']
        print(f"  Pairs dropped (not in Cell Ranger lookup): {skipped_no_cr:,}")
        print(f"  Pairs dropped (xf not in {xf_display}):       {skipped_xf:,}")
    print(f"Analyzed {len(rows):,} pairs")
    return pd.DataFrame(rows)


def print_summary(df):
    """Print a concise summary of the analysis."""
    n = len(df)
    if n == 0:
        print("No pairs to summarize")
        return
    polya = df['polya_found'].sum()
    terminal = (df['genomic_classification'] == 'terminal').sum()
    print(f"\n{'='*60}\nSUMMARY\n{'='*60}")
    print(f"Total pairs:       {n:,}")
    print(f"PolyA found:       {polya:,} ({polya/n*100:.1f}%)")
    print(f"Terminal:          {terminal:,} ({terminal/n*100:.1f}%)")
    print(f"Undetermined:      {n-terminal:,} ({(n-terminal)/n*100:.1f}%)")
    print(f"Unique barcodes:   {df['cell_barcode'].nunique():,}")
    print(f"Unique UMIs:       {df['umi'].nunique():,}")
    if polya > 0:
        polya_df = df[df['polya_found']]
        print(f"PolyA length mean: {polya_df['polya_length'].mean():.1f} bp")
        print(f"A fraction mean:   {polya_df['polya_a_fraction'].mean():.3f}")


def main():
    parser = argparse.ArgumentParser(description='BAM polyA structure analyzer (processing only)')
    parser.add_argument('bam_file', help='STAR-aligned BAM file')
    parser.add_argument('output_prefix', help='Output prefix (writes <prefix>_polya_structure_analysis.csv)')
    parser.add_argument('--sample-size', type=int, help='First N pairs only (for testing)')
    parser.add_argument('--min-polya-length', type=int, default=6)
    parser.add_argument('--min-a-fraction', type=float, default=0.67)
    parser.add_argument('--cellranger-bam',
                        help='Cell Ranger possorted_genome_bam.bam (must be '
                             'indexed) to source xf/CB/UB tags. When provided, '
                             'pairs are filtered to xf in --xf-allow.')
    parser.add_argument('--cellranger-chrom', default='mac239',
                        help='Chromosome in the Cell Ranger BAM to scan for the '
                             'xf/CB/UB lookup (default: mac239)')
    parser.add_argument('--xf-allow', default=','.join(str(v) for v in DEFAULT_XF_KEEP),
                        help='Comma-separated xf values to retain (default: '
                             '"17,25" = Cell Ranger "confidently counted" reads). '
                             'For mac239 you typically want "0,17,25" so LTR '
                             'multi-mappers (xf=0) are kept. Use "none" as a '
                             'value to also keep reads with no xf tag.')
    args = parser.parse_args()

    if not os.path.exists(args.bam_file):
        sys.exit(f"ERROR: BAM not found: {args.bam_file}")

    xf_keep = set()
    for tok in args.xf_allow.split(','):
        tok = tok.strip().lower()
        if tok == 'none':
            xf_keep.add(None)
        elif tok:
            try:
                xf_keep.add(int(tok))
            except ValueError:
                sys.exit(f"ERROR: --xf-allow value '{tok}' is not an integer or 'none'")
    if not xf_keep:
        sys.exit("ERROR: --xf-allow produced an empty set")

    os.makedirs(os.path.dirname(args.output_prefix) or '.', exist_ok=True)

    cr_lookup = None
    if args.cellranger_bam:
        if not os.path.exists(args.cellranger_bam):
            sys.exit(f"ERROR: Cell Ranger BAM not found: {args.cellranger_bam}")
        cr_lookup = load_cellranger_lookup(args.cellranger_bam,
                                           args.cellranger_chrom)

    df = analyze_read_structure(args.bam_file, args.sample_size,
                                args.min_polya_length, args.min_a_fraction,
                                cr_lookup=cr_lookup, xf_keep=xf_keep)
    if df.empty:
        sys.exit("ERROR: No pairs survived filtering")

    out_csv = f"{args.output_prefix}_polya_structure_analysis.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv}")
    print_summary(df)


if __name__ == '__main__':
    main()
