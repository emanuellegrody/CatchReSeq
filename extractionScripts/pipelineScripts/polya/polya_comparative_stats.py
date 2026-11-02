#!/usr/bin/env python3
"""
Comparative PolyA statistics generator: target host genes vs. viral genome.
Sweeps minimum length and mismatch tolerance, writes one TSV row per
(genome, min_length, max_mismatches) combination.
"""

import argparse
import itertools
import os
import sys
from collections import defaultdict

import pandas as pd

import polya_utils as pu


def scan_sequence(sequence, param_grid, stats, entity):
    """Run all (min_length, max_nonA) combos on one sequence and accumulate."""
    for (min_len, max_nonA) in param_grid:
        regions = pu.deduplicate_regions(
            pu.find_polyA_in_Seq(sequence, min_len, max_nonA))
        stats[(min_len, max_nonA)][entity]['count'] += len(regions)
        stats[(min_len, max_nonA)][entity]['bases'] += sum(e - s for s, e in regions)
    return len(sequence)


def main():
    parser = argparse.ArgumentParser(description="Generate comparative PolyA statistics.")
    parser.add_argument('-host', required=True, help='Host genome FASTA file')
    parser.add_argument('-viral', required=True, help='Viral genome FASTA file')
    parser.add_argument('-gtf', required=True, help='Host GTF file')
    parser.add_argument('-genes', required=True, help='Text file with target host gene names')
    parser.add_argument('-out', required=True, help='Output directory')
    parser.add_argument('--lengths', type=int, nargs='+', default=[6, 7, 8, 9, 10])
    parser.add_argument('--mismatches', type=int, nargs='+', default=[0, 1, 2])
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    param_grid = list(itertools.product(args.lengths, args.mismatches))
    stats = defaultdict(lambda: {'host':  {'count': 0, 'bases': 0},
                                 'viral': {'count': 0, 'bases': 0}})
    total_host_bases = 0
    total_viral_bases = 0

    print("Loading target genes...")
    target_genes = pu.load_target_genes(args.genes)
    print(f"  {len(target_genes):,} target gene names")

    print("\nProcessing host genome...")
    processed_genes = 0
    for gene_name, seq, _, _, _ in pu.iter_target_gene_seqs(
            args.host, args.gtf, target_genes, forward_only=True):
        total_host_bases += scan_sequence(seq, param_grid, stats, 'host')
        processed_genes += 1
    print(f"  Scanned {processed_genes:,} forward-strand target genes "
          f"({total_host_bases:,} bases)")

    print("\nProcessing viral genome...")
    viral_fa = pu.open_indexed_fasta(args.viral)
    for chrom in viral_fa.references:
        total_viral_bases += scan_sequence(
            viral_fa.fetch(chrom).upper(), param_grid, stats, 'viral')
    print(f"  Scanned viral genome ({total_viral_bases:,} bases)")

    print("\nCompiling statistics...")
    rows = []
    for (min_l, max_m) in param_grid:
        for entity in ('host', 'viral'):
            s = stats[(min_l, max_m)][entity]
            total_b = total_host_bases if entity == 'host' else total_viral_bases
            if total_b == 0:
                print(f"  Warning: total bases for {entity} is 0", file=sys.stderr)
                continue
            count = s['count']
            poly_bases = s['bases']
            rows.append({
                'Genome': entity.capitalize(),
                'Min_Length': min_l,
                'Max_Mismatches': max_m,
                'Total_Genome_Bases_Scanned': total_b,
                'PolyA_Count': count,
                'Total_PolyA_Bases': poly_bases,
                'PolyA_Seqs_per_1000_bp': round(count / total_b * 1000, 4),
                'Percent_Genome_Covered': round(poly_bases / total_b * 100, 6),
                'Avg_PolyA_Length': round(poly_bases / count, 2) if count else 0,
            })

    df = pd.DataFrame(rows)
    out_file = os.path.join(args.out, 'polya_comparative_stats.tsv')
    df.to_csv(out_file, sep='\t', index=False)
    print(f"\nWrote {out_file}")
    print("\nPreview:")
    print(df[['Genome', 'Min_Length', 'Max_Mismatches',
              'PolyA_Count', 'PolyA_Seqs_per_1000_bp']].to_string(index=False))


if __name__ == '__main__':
    main()
