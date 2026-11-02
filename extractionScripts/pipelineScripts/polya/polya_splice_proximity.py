#!/usr/bin/env python3
"""
Internal polyA sites relative to viral splice sites (HXB2 and SIVmac239).

Question: for a lentiviral genome, how many internal A-rich tracts (candidate
internal-priming sites) fall 200-500 bp UPSTREAM of any splice donor, versus
200-500 bp DOWNSTREAM of any splice acceptor?

The window is the CaRe/CatchReSeq read-reach window: a read priming at an
internal polyA extends ~200-500 bp back toward the 5' side, so a tract sitting
200-500 bp downstream of an acceptor can produce a spliced-looking read across
that junction, while a tract 200-500 bp upstream of a donor cannot.

PolyA tracts come from polya_utils.find_polyA_in_Seq (the same maximal A-rich
region finder used by polya_comparative_stats.py and bam_polya_host.py), swept
over min_length x max_nonA. Tracts within --terminal-window bp of the genome 3'
end are dropped as the terminal polyA tail, so every counted tract is internal.

Splice-site coordinates are 1-based in each virus's own reference numbering.
mac239 values match bamsort_splice_sites.py; HXB2 values are the canonical
K03455 sites (Purcell & Martin 1993).

Outputs (one row per genome x min_length x max_nonA):
    <out>/polya_splice_proximity_summary.tsv  counts + density per window class
    <out>/polya_splice_proximity_sites.tsv    every internal tract, annotated
    <out>/polya_splice_proximity_persite.tsv  counts broken out per donor/acceptor

Usage:
    python3 polya_splice_proximity.py -out OUTDIR
    python3 polya_splice_proximity.py -out OUTDIR --window 200 500 \
        --lengths 6 8 10 --mismatches 0 1 2 --anchor start
"""

import argparse
import itertools
import os
import sys

import pandas as pd

import polya_utils as pu

# -----------------------------------------------------------------------------
# Reference genomes and their canonical splice sites (1-based, own numbering)
# -----------------------------------------------------------------------------
# Override the genome root to run against a different copy of the references:
#   GENOME_ROOT=/my/genomes python3 polya_splice_proximity.py ...
GENOME_ROOT = os.environ.get('GENOME_ROOT', '/projects/b1042/GoyalLab/egrody/genomes')

HXB2_FASTA = os.path.join(GENOME_ROOT, 'GRCh38_HXB2',
                          'GRCh38_HXB2_inputs', 'HXB2.fa')
MAC239_FASTA = os.path.join(GENOME_ROOT, 'Mmul_10_mac239',
                            'inputs', 'mac239.fa')

GENOMES = {
    'HXB2': {
        'fasta': HXB2_FASTA,
        # Purcell & Martin 1993 / LANL HIV-1 K03455 numbering
        'donors': {'D1': 743, 'D2': 4962, 'D3': 5463, 'D4': 6045},
        'acceptors': {'A1': 4913, 'A2': 5389, 'A3': 5777, 'A4c': 5936,
                      'A4a': 5954, 'A4b': 5960, 'A5': 5976, 'A7': 8379},
    },
    'SIVmac239': {
        'fasta': MAC239_FASTA,
        # Same values used by bamsort_splice_sites.py
        'donors': {'D1': 431, 'D2': 4730, 'D3': 5217, 'D4': 6043},
        'acceptors': {'A1': 4658, 'A2': 5141, 'A3': 5746, 'A4': 5818,
                      'A5': 5959, 'A7': 8249},
    },
}

TERMINAL_WINDOW = 400  # matches TERMINAL_DIST_BP in bam_polya_virus.py


def load_genome_seq(fasta_path, chrom=None):
    """Return (chrom_name, uppercase sequence) for a single-contig viral FASTA."""
    fa = pu.open_indexed_fasta(fasta_path)
    refs = list(fa.references)
    if chrom is None:
        if len(refs) != 1:
            sys.exit(f"ERROR: {fasta_path} has {len(refs)} contigs; pass an "
                     f"explicit chromosome name")
        chrom = refs[0]
    elif chrom not in refs:
        sys.exit(f"ERROR: contig '{chrom}' not in {fasta_path}")
    return chrom, fa.fetch(chrom).upper()


def splice_dinucleotide(seq, pos_1based, kind):
    """
    Dinucleotide at a splice site under the convention used here: a donor
    position is the last exonic base, so it is followed by GT; an acceptor
    position is the first exonic base, so it is preceded by AG.
    """
    i = pos_1based - 1
    if kind == 'donor':
        return seq[i + 1:i + 3] if 0 <= i + 3 <= len(seq) else 'NA'
    return seq[i - 2:i] if i - 2 >= 0 else 'NA'


def snap_sites(seq, sites, kind, max_shift):
    """
    Put every site on one convention (donor = last exonic base, acceptor =
    first exonic base) by shifting to the nearest GT/AG within max_shift bp.
    The published HXB2 and mac239 coordinate tables differ from each other by
    1-2 bp on this point; snapping removes that drift so the two genomes are
    measured the same way. Sites with no dinucleotide in range are left alone.

    Returns (snapped_sites, shifts) where shifts maps name -> applied offset.
    """
    snapped, shifts = {}, {}
    for name, pos in sites.items():
        best = None
        for shift in sorted(range(-max_shift, max_shift + 1), key=abs):
            cand = pos + shift
            if not 1 <= cand <= len(seq):
                continue
            if splice_dinucleotide(seq, cand, kind) == ('GT' if kind == 'donor'
                                                        else 'AG'):
                best = shift
                break
        snapped[name] = pos + best if best is not None else pos
        shifts[name] = best
    return snapped, shifts


def anchor_position(start_0, end_excl, anchor):
    """1-based representative coordinate for an A-rich tract."""
    if anchor == 'start':
        return start_0 + 1
    if anchor == 'end':
        return end_excl
    return (start_0 + 1 + end_excl) // 2


def window_positions(sites, near, far, direction, genome_len):
    """
    Set of 1-based genome positions lying `near`..`far` bp on `direction`
    ('upstream' / 'downstream') of any site. Used for the union window size,
    which normalizes counts into tracts per kb of eligible sequence.
    """
    covered = set()
    for pos in sites:
        if direction == 'upstream':
            lo, hi = pos - far, pos - near
        else:
            lo, hi = pos + near, pos + far
        lo = max(lo, 1)
        hi = min(hi, genome_len)
        covered.update(range(lo, hi + 1))
    return covered


def matching_sites(polya_pos, sites, near, far, direction):
    """
    Names of sites for which polya_pos sits `near`..`far` bp on `direction`
    of the site. 'upstream' means the tract is 5' of the donor
    (donor - tract in [near, far]); 'downstream' means the tract is 3' of the
    acceptor (tract - acceptor in [near, far]).
    """
    hits = []
    for name, pos in sites.items():
        delta = (pos - polya_pos) if direction == 'upstream' else (polya_pos - pos)
        if near <= delta <= far:
            hits.append((name, delta))
    return hits


def analyze_genome(genome_name, cfg, param_grid, near, far, anchor,
                   terminal_window, snap_shift, chrom=None):
    """Return (summary_rows, site_rows, persite_rows) for one viral genome."""
    chrom, seq = load_genome_seq(cfg['fasta'], chrom)
    genome_len = len(seq)
    donors, donor_shifts = snap_sites(seq, cfg['donors'], 'donor', snap_shift)
    acceptors, acceptor_shifts = snap_sites(seq, cfg['acceptors'], 'acceptor',
                                            snap_shift)

    def describe(sites, shifts, kind):
        return ', '.join(
            f"{n}={p}{'' if shifts[n] in (0, None) else f'[{shifts[n]:+d}]'}"
            f" ({splice_dinucleotide(seq, p, kind)})"
            for n, p in sorted(sites.items(), key=lambda kv: kv[1]))

    print(f"\n{genome_name}: contig '{chrom}', {genome_len:,} bp")
    print(f"  donors:    " + describe(donors, donor_shifts, 'donor'))
    print(f"  acceptors: " + describe(acceptors, acceptor_shifts, 'acceptor'))
    unsnapped = ([n for n, s in donor_shifts.items() if s is None] +
                 [n for n, s in acceptor_shifts.items() if s is None])
    if unsnapped:
        print(f"  NOTE: no GT/AG within {snap_shift} bp for "
              f"{', '.join(unsnapped)}; annotated position used as given")

    up_window = window_positions(donors.values(), near, far,
                                 'upstream', genome_len)
    down_window = window_positions(acceptors.values(), near, far,
                                   'downstream', genome_len)
    terminal_cutoff = genome_len - terminal_window
    up_window = {p for p in up_window if p <= terminal_cutoff}
    down_window = {p for p in down_window if p <= terminal_cutoff}
    both_window = up_window & down_window
    print(f"  window bases (excluding terminal {terminal_window} bp): "
          f"upstream-of-donor {len(up_window):,}, "
          f"downstream-of-acceptor {len(down_window):,}, "
          f"overlapping both {len(both_window):,}")

    summary_rows, site_rows, persite_rows = [], [], []

    for min_len, max_nonA in param_grid:
        regions = pu.deduplicate_regions(
            pu.find_polyA_in_Seq(seq, min_len, max_nonA))

        n_total = len(regions)
        n_terminal = 0
        n_up = n_down = n_both = n_neither = 0
        per_donor = {n: 0 for n in donors}
        per_acceptor = {n: 0 for n in acceptors}

        for start_0, end_excl in regions:
            pos = anchor_position(start_0, end_excl, anchor)
            if pos > terminal_cutoff:
                n_terminal += 1
                continue

            up_hits = matching_sites(pos, donors, near, far, 'upstream')
            down_hits = matching_sites(pos, acceptors, near, far, 'downstream')

            if up_hits and down_hits:
                n_both += 1
            elif up_hits:
                n_up += 1
            elif down_hits:
                n_down += 1
            else:
                n_neither += 1

            for name, _ in up_hits:
                per_donor[name] += 1
            for name, _ in down_hits:
                per_acceptor[name] += 1

            site_rows.append({
                'Genome': genome_name,
                'Min_Length': min_len,
                'Max_Mismatches': max_nonA,
                'PolyA_Start': start_0 + 1,
                'PolyA_End': end_excl,
                'PolyA_Length': end_excl - start_0,
                'Anchor_Position': pos,
                'Sequence': seq[start_0:end_excl],
                'Upstream_Of_Donors': ';'.join(f"{n}:{d}" for n, d in up_hits),
                'Downstream_Of_Acceptors': ';'.join(f"{n}:{d}" for n, d in down_hits),
                'Class': ('both' if (up_hits and down_hits)
                          else 'upstream_of_donor' if up_hits
                          else 'downstream_of_acceptor' if down_hits
                          else 'neither'),
            })

        n_internal = n_total - n_terminal
        # Tracts counted per window include the 'both' overlap, so the two
        # densities stay comparable to their own window sizes.
        up_any = n_up + n_both
        down_any = n_down + n_both
        summary_rows.append({
            'Genome': genome_name,
            'Min_Length': min_len,
            'Max_Mismatches': max_nonA,
            'Genome_Length': genome_len,
            'PolyA_Total': n_total,
            'PolyA_Terminal_Excluded': n_terminal,
            'PolyA_Internal': n_internal,
            'Upstream_Of_Any_Donor': up_any,
            'Downstream_Of_Any_Acceptor': down_any,
            'Both': n_both,
            'Only_Upstream_Of_Donor': n_up,
            'Only_Downstream_Of_Acceptor': n_down,
            'Neither': n_neither,
            'Upstream_Window_bp': len(up_window),
            'Downstream_Window_bp': len(down_window),
            'Upstream_Per_kb': round(up_any / len(up_window) * 1000, 4)
                               if up_window else 0.0,
            'Downstream_Per_kb': round(down_any / len(down_window) * 1000, 4)
                                 if down_window else 0.0,
            'Genomewide_Internal_Per_kb': round(
                n_internal / max(terminal_cutoff, 1) * 1000, 4),
        })

        for name, pos in sorted(donors.items(), key=lambda kv: kv[1]):
            persite_rows.append({
                'Genome': genome_name, 'Min_Length': min_len,
                'Max_Mismatches': max_nonA, 'Site_Type': 'donor',
                'Site': name, 'Site_Position': pos,
                'Window_Start': max(pos - far, 1), 'Window_End': max(pos - near, 1),
                'PolyA_Count': per_donor[name],
            })
        for name, pos in sorted(acceptors.items(), key=lambda kv: kv[1]):
            persite_rows.append({
                'Genome': genome_name, 'Min_Length': min_len,
                'Max_Mismatches': max_nonA, 'Site_Type': 'acceptor',
                'Site': name, 'Site_Position': pos,
                'Window_Start': min(pos + near, genome_len),
                'Window_End': min(pos + far, genome_len),
                'PolyA_Count': per_acceptor[name],
            })

    return summary_rows, site_rows, persite_rows


def main():
    parser = argparse.ArgumentParser(
        description="Count internal polyA tracts upstream of splice donors vs "
                    "downstream of splice acceptors in HXB2 and SIVmac239.")
    parser.add_argument('-out', required=True, help='Output directory')
    parser.add_argument('--genomes', nargs='+', default=list(GENOMES),
                        choices=list(GENOMES),
                        help='Which viral genomes to analyze')
    parser.add_argument('--window', type=int, nargs=2, default=[200, 500],
                        metavar=('NEAR', 'FAR'),
                        help='Distance window from the splice site (default 200 500)')
    parser.add_argument('--lengths', type=int, nargs='+', default=[6, 7, 8, 9, 10],
                        help='Minimum polyA tract lengths to sweep')
    parser.add_argument('--mismatches', type=int, nargs='+', default=[0, 1, 2],
                        help='Maximum non-A bases tolerated within a tract')
    parser.add_argument('--anchor', choices=['start', 'end', 'mid'],
                        default='start',
                        help="Tract coordinate used for distance (default start, "
                             "the 5'-most base of the A-rich tract)")
    parser.add_argument('--terminal-window', type=int, default=TERMINAL_WINDOW,
                        help='Tracts within this many bp of the genome 3\' end '
                             'are the terminal polyA tail, not internal')
    parser.add_argument('--snap-shift', type=int, default=5,
                        help='Shift each annotated splice site up to this many '
                             'bp onto the nearest GT/AG so HXB2 and mac239 use '
                             'one coordinate convention (0 disables)')
    parser.add_argument('--report-length', type=int, default=8,
                        help='Parameter combination echoed in the stdout preview')
    parser.add_argument('--report-mismatches', type=int, default=2,
                        help='Parameter combination echoed in the stdout preview')
    args = parser.parse_args()

    near, far = args.window
    if near > far:
        sys.exit("ERROR: --window NEAR must be <= FAR")

    output_dir = args.out
    os.makedirs(output_dir, exist_ok=True)

    param_grid = list(itertools.product(args.lengths, args.mismatches))
    print(f"Window: {near}-{far} bp   anchor: {args.anchor}   "
          f"{len(param_grid)} polyA parameter combinations")

    summary_rows, site_rows, persite_rows = [], [], []
    for genome_name in args.genomes:
        s, si, ps = analyze_genome(genome_name, GENOMES[genome_name],
                                   param_grid, near, far, args.anchor,
                                   args.terminal_window, args.snap_shift)
        summary_rows += s
        site_rows += si
        persite_rows += ps

    summary = pd.DataFrame(summary_rows)
    sites = pd.DataFrame(site_rows)
    persite = pd.DataFrame(persite_rows)

    summary_file = os.path.join(output_dir, 'polya_splice_proximity_summary.tsv')
    sites_file = os.path.join(output_dir, 'polya_splice_proximity_sites.tsv')
    persite_file = os.path.join(output_dir, 'polya_splice_proximity_persite.tsv')
    summary.to_csv(summary_file, sep='\t', index=False)
    sites.to_csv(sites_file, sep='\t', index=False)
    persite.to_csv(persite_file, sep='\t', index=False)

    print(f"\nWrote {summary_file}")
    print(f"Wrote {sites_file}")
    print(f"Wrote {persite_file}")

    preview = summary[(summary.Min_Length == args.report_length) &
                      (summary.Max_Mismatches == args.report_mismatches)]
    if preview.empty:
        preview = summary
    print(f"\nPreview (min_length={args.report_length}, "
          f"max_mismatches={args.report_mismatches}):")
    print(preview[['Genome', 'Min_Length', 'Max_Mismatches', 'PolyA_Internal',
                   'Upstream_Of_Any_Donor', 'Downstream_Of_Any_Acceptor',
                   'Both', 'Neither', 'Upstream_Per_kb',
                   'Downstream_Per_kb']].to_string(index=False))


if __name__ == '__main__':
    main()
