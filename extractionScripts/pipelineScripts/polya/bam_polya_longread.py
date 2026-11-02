#!/usr/bin/env python3
"""
Long-read (ONT cDNA) internal-polyA classification for SIVmac239.

Long-read counterpart to bam_polya_virus.py. The question is the same one
Figure 2 Panels F/G ask of short reads: if the polyT primer anneals to genomic
A-rich tracts rather than to the polyA tail, molecules terminate at internal
tracts instead of the 3'LTR cleavage site, which would explain both the viral
read coverage away from the 3' end and the fractured long-read molecules.

Each molecule contributes one priming site, the aligned end carrying the
untemplated polyA/polyT clip. That is reference_end for a reverse-strand read
(the expected orientation for a polyT-primed 10x cDNA read; 98% of primary
alignments here) and reference_start for a forward one. It is the same
polyT-anchor side as Panel G's Read1 anchor, so the Read1 classifier applies:
nearest tract start in either direction, within D bp.

Two control end sets are classified alongside it, because the internal
percentage on its own rises with how loosely a tract is defined:
  other   - the molecule's opposite end, set by where reverse transcription
            stopped rather than by a primer, so it carries no priming signal
  random  - positions drawn uniformly along the full-length transcript

Coordinates: the long-read BAM is aligned to a 10279bp mac239 reference
carrying complete LTRs at both ends, while Panels F/G use the 9602bp reference,
which is exactly the full-length unspliced transcript. The relation is
10279-coordinate = 9602-coordinate + 553, with the 141bp R region duplicated
between the two LTRs. The 10279 FASTA is not kept on disk, so it is rebuilt
here from the 9602 one; the reconstruction reproduces every MD-tag reference
base in the BAM (556,777 bases, 0 mismatches).

Usage:
    python3 bam_polya_longread.py BAM VIRAL_FASTA_9602 --out-dir DIR \
        [--seurat-cellids CSV]
"""

import argparse
import collections
import csv
import os

import numpy as np
import pysam

from polya_utils import find_polyA_in_Seq

# Amplicon products, as Figure 2 "Long-read reported values" defines them: a
# 103bp gag fragment and a product starting 1bp before the env primer. Both
# have primer-defined ends, so neither reports on where priming happened.
AMPLICON_STARTS = {'amplicon_gag': 1440, 'amplicon_envprimer': 8797}
AMPLICON_TOLERANCE = 6

# Offset and R-region repeat length between the two mac239 references.
COORDINATE_OFFSET = 553
R_REGION_LENGTH = 141

# Terminal polyA cleavage sites: Panel G's (141, 9602) carried over by +553.
# 10155 is the 3'LTR cleavage site; 694 is its 5'LTR mirror, which minimap2 can
# pick for a read landing in the duplicated R region.
TERMINAL_SITES = (694, 10155)

# Full-length unspliced transcript, in 10279 coordinates: the 9602bp reference
# start to the 3'LTR cleavage site. The random control is drawn from this span.
TRANSCRIPT_START, TRANSCRIPT_END = COORDINATE_OFFSET, 10155

MIN_LENGTHS = (6, 7, 8, 9, 10)
MAX_NONAS = (0, 1, 2)
MIN_DISTS = (50, 100, 150)

# Reference tract definition for the per-molecule table. 10A_1M is deliberately
# not used: mac239 contains exactly one such tract, so every rate computed
# against it collapses.
REFERENCE_A, REFERENCE_M, REFERENCE_D = 8, 1, 100

N_RANDOM = 20000
RANDOM_SEED = 0


def load_sequence(fasta_path):
    """Read a single-contig FASTA into one uppercase string."""
    with open(fasta_path) as handle:
        return ''.join(line.strip() for line in handle
                       if not line.startswith('>')).upper()


def build_longread_reference(sequence_9602):
    """
    Rebuild the 10279bp reference from the 9602bp one.

    The 9602bp reference spans the transcript, so it starts at the 5'LTR R
    region and ends at the 3'LTR cleavage site. The longer reference adds the
    U3 upstream of that start and the U5 downstream of that end; both are
    present inside the 9602bp sequence because the LTRs are duplicated.
    """
    if len(sequence_9602) != 9602:
        raise ValueError(f'expected a 9602bp reference, got {len(sequence_9602)}bp')
    upstream_u3 = sequence_9602[8908:9461]
    downstream_u5 = sequence_9602[R_REGION_LENGTH:R_REGION_LENGTH + 124]
    reference = upstream_u3 + sequence_9602 + downstream_u5
    if len(reference) != 10279:
        raise ValueError(f'reconstruction is {len(reference)}bp, expected 10279bp')
    return reference


def verify_reference(reference, bam_path, contig):
    """
    Check the reconstruction against the MD tags of every aligned read. Returns
    (matching bases, mismatching bases, positions covered).
    """
    match = mismatch = 0
    covered = set()
    with pysam.AlignmentFile(bam_path) as bam:
        for read in bam:
            if read.is_unmapped or read.reference_name != contig:
                continue
            try:
                pairs = read.get_aligned_pairs(matches_only=True, with_seq=True)
            except ValueError:
                return None
            for _, refpos, refbase in pairs:
                if refbase is None:
                    continue
                covered.add(refpos)
                if reference[refpos] == refbase.upper():
                    match += 1
                else:
                    mismatch += 1
    return match, mismatch, len(covered)


def build_molecules(bam_path, contig='mac239'):
    """
    One molecule per (CB, UMI), following bamsort_isoform_longread.py so the
    molecule set matches cDNA_isoforms.csv: group alignments by read name
    (pooling supplementary records, as --merge-supplementary does), take CB/UMI
    from the primary's tags, and keep the first read seen for a (CB, UMI) pair.
    """
    alignments_by_name = collections.defaultdict(list)
    with pysam.AlignmentFile(bam_path) as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary:
                continue
            if read.reference_name != contig:
                continue
            alignments_by_name[read.query_name].append(read)

    molecules = {}
    n_duplicate = 0
    for alignments in alignments_by_name.values():
        primaries = [r for r in alignments if not r.is_supplementary] or alignments
        representative = primaries[0]
        if not representative.has_tag('CB'):
            continue
        cell_barcode = representative.get_tag('CB')
        if representative.has_tag('UB'):
            umi = representative.get_tag('UB')
        elif representative.has_tag('UR'):
            umi = representative.get_tag('UR')
        else:
            umi = ''
        key = (cell_barcode, umi)
        if key in molecules:
            n_duplicate += 1
            continue

        blocks = []
        for alignment in alignments:
            blocks.extend(alignment.get_blocks())
        blocks.sort()
        ref_start = min(start for start, _ in blocks)
        ref_end = max(end for _, end in blocks)

        group = 'untargeted'
        for name, position in AMPLICON_STARTS.items():
            if abs(ref_start - position) <= AMPLICON_TOLERANCE:
                group = name

        is_reverse = representative.is_reverse
        molecules[key] = dict(
            cellID=cell_barcode,
            UMI=umi,
            group=group,
            is_reverse=is_reverse,
            ref_start=ref_start,
            ref_end=ref_end,
            aligned_bases=sum(end - start for start, end in blocks),
            priming_end=ref_end if is_reverse else ref_start,
            other_end=ref_start if is_reverse else ref_end,
        )
    return list(molecules.values()), n_duplicate


def tract_starts(sequence, min_length, max_nonA):
    """
    Maximal A-rich region starts, using the same find_polyA_in_Seq as
    bam_polya_virus.py and the same nested-interval removal, so tract
    definitions are identical to Panels F/G.
    """
    regions = find_polyA_in_Seq(sequence, min_length, max_nonA)
    kept, max_end = [], -1
    for start, end in sorted(regions, key=lambda region: (region[0], -region[1])):
        if end <= max_end:
            continue
        kept.append(start)
        max_end = end
    return kept


def classify(positions, starts, min_dist):
    """
    Panel G's Read1 classifier: for each position take the nearest tract start
    in either direction, call it unexplained beyond min_dist, terminal if that
    start is a terminal cleavage site, internal otherwise. Returns a label and
    the signed distance per position.
    """
    candidates = np.array(sorted(set(list(starts) + list(TERMINAL_SITES))))
    terminal = set(TERMINAL_SITES)
    labels, distances, nearest = [], [], []
    for position in positions:
        index = np.searchsorted(candidates, position, side='right')
        options = []
        if index < candidates.size:
            options.append(int(candidates[index]))
        if index > 0:
            options.append(int(candidates[index - 1]))
        best = min(options, key=lambda start: abs(start - position))
        distance = abs(best - position)
        if distance > min_dist:
            label = 'unexplained'
        elif best in terminal:
            label = 'terminal'
        else:
            label = 'internal'
        labels.append(label)
        distances.append(distance)
        nearest.append(best)
    return labels, distances, nearest


def main():
    parser = argparse.ArgumentParser(
        description='Classify long-read priming sites as internal or terminal polyA')
    parser.add_argument('bam', help='long-read BAM aligned to the 10279bp mac239')
    parser.add_argument('fasta', help='9602bp mac239 FASTA')
    parser.add_argument('--out-dir', required=True)
    parser.add_argument('--contig', default='mac239')
    parser.add_argument('--seurat-cellids', default=None,
                        help='CSV with a cellID column; molecules outside it are '
                             'dropped, as Panels F/G drop short reads outside the '
                             'Seurat-filtered set')
    parser.add_argument('--group', default='untargeted',
                        choices=['untargeted', 'all'],
                        help='untargeted excludes the two amplicon products')
    parser.add_argument('--skip-verify', action='store_true',
                        help='skip the MD-tag check of the rebuilt reference')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    sequence = build_longread_reference(load_sequence(args.fasta))
    print(f'rebuilt the long-read reference: {len(sequence)}bp')
    if not args.skip_verify:
        result = verify_reference(sequence, args.bam, args.contig)
        if result is None:
            print('  BAM carries no MD tags; reconstruction not verified')
        else:
            match, mismatch, covered = result
            print(f'  MD-tag check: {match} bases match, {mismatch} mismatch, '
                  f'{covered}/{len(sequence)} positions covered')
            if mismatch:
                raise SystemExit('reconstruction disagrees with the BAM; aborting')

    molecules, n_duplicate = build_molecules(args.bam, args.contig)
    print(f'molecules (CB,UMI): {len(molecules)} '
          f'({n_duplicate} duplicate-key reads skipped)')
    print(f'  reverse-strand primary: '
          f'{sum(1 for m in molecules if m["is_reverse"])}')

    if args.seurat_cellids:
        with open(args.seurat_cellids) as handle:
            # Long-read barcodes carry no "-1" suffix, Seurat cellIDs do.
            keep = {row['cellID'].split('-')[0]
                    for row in csv.DictReader(handle)}
        n_before = len(molecules)
        molecules = [m for m in molecules if m['cellID'] in keep]
        print(f'  Seurat cellID filter: {n_before} -> {len(molecules)} '
              f'({100 * len(molecules) / max(n_before, 1):.1f}% kept)')

    print('  group: ' + str(dict(collections.Counter(
        m['group'] for m in molecules))))
    if args.group == 'untargeted':
        molecules = [m for m in molecules if m['group'] == 'untargeted']
    print(f'  classifying {len(molecules)} molecules in '
          f'{len({m["cellID"] for m in molecules})} cells')

    generator = np.random.default_rng(RANDOM_SEED)
    end_sets = {
        'priming': [m['priming_end'] for m in molecules],
        'other': [m['other_end'] for m in molecules],
        'random': [int(p) for p in generator.integers(
            TRANSCRIPT_START, TRANSCRIPT_END + 1, size=N_RANDOM)],
    }

    sweep_rows = []
    reference_starts = None
    for min_length in MIN_LENGTHS:
        for max_nonA in MAX_NONAS:
            starts = tract_starts(sequence, min_length, max_nonA)
            if not starts:
                print(f'  skipping {min_length}A_{max_nonA}M: no tracts in genome')
                continue
            if (min_length, max_nonA) == (REFERENCE_A, REFERENCE_M):
                reference_starts = starts
            for min_dist in MIN_DISTS:
                for end_set, positions in end_sets.items():
                    labels, _, _ = classify(positions, starts, min_dist)
                    counts = collections.Counter(labels)
                    total = len(labels)
                    sweep_rows.append(dict(
                        end_set=end_set,
                        parameters=f'{min_length}A_{max_nonA}M_{min_dist}D',
                        A=min_length, M=max_nonA, D=min_dist,
                        n_tracts=len(starts), n_positions=total,
                        explained_terminal=counts['terminal'],
                        explained_internal=counts['internal'],
                        unexplained=counts['unexplained']))

    sweep_csv = os.path.join(args.out_dir, 'longread_polya_sweep.csv')
    with open(sweep_csv, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sweep_rows[0].keys()))
        writer.writeheader()
        writer.writerows(sweep_rows)
    print(f'wrote {sweep_csv} ({len(sweep_rows)} rows)')

    if reference_starts is None:
        raise SystemExit(f'{REFERENCE_A}A_{REFERENCE_M}M yields no tracts')

    # Per-molecule table for the distance panel. Distances are unbounded here so
    # the plot can show the full distribution; the label applies REFERENCE_D.
    labels, priming_distances, nearest = classify(
        [m['priming_end'] for m in molecules], reference_starts, REFERENCE_D)
    _, other_distances, _ = classify(
        [m['other_end'] for m in molecules], reference_starts, 10 ** 9)
    _, priming_distances_full, _ = classify(
        [m['priming_end'] for m in molecules], reference_starts, 10 ** 9)

    molecules_csv = os.path.join(args.out_dir, 'longread_polya_molecules.csv')
    with open(molecules_csv, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['cellID', 'UMI', 'group', 'is_reverse', 'ref_start',
                         'ref_end', 'aligned_bases', 'priming_end', 'other_end',
                         'nearest_tract', 'dist_priming', 'dist_other',
                         'label_priming'])
        for molecule, label, nearest_tract, dist_priming, dist_other in zip(
                molecules, labels, nearest, priming_distances_full,
                other_distances):
            writer.writerow([
                molecule['cellID'], molecule['UMI'], molecule['group'],
                int(molecule['is_reverse']), molecule['ref_start'],
                molecule['ref_end'], molecule['aligned_bases'],
                molecule['priming_end'], molecule['other_end'],
                nearest_tract, dist_priming, dist_other, label])
    print(f'wrote {molecules_csv} ({len(molecules)} molecules)')

    # Random-position distances go in their own file: there are N_RANDOM of them
    # and they carry no molecule identity.
    _, random_distances, _ = classify(
        end_sets['random'], reference_starts, 10 ** 9)
    random_csv = os.path.join(args.out_dir, 'longread_polya_random.csv')
    with open(random_csv, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['position', 'dist_random'])
        writer.writerows(zip(end_sets['random'], random_distances))
    print(f'wrote {random_csv} ({N_RANDOM} positions)')

    reference_key = f'{REFERENCE_A}A_{REFERENCE_M}M_{REFERENCE_D}D'
    print(f'\nreference setting {reference_key} '
          f'({len(reference_starts)} tracts):')
    for end_set in ('priming', 'other', 'random'):
        row = next(r for r in sweep_rows
                   if r['end_set'] == end_set and r['parameters'] == reference_key)
        total = row['n_positions']
        print(f'  {end_set:>8}: internal {row["explained_internal"]:5d} '
              f'({100 * row["explained_internal"] / total:5.1f}%)  '
              f'terminal {row["explained_terminal"]:5d} '
              f'({100 * row["explained_terminal"] / total:5.1f}%)  '
              f'unexplained {row["unexplained"]:5d} '
              f'({100 * row["unexplained"] / total:5.1f}%)')
    print(f'  median distance to nearest tract: '
          f'priming {int(np.median(priming_distances_full))}bp, '
          f'other {int(np.median(other_distances))}bp, '
          f'random {int(np.median(random_distances))}bp')


if __name__ == '__main__':
    main()
