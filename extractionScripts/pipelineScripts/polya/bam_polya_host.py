#!/usr/bin/env python3
"""
Host BAM dual-anchor polyA classification (Read1 vs Read2).

Operates on the STAR paired-end alignment of host reads against the joint
Mmul_10+mac239 genome. For each Read1 record on a forward-strand target
protein-coding gene:

- read1_end  = priming site (Read1 aligns reverse-complemented on + strand
               polyT-captured fragments, so reference_end is the polyT anchor).
- read2_pos  = mate's reference_start (upstream of priming by ~insert length).

For every (min_length A, max_nonA M, min_dist D) combination, classify each
pair using BOTH anchors:

- Read1 anchor: closest polyA-tract start (or gene_end terminal site) within
                D bp of read1_end IN EITHER DIRECTION  ->  terminal/internal.
- Read2 anchor: closest polyA-tract start (or gene_end) DOWNSTREAM of
                read2_pos within D bp                   ->  terminal/internal.

When --cellranger-bam is supplied, each per-chromosome worker pre-scans the
Cell Ranger possorted BAM's matching chromosome and keeps only reads whose
query_name carries xf in {17, 25}.

Output: tidy CSV with rows {anchor, min_length, max_nonA, min_dist,
terminal, internal, unexplained, total_reads}, parallelized across chromosomes.
"""

import argparse
import bisect
import csv
import itertools
import multiprocessing as mp
import os
import sys
import time
from functools import partial

import csv

import pysam

import polya_utils as pu


XF_KEEP = {17, 25}


def load_seurat_cellids(path):
    """Load Seurat cellID whitelist from a CSV with a 'cellID' column (the
    format written by UMAPcoords/*UMAP_coords.csv). Returns a set of strings
    matching the Cell Ranger `CB` tag format (16 bp + '-1' suffix)."""
    cellids = set()
    with open(path, newline='') as f:
        rd = csv.DictReader(f)
        if 'cellID' not in rd.fieldnames:
            sys.exit(f"ERROR: --seurat-cellids file {path} has no 'cellID' "
                     f"column (found {rd.fieldnames})")
        for row in rd:
            cb = row['cellID'].strip()
            if cb:
                cellids.add(cb)
    return cellids


def load_cellranger_tag_lookup(cr_bam_path, chrom, seurat_cellids=None):
    """Return {query_name -> (CB, UB)} for reads on `chrom` of the Cell Ranger
    BAM whose xf tag is in XF_KEEP. If seurat_cellids is given, also require CB
    to be in that whitelist. UB defaults to None if absent. Empty dict if chrom
    is absent."""
    lookup = {}
    with pysam.AlignmentFile(cr_bam_path, 'rb') as bam:
        if chrom not in bam.references:
            return lookup
        for r in bam.fetch(chrom):
            if r.is_secondary or r.is_supplementary:
                continue
            if not r.has_tag('xf') or r.get_tag('xf') not in XF_KEEP:
                continue
            if not r.has_tag('CB'):
                continue
            cb = r.get_tag('CB')
            if seurat_cellids is not None and cb not in seurat_cellids:
                continue
            ub = r.get_tag('UB') if r.has_tag('UB') else None
            lookup[r.query_name] = (cb, ub)
    return lookup


# ---------------------------------------------------------------------------
# Reference build
# ---------------------------------------------------------------------------

def build_gene_lookups(host_fasta, gtf, target_genes, lengths, mismatches):
    """
    Returns (gene_table, polya_lookups) where:
      gene_table   = {chrom: [(gene_start_0, gene_end, gene_name), ...] sorted by gene_start_0}
      polya_lookups = {(A, M): {gene_name: sorted_list_of_genome_polyA_starts_plus_gene_end}}
                      (gene_end is always appended as the "terminal site" for that gene)
    Forward-strand target genes only; identical to the convention used by
    polya_comparative_stats.py and fasta_polya_host.py.
    """
    fasta = pu.open_indexed_fasta(host_fasta)
    db    = pu.open_or_build_gtf_db(gtf)
    fasta_chroms = set(fasta.references)

    # First pass: collect (chrom, start_0, end, gene_name, sequence) for every
    # forward-strand target gene. Re-using iter_target_gene_seqs would cost a
    # second FASTA fetch later -- do it once here.
    gene_records = []  # list of (chrom, start_0, end, gene_name, seq)
    for gene in db.features_of_type('gene'):
        gene_name = gene.attributes.get('gene_name', [gene.id])[0]
        if gene_name not in target_genes or gene.strand == '-':
            continue
        chrom = gene.chrom
        if chrom not in fasta_chroms:
            alt = chrom[3:] if chrom.startswith('chr') else f"chr{chrom}"
            if alt not in fasta_chroms:
                continue
            chrom = alt
        start_0 = gene.start - 1
        end     = gene.end
        try:
            seq = fasta.fetch(chrom, start_0, end).upper()
        except KeyError:
            continue
        gene_records.append((chrom, start_0, end, gene_name, seq))

    # Build per-chrom sorted gene tables (drop seq) for read-to-gene lookup.
    gene_table = {}
    for chrom, start_0, end, gene_name, _seq in gene_records:
        gene_table.setdefault(chrom, []).append((start_0, end, gene_name))
    for chrom in gene_table:
        gene_table[chrom].sort(key=lambda r: r[0])

    # Build per-(A, M) per-gene polyA-start lists, with gene_end appended.
    polya_lookups = {}
    for (A, M) in itertools.product(lengths, mismatches):
        per_gene = {}
        for chrom, start_0, end, gene_name, seq in gene_records:
            regions = pu.deduplicate_regions(pu.find_polyA_in_Seq(seq, A, M))
            starts = sorted(r[0] + start_0 for r in regions)
            # Terminal site is gene_end; append once and keep sorted.
            starts.append(end)
            starts.sort()
            per_gene[gene_name] = (starts, end)  # (sorted starts incl. terminal, terminal_site)
        polya_lookups[(A, M)] = per_gene
        print(f"  A={A} M={M}: {len(per_gene):,} genes "
              f"(mean polyA tracts/gene: "
              f"{sum(len(v[0]) for v in per_gene.values()) / max(len(per_gene), 1):.1f})",
              flush=True)

    return gene_table, polya_lookups


# ---------------------------------------------------------------------------
# Per-pair classification
# ---------------------------------------------------------------------------

def classify_r1(r1_end, sorted_starts, terminal_site, D):
    """Closest start within D in either direction. Returns 'terminal',
    'internal', or 'unexplained'."""
    idx = bisect.bisect_left(sorted_starts, r1_end)
    candidates = []
    if idx < len(sorted_starts):
        candidates.append(sorted_starts[idx])
    if idx > 0:
        candidates.append(sorted_starts[idx - 1])
    if not candidates:
        return 'unexplained'
    closest = min(candidates, key=lambda s: abs(s - r1_end))
    if abs(closest - r1_end) > D:
        return 'unexplained'
    return 'terminal' if closest == terminal_site else 'internal'


def classify_r2(r2_pos, sorted_starts, terminal_site, D):
    """Closest start strictly DOWNSTREAM of r2_pos within D bp."""
    idx = bisect.bisect_right(sorted_starts, r2_pos)
    if idx >= len(sorted_starts):
        return 'unexplained'
    closest = sorted_starts[idx]
    if closest - r2_pos > D:
        return 'unexplained'
    return 'terminal' if closest == terminal_site else 'internal'


def find_gene_at(gene_list, pos):
    """Return gene_name whose [start_0, end) contains pos, or None.
    gene_list is sorted by start_0; uses bisect. Assumes non-overlapping or
    keeps the first hit."""
    if not gene_list:
        return None
    starts = [g[0] for g in gene_list]
    idx = bisect.bisect_right(starts, pos) - 1
    if idx < 0:
        return None
    s, e, name = gene_list[idx]
    if s <= pos < e:
        return name
    return None


# ---------------------------------------------------------------------------
# Worker
# ---------------------------------------------------------------------------

def _read_passes_static_filters(read):
    """Paired, primary, mapped, R1 reverse, same-chrom mate."""
    if not read.is_paired or read.is_unmapped or read.mate_is_unmapped:
        return False
    if read.is_secondary or read.is_supplementary:
        return False
    if not read.is_read1 or not read.is_reverse:
        return False
    if read.next_reference_id != read.reference_id:
        return False
    return True


def classify_chromosome(chrom, bam_file, gene_table, polya_lookups, dist_grid,
                        cr_bam_file=None, seurat_cellids=None,
                        min_umi_reads=1):
    """Per-anchor counts. Returns dict keyed by (anchor, A, M, D).
    If cr_bam_file is provided, pre-scan that chrom for query_names whose xf is
    in XF_KEEP (and, if seurat_cellids is given, whose CB is in the whitelist);
    skip any read whose name isn't in the lookup. If min_umi_reads > 1, a first
    pass over the STAR BAM counts reads per (CB, UB) and only UMIs with at
    least `min_umi_reads` surviving reads are kept for the classification pass."""
    anchors = ('Read1', 'Read2')
    out = {(anchor, A, M, D): {'terminal': 0, 'internal': 0,
                               'unexplained': 0, 'total': 0}
           for anchor in anchors
           for (A, M) in polya_lookups
           for D in dist_grid}

    cr_lookup = None
    if cr_bam_file is not None:
        cr_lookup = load_cellranger_tag_lookup(cr_bam_file, chrom,
                                               seurat_cellids)
        if not cr_lookup:
            # No xf/CB-passing reads on this chrom in Cell Ranger BAM -> nothing to count
            return out

    bam = pysam.AlignmentFile(bam_file, 'rb')
    try:
        if chrom not in bam.references:
            return out
        gene_list = gene_table.get(chrom, [])
        if not gene_list:
            return out
        am_pairs = list(polya_lookups.keys())

        # Pass 1 (only when UMI filter active): count reads per (CB, UB) that
        # pass all per-read filters AND land within a target gene. Singletons
        # are dropped on the classification pass.
        umi_keep = None
        if cr_lookup is not None and min_umi_reads > 1:
            umi_counts = {}
            for read in bam.fetch(chrom):
                if not _read_passes_static_filters(read):
                    continue
                entry = cr_lookup.get(read.query_name)
                if entry is None:
                    continue
                r1_end = read.reference_end
                if r1_end is None:
                    continue
                if find_gene_at(gene_list, r1_end) is None:
                    continue
                umi_counts[entry] = umi_counts.get(entry, 0) + 1
            umi_keep = {k for k, v in umi_counts.items() if v >= min_umi_reads}
            if not umi_keep:
                return out

        # Pass 2: classify
        for read in bam.fetch(chrom):
            if not _read_passes_static_filters(read):
                continue
            if cr_lookup is not None:
                entry = cr_lookup.get(read.query_name)
                if entry is None:
                    continue
                if umi_keep is not None and entry not in umi_keep:
                    continue
            r1_end = read.reference_end
            r2_pos = read.next_reference_start
            if r1_end is None or r2_pos is None:
                continue
            gene_name = find_gene_at(gene_list, r1_end)
            if gene_name is None:
                continue

            for (A, M) in am_pairs:
                gene_polya = polya_lookups[(A, M)].get(gene_name)
                if gene_polya is None:
                    continue
                sorted_starts, terminal_site = gene_polya
                for D in dist_grid:
                    cls1 = classify_r1(r1_end, sorted_starts, terminal_site, D)
                    out[('Read1', A, M, D)]['total']  += 1
                    out[('Read1', A, M, D)][cls1]     += 1
                    cls2 = classify_r2(r2_pos, sorted_starts, terminal_site, D)
                    out[('Read2', A, M, D)]['total']  += 1
                    out[('Read2', A, M, D)][cls2]     += 1
    finally:
        bam.close()
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Host BAM dual-anchor polyA classification.")
    parser.add_argument('bam_file', help='STAR paired-end indexed BAM (.bai required)')
    parser.add_argument('gene_file', help='Target host gene name list')
    parser.add_argument('genome_file', help='Host genome FASTA (indexed)')
    parser.add_argument('gtf_file', help='Host GTF')
    parser.add_argument('--out', required=True, help='Output CSV path')
    parser.add_argument('--lengths', type=int, nargs='+', default=[6, 7, 8, 9, 10])
    parser.add_argument('--mismatches', type=int, nargs='+', default=[0, 1, 2])
    parser.add_argument('--dists-r1', type=int, nargs='+',
                        default=[50, 100, 150, 200])
    parser.add_argument('--dists-r2', type=int, nargs='+',
                        default=[250, 400, 500, 600])
    parser.add_argument('--cellranger-bam',
                        help='Cell Ranger possorted_genome_bam.bam (must be '
                             'indexed) to source xf/CB/UB tags. When provided, '
                             f'reads are kept only when xf in {sorted(XF_KEEP)}.')
    parser.add_argument('--seurat-cellids',
                        help='CSV with a "cellID" column listing the Seurat-'
                             'filtered cell barcodes (CB format, e.g. '
                             'XXXXXXXXXXXXXXXX-1). Requires --cellranger-bam.')
    parser.add_argument('--min-umi-reads', type=int, default=1,
                        help='Drop UMIs (CB, UB) with fewer than this many '
                             'surviving reads on a given chrom. Requires '
                             '--cellranger-bam. Default 1 = no UMI filter.')
    slurm_cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', mp.cpu_count()))
    parser.add_argument('--cpus', type=int, default=slurm_cpus)
    args = parser.parse_args()

    if not os.path.exists(f"{args.bam_file}.bai"):
        sys.exit(f"ERROR: BAM index missing. Run 'samtools index {args.bam_file}'.")
    if args.cellranger_bam:
        if not os.path.exists(args.cellranger_bam):
            sys.exit(f"ERROR: Cell Ranger BAM not found: {args.cellranger_bam}")
        if not os.path.exists(f"{args.cellranger_bam}.bai"):
            sys.exit(f"ERROR: Cell Ranger BAM index missing: "
                     f"{args.cellranger_bam}.bai")
    if (args.seurat_cellids or args.min_umi_reads > 1) and not args.cellranger_bam:
        sys.exit("ERROR: --seurat-cellids and --min-umi-reads > 1 require "
                 "--cellranger-bam (we need CB/UB tags from the Cell Ranger BAM).")

    seurat_cellids = None
    if args.seurat_cellids:
        if not os.path.exists(args.seurat_cellids):
            sys.exit(f"ERROR: --seurat-cellids file not found: {args.seurat_cellids}")
        seurat_cellids = load_seurat_cellids(args.seurat_cellids)
        print(f"Loaded {len(seurat_cellids):,} Seurat cellIDs from {args.seurat_cellids}",
              flush=True)

    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)
    t0 = time.perf_counter()

    target_genes = pu.load_target_genes(args.gene_file)
    print(f"Loaded {len(target_genes):,} target genes", flush=True)

    # Single dist grid passed to the worker; we keep R1 and R2 separated by
    # *anchor*, not D, so we run the union here and slice per-anchor in output.
    dist_grid = sorted(set(args.dists_r1) | set(args.dists_r2))
    print(f"Dist grid: {dist_grid}", flush=True)

    print("\nBuilding gene table + polyA lookups...", flush=True)
    gene_table, polya_lookups = build_gene_lookups(
        args.genome_file, args.gtf_file, target_genes, args.lengths, args.mismatches)
    n_genes_total = sum(len(v) for v in gene_table.values())
    print(f"  Total + strand target genes mapped: {n_genes_total:,}", flush=True)

    print(f"\nClassifying paired reads with {args.cpus} workers across chromosomes...",
          flush=True)
    with pysam.AlignmentFile(args.bam_file, 'rb') as bam:
        chroms = list(bam.references)

    worker = partial(classify_chromosome,
                     bam_file=args.bam_file,
                     gene_table=gene_table,
                     polya_lookups=polya_lookups,
                     dist_grid=dist_grid,
                     cr_bam_file=args.cellranger_bam,
                     seurat_cellids=seurat_cellids,
                     min_umi_reads=args.min_umi_reads)
    if args.cellranger_bam:
        print(f"xf filter ACTIVE (Cell Ranger BAM): {args.cellranger_bam}",
              flush=True)
        if seurat_cellids is not None:
            print(f"Seurat cellID filter ACTIVE ({len(seurat_cellids):,} cells)",
                  flush=True)
        if args.min_umi_reads > 1:
            print(f"UMI read-count filter ACTIVE: keep (CB,UB) with >= "
                  f"{args.min_umi_reads} reads", flush=True)

    # Pre-populate global totals (plain dict, picklable).
    totals = {(anchor, A, M, D): {'terminal': 0, 'internal': 0,
                                  'unexplained': 0, 'total': 0}
              for anchor in ('Read1', 'Read2')
              for A in args.lengths
              for M in args.mismatches
              for D in dist_grid}

    with mp.Pool(processes=args.cpus) as pool:
        for chunk in pool.imap_unordered(worker, chroms):
            for k, v in chunk.items():
                if k not in totals:
                    continue
                for stat, n in v.items():
                    totals[k][stat] += n

    # Write tidy CSV. For Read1 keep only D values from --dists-r1, for Read2
    # only --dists-r2.
    r1_set = set(args.dists_r1)
    r2_set = set(args.dists_r2)
    with open(args.out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['anchor', 'min_length', 'max_nonA', 'min_dist',
                    'terminal', 'internal', 'unexplained', 'total_reads'])
        for anchor in ('Read1', 'Read2'):
            allowed = r1_set if anchor == 'Read1' else r2_set
            for A in args.lengths:
                for M in args.mismatches:
                    for D in dist_grid:
                        if D not in allowed:
                            continue
                        t = totals[(anchor, A, M, D)]
                        w.writerow([anchor, A, M, D,
                                    t['terminal'], t['internal'],
                                    t['unexplained'], t['total']])

    print(f"\nWrote {args.out}", flush=True)
    print(f"Total runtime: {time.perf_counter()-t0:.1f}s", flush=True)


if __name__ == '__main__':
    main()
