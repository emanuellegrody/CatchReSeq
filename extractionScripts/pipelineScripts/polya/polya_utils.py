"""
Shared utilities for the polyA analysis pipeline.

Centralizes the maximal A-rich region finder, region deduplication, and
GTF/FASTA loading so bam_polya_virus.py, fasta_polya_host.py, and
polya_comparative_stats.py operate on a single canonical implementation.

The find_polyA_in_Seq implementation is vectorized via NumPy and produces
results identical to the original pure-Python algorithm.
"""

import os
import sys
import numpy as np
import pysam
import gffutils

A_BYTE = ord('A')


def load_target_genes(gene_file):
    """Read newline-separated gene names into a set; drop blanks and # comments."""
    if not os.path.exists(gene_file):
        raise FileNotFoundError(f"Gene file not found: {gene_file}")
    genes = set()
    with open(gene_file) as f:
        for line in f:
            name = line.strip()
            if name and not name.startswith('#'):
                genes.add(name)
    return genes


def find_polyA_in_Seq(sequence, min_length, max_nonA):
    """
    Maximal A-rich regions: for each A-start i, walk forward while
    accumulated non-A count <= max_nonA, keeping the furthest A-end with
    run length >= min_length.

    Returns list of (start, end_exclusive) — same convention as the
    original Python implementation.
    """
    if not sequence:
        return []

    arr = np.frombuffer(sequence.encode('ascii', errors='replace'),
                        dtype=np.uint8)
    is_a = (arr == A_BYTE)
    n = arr.size

    a_pos = np.flatnonzero(is_a)
    if a_pos.size == 0:
        return []

    nonA_pos = np.flatnonzero(~is_a)

    # For each A-start i, find termination = position of (max_nonA+1)-th
    # non-A at or after i, or n if there isn't one. Inner loop in the
    # original increments mismatches only at non-A positions and breaks
    # the next iteration -- so termination IS that limiting non-A.
    k = np.searchsorted(nonA_pos, a_pos)            # non-A's strictly before i
    term_idx = k + max_nonA                         # index into nonA_pos
    has_term = term_idx < nonA_pos.size
    if nonA_pos.size > 0:
        clipped = np.minimum(term_idx, nonA_pos.size - 1)
        termination = np.where(has_term, nonA_pos[clipped], n)
    else:
        termination = np.full(a_pos.shape, n, dtype=np.int64)

    # Largest A position <= termination-1 (which is strictly before the
    # limiting non-A, so guaranteed to be A or undefined). prev_a built once.
    prev_a = np.where(is_a, np.arange(n), -1)
    np.maximum.accumulate(prev_a, out=prev_a)

    last_valid_end = prev_a[termination - 1]
    run_len = last_valid_end - a_pos + 1
    keep = run_len >= min_length
    if not keep.any():
        return []

    starts = a_pos[keep]
    ends_excl = last_valid_end[keep] + 1
    return list(zip(starts.tolist(), ends_excl.tolist()))


def deduplicate_regions(regions):
    """Drop nested intervals: keep regions whose end strictly exceeds the
    running max end. Sort by (start asc, end desc) so the longest extension
    wins at each start."""
    if not regions:
        return []
    regions = sorted(regions, key=lambda r: (r[0], -r[1]))
    dedup = []
    max_end = -1
    for r in regions:
        if r[1] <= max_end:
            continue
        dedup.append(r)
        max_end = r[1]
    return dedup


def open_indexed_fasta(fasta_path):
    """pysam.FastaFile with a friendly error if the .fai is missing."""
    try:
        return pysam.FastaFile(fasta_path)
    except OSError as e:
        sys.stderr.write(
            f"ERROR: Cannot open FASTA index for {fasta_path}.\n"
            f"       Run 'samtools faidx {fasta_path}' first.\n"
            f"       ({e})\n"
        )
        sys.exit(1)


def open_or_build_gtf_db(gtf_path):
    """Open the gffutils SQLite DB beside the GTF; build it once if missing."""
    db_file = f"{gtf_path}.db"
    if os.path.exists(db_file):
        return gffutils.FeatureDB(db_file)
    print(f"Building GTF database at {db_file} (one-time setup)...",
          file=sys.stderr, flush=True)
    return gffutils.create_db(gtf_path, dbfn=db_file, force=True, verbose=False)


def iter_target_gene_seqs(host_fasta_path, gtf_path, target_genes,
                          forward_only=True):
    """
    Yield (gene_name, sequence_uppercase, gene_start_0based, gene_end, strand)
    for every gene whose `gene_name` is in target_genes.

    - Repairs `chr`-prefix mismatches between GTF and FASTA when possible.
    - Skips reverse-strand genes when forward_only=True (matches Amir's
      original logic: polyA tracts on '+' strand only).
    """
    fasta = open_indexed_fasta(host_fasta_path)
    db = open_or_build_gtf_db(gtf_path)
    fasta_chroms = set(fasta.references)

    for gene in db.features_of_type('gene'):
        gene_name = gene.attributes.get('gene_name', [gene.id])[0]
        if gene_name not in target_genes:
            continue
        if forward_only and gene.strand == '-':
            continue

        chrom = gene.chrom
        if chrom not in fasta_chroms:
            alt = chrom[3:] if chrom.startswith('chr') else f"chr{chrom}"
            if alt not in fasta_chroms:
                continue
            chrom = alt

        start_0 = gene.start - 1
        end = gene.end
        try:
            seq = fasta.fetch(chrom, start_0, end).upper()
        except KeyError:
            continue

        yield gene_name, seq, start_0, end, gene.strand
