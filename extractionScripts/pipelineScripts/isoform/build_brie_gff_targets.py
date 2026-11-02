#!/usr/bin/env python3
"""
Build a small BRIE2-compatible GFF3 from the EGS024 host coordinates CSV.

Reads bamsort_coordinates_host.csv (or any CSV with the same columns):
  Target, Chr, Start, End,
  Strand, Upstream_Start, Upstream_End, Downstream_Start, Downstream_End

Where:
  Start/End          = the CaRe-targeted cassette exon
  Strand             = '+' or '-'
  Upstream_*         = the constitutive exon immediately 5' of the cassette
                       (in transcript order, not genomic)
  Downstream_*       = the constitutive exon immediately 3' of the cassette
                       (in transcript order, not genomic)

Rows missing any BRIE2 column are skipped. Coordinates are 1-based, inclusive.

For each row, two isoforms are emitted:
  iso1 (inclusion): upstream + cassette + downstream
  iso2 (exclusion): upstream + downstream

Usage:
  python build_brie_gff_targets.py <coords.csv> <out.gff3>
"""

import csv
import sys


REQUIRED_BRIE_COLS = (
    "Strand",
    "Upstream_Start", "Upstream_End",
    "Downstream_Start", "Downstream_End",
)


def write_event(out, row):
    gene_id = row["Target"]
    chrom = row["Chr"]
    strand = row["Strand"]

    up_s, up_e = int(row["Upstream_Start"]), int(row["Upstream_End"])
    cs_s, cs_e = int(row["Start"]),          int(row["End"])
    dn_s, dn_e = int(row["Downstream_Start"]), int(row["Downstream_End"])

    gene_start = min(up_s, cs_s, dn_s)
    gene_end   = max(up_e, cs_e, dn_e)

    inc_id = f"{gene_id}.in"
    exc_id = f"{gene_id}.ex"

    def line(feature, start, end, attrs):
        out.write("\t".join([
            chrom, "manual", feature,
            str(start), str(end),
            ".", strand, ".", attrs
        ]) + "\n")

    line("gene", gene_start, gene_end,
         f"ID={gene_id};gene_id={gene_id};gene_name={gene_id}")

    line("mRNA", gene_start, gene_end, f"ID={inc_id};Parent={gene_id}")
    line("exon", up_s, up_e, f"ID={inc_id}.1;Parent={inc_id}")
    line("exon", cs_s, cs_e, f"ID={inc_id}.2;Parent={inc_id}")
    line("exon", dn_s, dn_e, f"ID={inc_id}.3;Parent={inc_id}")

    line("mRNA", gene_start, gene_end, f"ID={exc_id};Parent={gene_id}")
    line("exon", up_s, up_e, f"ID={exc_id}.1;Parent={exc_id}")
    line("exon", dn_s, dn_e, f"ID={exc_id}.2;Parent={exc_id}")


def main():
    if len(sys.argv) != 3:
        print("Usage: python build_brie_gff_targets.py <coords.csv> <out.gff3>")
        sys.exit(1)

    in_csv, out_gff = sys.argv[1], sys.argv[2]

    n_written = 0
    with open(in_csv) as f, open(out_gff, "w") as out:
        rows = (line for line in f
                if line.strip() and not line.lstrip().startswith("#"))
        reader = csv.DictReader(rows)
        out.write("##gff-version 3\n")
        for row in reader:
            if not all(row.get(c, "").strip() for c in REQUIRED_BRIE_COLS):
                print(f"Skipping {row.get('Target', '?')}: missing BRIE2 columns")
                continue
            write_event(out, row)
            print(f"Wrote {row['Target']}: {row['Chr']}:{row['Start']}-{row['End']} ({row['Strand']})")
            n_written += 1

    print(f"\n{n_written} event(s) written to {out_gff}")


if __name__ == "__main__":
    main()
