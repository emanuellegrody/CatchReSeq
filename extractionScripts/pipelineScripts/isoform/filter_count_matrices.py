#!/usr/bin/env python3
"""
Filter barcodes for BRIE2 based on cellIDs expressing target genes.

For each gene column in the CSV, finds matching sample directories in counts_dir,
filters the barcode list to only include cellIDs with nonzero counts for that gene,
and writes filtered barcode files to the output directory.

Usage:
    python filter_barcodes.py <csv_file> <condition> <counts_dir> <output_dir>

    csv_file:    Path to CSV with 'cellID' column and gene columns (e.g., Kras, Dgkd)
    condition:   'WT' or 'KO'
    counts_dir:  Path to Cell Ranger counts directory containing sample folders
    output_dir:  Path to BRIE2 output directory
"""

import sys
import os
import gzip
import pandas as pd


def find_matching_dirs(gene, condition, counts_dir):
    """Find sample directories matching a gene name and condition (WT/KO)."""
    matches = []
    for dirname in os.listdir(counts_dir):
        dirpath = os.path.join(counts_dir, dirname)
        if not os.path.isdir(dirpath):
            continue
        # Match condition
        if condition == "WT" and "_WT_" not in dirname:
            continue
        if condition == "KO" and "_KO_" not in dirname:
            continue
        # Match gene: directory name ends with the gene name or gene name
        # is a prefix of the suffix (e.g., "Slk" matches "Slkv3", "Slkv4")
        suffix = dirname.split("_")[-1]
        if suffix.startswith(gene):
            matches.append(dirname)
    return sorted(matches)


def read_barcodes(barcode_path):
    """Read barcodes from a gzipped or plain text file."""
    if barcode_path.endswith(".gz"):
        with gzip.open(barcode_path, "rt") as f:
            return [line.strip() for line in f]
    else:
        with open(barcode_path, "r") as f:
            return [line.strip() for line in f]


def main():
    if len(sys.argv) != 5:
        print("Usage: python filter_barcodes.py <csv_file> <condition> <counts_dir> <output_dir>")
        sys.exit(1)

    csv_file = sys.argv[1]
    condition = sys.argv[2]
    counts_dir = sys.argv[3]
    output_dir = sys.argv[4]

    df = pd.read_csv(csv_file)

    if "cellID" not in df.columns:
        print("Error: CSV must contain a 'cellID' column.")
        sys.exit(1)

    gene_columns = [c for c in df.columns if c != "cellID"]
    print(f"Found {len(gene_columns)} gene columns: {gene_columns}")
    print(f"Condition: {condition}")
    print(f"Total rows in CSV: {len(df)}")

    for gene in gene_columns:
        # Get cellIDs with nonzero counts for this gene
        mask = pd.to_numeric(df[gene], errors="coerce").fillna(0) > 0
        cell_ids = set(df.loc[mask, "cellID"].astype(str).str.strip())
        print(f"\n--- {gene} ---")
        print(f"  cellIDs with nonzero counts: {len(cell_ids)}")

        # Find matching sample directories
        matching_dirs = find_matching_dirs(gene, condition, counts_dir)
        if not matching_dirs:
            print(f"  WARNING: No matching directories found for gene={gene}, condition={condition}")
            continue

        print(f"  Matching directories: {matching_dirs}")

        for sample_dir in matching_dirs:
            barcode_path = os.path.join(
                counts_dir, sample_dir, "outs",
                "filtered_feature_bc_matrix", "barcodes.tsv.gz"
            )
            if not os.path.exists(barcode_path):
                print(f"  WARNING: Barcode file not found: {barcode_path}")
                continue

            all_barcodes = read_barcodes(barcode_path)

            # Filter: keep barcodes that appear in the cellID list
            # Handle possible suffix mismatches (e.g., "-1" suffix in 10x barcodes)
            filtered = [bc for bc in all_barcodes if bc in cell_ids or bc.split("-")[0] in cell_ids]

            # Also try matching without the -1 suffix from the cellID side
            if len(filtered) == 0:
                cell_ids_stripped = {c.split("-")[0] for c in cell_ids}
                filtered = [bc for bc in all_barcodes if bc.split("-")[0] in cell_ids_stripped]

            short_name = sample_dir.replace("GRCm39_EGL002_", "")
            out_subdir = os.path.join(output_dir, short_name)
            os.makedirs(out_subdir, exist_ok=True)

            out_path = os.path.join(out_subdir, "filtered_barcodes.tsv")
            with open(out_path, "w") as f:
                for bc in filtered:
                    f.write(bc + "\n")

            print(f"  {sample_dir} -> {short_name}/filtered_barcodes.tsv: "
                  f"{len(filtered)}/{len(all_barcodes)} barcodes kept")


if __name__ == "__main__":
    main()