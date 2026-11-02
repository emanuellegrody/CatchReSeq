#!/usr/bin/env python3
"""Correct the PID8027 library labels in our HIVseq outputs.

GEO GSE266455 attaches the SRA runs of six PID8027 GSMs to the wrong sample: the
runs filed under replicate N are the library published as replicate N+3 (and
vice versa), for N = 1, 2, 3. Our pipeline followed that linkage faithfully, so
six of our eight PID8027 libraries carry the wrong week and four the wrong
Primer/NOTHING chemistry. Established by cell-barcode overlap against the
published matrices (correct pairing J ~ 0.97, every other pairing J ~ 0.005) and
corroborated by HIV UMI totals; see barcode_overlap_matrix.py.

The corrected names are a permutation of the names already on disk - three
2-cycles - so this is a swap, not a rename, and it runs in two phases through a
temporary suffix. Being an involution, running it twice returns the original
state.

Prints what it would do and exits; pass --apply to perform the swap.
"""

import argparse
import csv
from pathlib import Path
import os
import sys

# Each pair swaps: the data currently under the left name is the library named
# on the right, and vice versa.
SWAPS = [
    ("PID8027week0replicate1NOTHING",  "PID8027week45replicate4NOTHING"),
    ("PID8027week0replicate2NOTHING",  "PID8027week0replicate5Primer"),
    ("PID8027week0replicate6Primer",   "PID8027week45replicate3NOTHING"),
]

RENAME_MAP = {}
for left, right in SWAPS:
    RENAME_MAP[left] = right
    RENAME_MAP[right] = left

TMP_SUFFIX = ".__pid8027swap__"

# Directories whose entry names embed a sample name. publicData/HIVseq/matrices
# is deliberately absent: those are the published files and already correct.
TARGET_DIRS = [
    "/projects/b1042/GoyalLab/egrody/publicData/HIVseq/counts",
    "/projects/b1042/GoyalLab/egrody/extractedData/HIVseq/counts",
    "/projects/b1042/GoyalLab/egrody/extractedData/HIVseq/bamsort",
    "/projects/b1042/GoyalLab/egrody/extractedData/HIVseq/bamsort/splice",
    "/projects/b1042/GoyalLab/egrody/extractedData/HIVseq/hiv",
    "/projects/b1042/GoyalLab/egrody/extractedData/HIVseq/cladeB_comparison",
    "/projects/b1042/GoyalLab/egrody/extractedData/HIVseq/Seurat/RDS",
    "/projects/b1042/GoyalLab/egrody/extractedData/HIVseq/Seurat/RDS_publishedQC",
]

# This script lives in extractionScripts/pipelineScripts/, so the repo root is
# two levels up. Resolving it this way keeps the manifest path below working
# wherever the repo is cloned.
REPO_ROOT = Path(__file__).resolve().parents[3]

# CSVs holding a sample name in a column rather than in the filename.
TARGET_CSVS = [
    ("/projects/b1042/GoyalLab/egrody/extractedData/HIVseq/bamsort/"
     "all_samples_alignment_summary.csv", "sample"),
    ("/projects/b1042/GoyalLab/egrody/extractedData/HIVseq/cladeB_comparison/"
     "comparison_summary.csv", "sample"),
    ("/projects/b1042/GoyalLab/egrody/publicData/HIVseq/scripts/"
     "HIVseq_manifest.csv", "SAMPLE"),
    (str(REPO_ROOT / "extractionScripts" / "submissionScripts" /
         "PublicData" / "HIVseq" / "samples_manifest.csv"), "SAMPLE"),
]


def corrected(name):
    """Apply the swap to any entry name embedding one of the six sample names."""
    for old, new in RENAME_MAP.items():
        if old in name:
            return name.replace(old, new)
    return None


def plan_dir(directory):
    if not os.path.isdir(directory):
        return []
    moves = []
    for entry in sorted(os.listdir(directory)):
        new = corrected(entry)
        if new is not None and new != entry:
            moves.append((entry, new))
    return moves


def nested_moves(root, old_sample, new_sample):
    """Entries below a library directory that also embed its sample name - the
    cellranger .mri.tgz and the per-sample STAR outputs under hiv/. Bottom-up so
    a directory is renamed only after its contents.

    Every entry here carries this one library's name, so the corrected names
    cannot collide with one another and no temporary phase is needed.
    """
    moves = []
    for parent, dirs, files in os.walk(root, topdown=False):
        for entry in sorted(dirs) + sorted(files):
            if old_sample in entry:
                moves.append((os.path.join(parent, entry),
                              os.path.join(parent,
                                           entry.replace(old_sample,
                                                         new_sample))))
    return moves


def sample_of(entry):
    """Which of the six swapped sample names this entry name embeds."""
    for old in RENAME_MAP:
        if old in entry:
            return old
    return None


def swap_dir(directory, moves, apply_changes):
    """Two-phase swap through a temporary name, so the two halves of a 2-cycle
    never collide. Nested entries are renamed while the library sits at its
    temporary name, which keeps each library's rename self-contained and makes
    the dry run report exactly what would run."""
    nested_count = 0
    staged = []
    for old, new in moves:
        src = os.path.join(directory, old)
        tmp = os.path.join(directory, old + TMP_SUFFIX)
        if apply_changes:
            os.rename(src, tmp)
        staged.append((tmp if apply_changes else src, old, new))

    for path, old, new in staged:
        old_sample, new_sample = sample_of(old), sample_of(new)
        if os.path.isdir(path) and old_sample:
            inner = nested_moves(path, old_sample, new_sample)
            nested_count += len(inner)
            if apply_changes:
                for src, dst in inner:
                    if os.path.exists(dst):
                        raise RuntimeError(f"refusing to overwrite {dst}")
                    os.rename(src, dst)

    for path, old, new in staged:
        dst = os.path.join(directory, new)
        if apply_changes:
            if os.path.exists(dst):
                raise RuntimeError(f"refusing to overwrite {dst}")
            os.rename(path, dst)
    return nested_count


def plan_csv(path, column):
    if not os.path.exists(path):
        return None
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or column not in rows[0]:
        return None
    changed = sum(1 for r in rows if corrected(r[column]) not in (None, r[column]))
    return rows, changed


def rewrite_csv(path, column, rows, apply_changes):
    for row in rows:
        new = corrected(row[column])
        if new is not None:
            row[column] = new
    if apply_changes:
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()),
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true",
                        help="perform the swap (default is a dry run)")
    args = parser.parse_args()

    print("PID8027 label correction -",
          "APPLYING" if args.apply else "DRY RUN (pass --apply to perform)")
    print("\nswap pairs:")
    for left, right in SWAPS:
        print(f"  {left:<32} <-> {right}")

    total, nested_total = 0, 0
    for directory in TARGET_DIRS:
        if not os.path.isdir(directory):
            print(f"\n[missing] {directory}")
            continue
        moves = plan_dir(directory)
        nested = swap_dir(directory, moves, args.apply)
        print(f"\n{directory}  ({len(moves)} entries, "
              f"{nested} nested)")
        for old, new in moves[:4]:
            print(f"    {old}\n      -> {new}")
        if len(moves) > 4:
            print(f"    ... and {len(moves) - 4} more")
        total += len(moves)
        nested_total += nested

    print(f"\n--- {total} top-level + {nested_total} nested entries ---\n")

    for path, column in TARGET_CSVS:
        result = plan_csv(path, column)
        if result is None:
            print(f"[skip] {path} (missing, or no '{column}' column)")
            continue
        rows, changed = result
        print(f"[csv ] {path}: {changed} rows in '{column}'")
        rewrite_csv(path, column, rows, args.apply)

    if not args.apply:
        print("\nNothing was changed. Re-run with --apply.")
    else:
        print("\nDone. This script is its own inverse: running it again reverts.")
        print("Still stale and NOT handled here: orig.ident inside the "
              "per-sample Seurat RDS objects, and the integrated objects, "
              "whose condition membership was built from the wrong labels.")


if __name__ == "__main__":
    main()
