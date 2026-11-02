#!/usr/bin/env python3
"""Correct the PID1052 timepoint labels in our HIVseq outputs.

Replicates 2, 4 and 6 are PID1052's ART-suppressed timepoint, and its week is
70. Four sources disagree, three of them wrong:

  * GEO GSE266455 titles all six PID1052 GSMs "week0", and the Series
    overall_design says suppression was at "week 24 or week 45" and never
    mentions week 70. Together these make PID1052 look viremia-only. It is not;
  * the published matrices name replicates 2, 4 and 6 as W24
    (PID1052W24w2NOTHING, PID1052W24w4HIV, PID1052W24w6HIV), which is right that
    they are the suppressed timepoint but wrong in the week: W24 is PID8026's;
  * the paper Methods are authoritative (Frouard et al., Nat Commun 2026,
    doi 10.1038/s41467-026-68797-3): "Samples were collected prior to ART
    initiation (Week 0) and following ART suppression (Week 24: PID8026;
    Week 45: PID8027; and Week 70: PID1052)".

That replicates 2, 4 and 6 are the suppressed timepoint at all rests on three
further grounds:

  * each of our runs pairs with its published matrix at 96-98% barcode overlap
    (logs/comparison_bcmatch_thr1.csv), and the pairing is replicate-for-
    replicate, so only the week token is in question;
  * split this way our HIV+ counts tracked the paper's at the time of the
    correction: the paper reports 8 HIV RNA+ cells viremic (replicates 1, 3, 5)
    and 14 on ART (replicates 2, 4, 6), and ours came to 9 and 13. Split by the
    GEO titles, all six are "week0" and neither number can be recovered. The
    pipeline has changed since, so do not read 9 and 13 as current;
  * the paper describes three viremic PID1052 libraries - 4,263 CD4 T cells for
    NOTHING and 8,772 for Primer, summing to the 13,035 it reports viremic - not
    six.

Not needed for a run started from runAll_publicData.sh: samples_manifest.csv
already carries week70 in SAMPLE, with GEO's name kept in ORIGINAL_GEO_SAMPLE,
and adtgate_hivseq.R takes the week from our own names. This script is only for
repairing output already built under GEO's week0 names.

Unlike the PID8027 correction this is a directed rename, not a permutation: the
three target names do not exist on disk, so there is no collision and no
temporary phase. It is therefore NOT its own inverse; pass --revert to undo.

Prints what it would do and exits; pass --apply to perform the rename.
"""

import argparse
import csv
import re
from pathlib import Path
import os

# The data currently under the left name is the library named on the right.
RENAMES = [
    ("PID1052week0replicate2NOTHING", "PID1052week70replicate2NOTHING"),
    ("PID1052week0replicate4Primer",  "PID1052week70replicate4Primer"),
    ("PID1052week0replicate6Primer",  "PID1052week70replicate6Primer"),
]

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

# CSVs holding a sample name in one or more columns rather than in the filename,
# as (path, name columns, week column or None). Every name column is rewritten,
# not just the first: the ADT tables carry the name twice, once as `sample` and
# once as the `cellID` prefix, and patching only one leaves the two keys
# pointing at different libraries. A week column holds the bare number, so it is
# derived from the corrected name rather than substituted into. The manifests'
# ORIGINAL_GEO_SAMPLE is deliberately absent: it records GEO's own name.
TARGET_CSVS = [
    ("/projects/b1042/GoyalLab/egrody/extractedData/HIVseq/bamsort/"
     "all_samples_alignment_summary.csv", ["sample"], None),
    ("/projects/b1042/GoyalLab/egrody/extractedData/HIVseq/cladeB_comparison/"
     "comparison_summary.csv", ["sample"], None),
    ("/projects/b1042/GoyalLab/egrody/publicData/HIVseq/scripts/"
     "HIVseq_manifest.csv", ["SAMPLE"], None),
    (str(REPO_ROOT / "extractionScripts" / "submissionScripts" /
         "PublicData" / "HIVseq" / "samples_manifest.csv"), ["SAMPLE"], None),
    ("/projects/b1042/GoyalLab/egrody/extractedData/HIVseq/ADT/"
     "CD4T_cellIDs.csv", ["sample", "cellID"], "week"),
    ("/projects/b1042/GoyalLab/egrody/extractedData/HIVseq/ADT/"
     "CD4T_gate_summary.csv", ["sample"], "week"),
]

WEEK_RE = re.compile(r"week(\d+)replicate")


def build_map(revert):
    return {new: old for old, new in RENAMES} if revert else dict(RENAMES)


def corrected(name, mapping):
    """Apply the rename to any entry name embedding one of the three names."""
    for old, new in mapping.items():
        if old in name:
            return name.replace(old, new)
    return None


def plan_dir(directory, mapping):
    if not os.path.isdir(directory):
        return []
    return [(e, corrected(e, mapping)) for e in sorted(os.listdir(directory))
            if corrected(e, mapping) not in (None, e)]


def nested_moves(root, old_sample, new_sample):
    """Entries below a library directory that also embed its sample name - the
    cellranger .mri.tgz and the per-sample STAR outputs under hiv/. Bottom-up so
    a directory is renamed only after its contents."""
    moves = []
    for parent, dirs, files in os.walk(root, topdown=False):
        for entry in sorted(dirs) + sorted(files):
            if old_sample in entry:
                moves.append((os.path.join(parent, entry),
                              os.path.join(parent, entry.replace(old_sample,
                                                                 new_sample))))
    return moves


def sample_of(entry, mapping):
    for old in mapping:
        if old in entry:
            return old
    return None


def rename_dir(directory, moves, mapping, apply_changes):
    """Rename nested entries first, then the library directory itself. No
    temporary phase: old and new names are disjoint sets."""
    nested_count = 0
    for old, new in moves:
        src = os.path.join(directory, old)
        dst = os.path.join(directory, new)
        old_sample = sample_of(old, mapping)
        new_sample = corrected(old_sample, mapping) if old_sample else None
        if os.path.isdir(src) and old_sample:
            inner = nested_moves(src, old_sample, new_sample)
            nested_count += len(inner)
            if apply_changes:
                for isrc, idst in inner:
                    if os.path.exists(idst):
                        raise RuntimeError(f"refusing to overwrite {idst}")
                    os.rename(isrc, idst)
        if apply_changes:
            if os.path.exists(dst):
                raise RuntimeError(f"refusing to overwrite {dst}")
            os.rename(src, dst)
    return nested_count


def plan_csv(path, columns, mapping):
    """Rows of `path`, and how many of them any of `columns` would change."""
    if not os.path.exists(path):
        return None
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    present = [c for c in columns if rows and c in rows[0]]
    if not present:
        return None
    changed = sum(1 for r in rows
                  if any(corrected(r[c], mapping) not in (None, r[c])
                         for c in present))
    return rows, present, changed


def rewrite_csv(path, columns, week_column, rows, mapping, apply_changes):
    for row in rows:
        touched = False
        for column in columns:
            new = corrected(row[column], mapping)
            if new is not None:
                row[column] = new
                touched = True
        # The week is read back off whichever name was just corrected, so it
        # cannot drift from it, and a revert carries it back to 0.
        if touched and week_column and week_column in row:
            match = WEEK_RE.search(row[columns[0]])
            if match:
                row[week_column] = match.group(1)
    if apply_changes:
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()),
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true",
                        help="perform the rename (default is a dry run)")
    parser.add_argument("--revert", action="store_true",
                        help="rename week70 back to week0")
    args = parser.parse_args()
    mapping = build_map(args.revert)

    print("PID1052 label correction -",
          "APPLYING" if args.apply else "DRY RUN (pass --apply to perform)",
          "(REVERT)" if args.revert else "")
    print("\nrenames:")
    for old, new in mapping.items():
        print(f"  {old:<32} -> {new}")

    total, nested_total = 0, 0
    for directory in TARGET_DIRS:
        if not os.path.isdir(directory):
            print(f"\n[missing] {directory}")
            continue
        moves = plan_dir(directory, mapping)
        nested = rename_dir(directory, moves, mapping, args.apply)
        print(f"\n{directory}  ({len(moves)} entries, {nested} nested)")
        for old, new in moves:
            print(f"    {old}\n      -> {new}")
        total += len(moves)
        nested_total += nested

    print(f"\n--- {total} top-level + {nested_total} nested entries ---\n")

    for path, columns, week_column in TARGET_CSVS:
        result = plan_csv(path, columns, mapping)
        if result is None:
            print(f"[skip] {path} (missing, or none of "
                  f"{', '.join(columns)} present)")
            continue
        rows, present, changed = result
        week_note = f" (+ '{week_column}')" if week_column else ""
        print(f"[csv ] {path}: {changed} rows in "
              f"{', '.join(repr(c) for c in present)}{week_note}")
        rewrite_csv(path, present, week_column, rows, mapping, args.apply)

    if not args.apply:
        print("\nNothing was changed. Re-run with --apply.")
    else:
        print("\nDone. Pass --revert --apply to undo.")
        print("Still stale and NOT handled here: the three per-sample Seurat "
              "objects, whose contents were built under the old name, and any "
              "integrated object that included them. Delete and rebuild "
              "through the Figure3 Extraction chunk.")


if __name__ == "__main__":
    main()
