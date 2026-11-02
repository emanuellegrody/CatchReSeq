#!/usr/bin/env python3
"""Match our cellranger runs to the published GSE266455 matrices by cell-barcode
overlap rather than by name, to check the name-based mapping.

Two libraries of the same 10X emulsion share nearly all called barcodes; two
different emulsions share only the few hundred expected by chance. Printed as a
Jaccard matrix per participant, with the best published match for each of ours.
"""

import gzip
import os
import sys
from collections import defaultdict

COUNTS_DIR = "/projects/b1042/GoyalLab/egrody/publicData/HIVseq/counts"
MATRIX_DIR = "/projects/b1042/GoyalLab/egrody/publicData/HIVseq/matrices"


def read_barcodes(path):
    with gzip.open(path, "rt") as handle:
        return {line.strip() for line in handle if line.strip()}


def main():
    participants = sys.argv[1:] or ["PID1052", "PID145", "PID8026", "PID8027"]

    ours, published = defaultdict(list), defaultdict(list)
    for d in sorted(os.listdir(COUNTS_DIR)):
        if not d.startswith("GRCh38_HXB2_"):
            continue
        sample = d[len("GRCh38_HXB2_"):]
        pid = sample.split("week")[0]
        path = os.path.join(COUNTS_DIR, d, "outs",
                            "filtered_feature_bc_matrix", "barcodes.tsv.gz")
        if os.path.exists(path):
            ours[pid].append((sample, path))
    for d in sorted(os.listdir(MATRIX_DIR)):
        path = os.path.join(MATRIX_DIR, d, "barcodes.tsv.gz")
        if os.path.exists(path):
            published[d.split("W")[0]].append((d, path))

    for pid in participants:
        mine, theirs = ours.get(pid, []), published.get(pid, [])
        if not mine or not theirs:
            continue
        mine_bc = {s: read_barcodes(p) for s, p in mine}
        theirs_bc = {s: read_barcodes(p) for s, p in theirs}

        print(f"\n===== {pid}: Jaccard of called-barcode sets =====")
        width = max(len(s) for s in mine_bc) + 2
        print("".ljust(width) + "".join(t.split("W", 1)[1].ljust(10)
                                        for t in theirs_bc))
        for s, a in mine_bc.items():
            cells = []
            for t, b in theirs_bc.items():
                j = len(a & b) / len(a | b)
                cells.append((f"{j:.3f}" if j >= 0.01 else ".").ljust(10))
            print(s.ljust(width) + "".join(cells))

        print(f"--- best published match per sample ({pid}) ---")
        for s, a in mine_bc.items():
            scored = sorted(((len(a & b) / len(a | b), t)
                             for t, b in theirs_bc.items()), reverse=True)
            best_j, best = scored[0]
            runner_j, runner = scored[1] if len(scored) > 1 else (0.0, "-")
            flag = "" if best_j > 0.5 else "   <-- NO CONFIDENT MATCH"
            print(f"  {s:<34} -> {best:<22} J={best_j:.3f} "
                  f"(next {runner} {runner_j:.3f}){flag}")


if __name__ == "__main__":
    main()
