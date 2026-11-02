#!/usr/bin/env python3
"""
Per-FOV fraction of cells called infected in the HCR reference.

One point per FOV, one column per imaging group, bar at the group median.
A cell is infected when n_US + n_SS > 0: at least one nef spot colocalized
with an env spot. Figure2.Rmd applies the same rule to build the HCR
reference for Panels L and M.

The denominator is manually masked cells, not cells in the field, so the
columns compare groups under a fixed rule rather than measuring infection rate.

Usage:
    python hcr_positivity_panel.py --output_dir /path/to/plotData/Figure2/positivityRule
"""

import os
import csv
import argparse
import collections

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_TABLE_DIR = "/projects/b1042/GoyalLab/egrody/extractedData/Figure2/bigFISH"
DEFAULT_REP1 = "manual_batchthr"
DEFAULT_REP2 = "manual_replicate2_batchthr"

GROUP_COLORS = collections.OrderedDict([
    ("EGS023 rep1", "#5C6BC0"),
    ("EGS027 0820", "#40BFEF"),
    ("EGS027 0821", "#FF8A65"),
])


def read_cell_stats(path):
    rows = []
    with open(path) as handle:
        for row in csv.DictReader(handle):
            rows.append({
                "fov": row["fov"],
                "n_US": int(row["n_US"]),
                "n_SS": int(row["n_SS"]),
            })
    return rows


def group_of_fov(fov, replicate):
    # Replicate 2 splits by acquisition date. well2 of 20260821 lands in the 0821
    # column and carries a far higher spot load than any well1 field.
    if replicate == "rep1":
        return "EGS023 rep1"
    return "EGS027 " + fov.split("_")[0][4:]


def percent_positive_per_fov(rows, replicate):
    per_fov_total = collections.Counter()
    per_fov_positive = collections.Counter()
    fov_group = {}

    for row in rows:
        fov = row["fov"]
        fov_group[fov] = group_of_fov(fov, replicate)
        per_fov_total[fov] += 1
        if row["n_US"] + row["n_SS"] > 0:
            per_fov_positive[fov] += 1

    grouped = collections.defaultdict(list)
    for fov, total in per_fov_total.items():
        grouped[fov_group[fov]].append(
            (100.0 * per_fov_positive[fov] / total, total))
    return grouped


def draw_panel(grouped, output_path, seed=0):
    # Point area scales with cells masked. Small FOVs sit near 100% because only
    # positive-looking cells were masked in them, so the unweighted median runs
    # well above the cells-weighted rate; the sizes keep that visible.
    rng = np.random.RandomState(seed)
    figure, axis = plt.subplots(figsize=(6.9, 5.0))

    labels = [label for label in GROUP_COLORS if label in grouped]
    for position, label in enumerate(labels):
        values = np.asarray([v for v, _ in grouped[label]], dtype=float)
        counts = np.asarray([n for _, n in grouped[label]], dtype=float)
        jitter = rng.uniform(-0.18, 0.18, size=values.size)
        axis.scatter(position + jitter, values, s=4.0 * np.sqrt(counts),
                     alpha=0.7, color=GROUP_COLORS[label], linewidths=0)
        axis.hlines(np.median(values), position - 0.32, position + 0.32,
                    color="black", linewidth=2.5, zorder=3)
        weighted = 100.0 * np.sum(values / 100.0 * counts) / np.sum(counts)
        axis.hlines(weighted, position - 0.32, position + 0.32,
                    color="black", linewidth=1.6, linestyle=":", zorder=3)

    axis.set_xticks(range(len(labels)))
    axis.set_xticklabels(labels, rotation=12, ha="right")
    axis.set_ylabel("% cells positive  (n_US + n_SS > 0)")
    axis.set_title("fraction of cells carrying colocalized signal")
    axis.legend(handles=[
        plt.Line2D([0], [0], color="black", linewidth=2.5, label="median FOV"),
        plt.Line2D([0], [0], color="black", linewidth=1.6, linestyle=":",
                   label="cells-weighted"),
    ], loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8, frameon=False)
    axis.set_ylim(-2, 104)
    axis.yaxis.grid(True, color="#DDDDDD", linewidth=0.8)
    axis.set_axisbelow(True)
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)

    figure.tight_layout()
    figure.savefig(output_path, dpi=200)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--table_dir", default=DEFAULT_TABLE_DIR)
    parser.add_argument("--rep1", default=DEFAULT_REP1)
    parser.add_argument("--rep2", default=DEFAULT_REP2)
    parser.add_argument("--output_dir", required=True)
    arguments = parser.parse_args()

    os.makedirs(arguments.output_dir, exist_ok=True)

    grouped = collections.defaultdict(list)
    for replicate, directory in (("rep1", arguments.rep1), ("rep2", arguments.rep2)):
        table_path = os.path.join(arguments.table_dir, directory, "all_cells_stats.csv")
        if not os.path.isfile(table_path):
            parser.error("no all_cells_stats.csv at " + table_path)
        for label, values in percent_positive_per_fov(
                read_cell_stats(table_path), replicate).items():
            grouped[label].extend(values)

    panel_path = os.path.join(arguments.output_dir, "panel_percent_positive.png")
    draw_panel(grouped, panel_path)

    table_path = os.path.join(arguments.output_dir, "panel_percent_positive.tsv")
    with open(table_path, "w") as handle:
        handle.write("group\tn_fov\tn_cells\tmedian_pct\tweighted_pct\n")
        for label in GROUP_COLORS:
            if label not in grouped:
                continue
            values = np.asarray([v for v, _ in grouped[label]], dtype=float)
            counts = np.asarray([n for _, n in grouped[label]], dtype=float)
            weighted = 100.0 * np.sum(values / 100.0 * counts) / np.sum(counts)
            handle.write("{}\t{}\t{}\t{:.1f}\t{:.1f}\n".format(
                label, values.size, int(counts.sum()), np.median(values), weighted))
            print("{:<14} {:3d} FOVs {:5d} cells   median {:5.1f}%   weighted {:5.1f}%".format(
                label, values.size, int(counts.sum()), np.median(values), weighted))

    print("wrote " + panel_path)
    print("wrote " + table_path)


if __name__ == "__main__":
    main()
