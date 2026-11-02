#!/usr/bin/env python3
"""
Per-FOV image check on the HCR positivity rule.

Two panels over the same FarRed (nef) plane: cells called positive by
n_total_spots > 1 on the left, by n_US + n_SS > 0 on the right, with cells the
colocalization rule drops shown in orange.

Cell outlines come from the ROIs2Label mask; its label values are the cell_id
in all_cells_stats.csv.

Detection ran on the max projection, so a dropped cell may still have a spot on
another z. Use --projection max to render what was detected on.

Usage:
    python hcr_positivity_overlay.py --fov_list fovs.tsv --output_dir /path/to/out
"""

import os
import csv
import argparse

import numpy as np
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from skimage.segmentation import find_boundaries
from skimage.morphology import binary_dilation, disk


DEFAULT_RAW_DIR = "/projects/b1042/GoyalLab/egrody/rawData"
DEFAULT_TABLE_DIR = "/projects/b1042/GoyalLab/egrody/extractedData/Figure2/bigFISH"
ARM_OF_SAMPLE = {"EGS023": "manual_batchthr", "EGS027": "manual_replicate2_batchthr"}

COLOR_POSITIVE = "#2ECC40"
COLOR_DROPPED = "#FF851B"
COLOR_NEGATIVE = "#8899A6"


def load_stats(table_dir, arm):
    stats = {}
    path = os.path.join(table_dir, arm, "all_cells_stats.csv")
    with open(path) as handle:
        for row in csv.DictReader(handle):
            stats[(row["fov"], int(row["cell_id"]))] = (
                int(row["n_US"]), int(row["n_SS"]), int(row["n_total_spots"]))
    return stats


def find_mask(raw_dir, sample, fov):
    mask_dir = os.path.join(raw_dir, sample, "masks")
    prefix = "ROIs2Label_" + fov + "_Camera"
    for name in sorted(os.listdir(mask_dir)):
        if name.startswith(prefix) and name.endswith(".tif"):
            return os.path.join(mask_dir, name)
    return None


def find_channel_planes(raw_dir, sample, fov, channel):
    planes = {}
    sample_dir = os.path.join(raw_dir, sample)
    for date in sorted(os.listdir(sample_dir)):
        for leaf in (os.path.join(sample_dir, date, "TIFF"),
                     os.path.join(sample_dir, date, "TIFF", "Export")):
            if not os.path.isdir(leaf):
                continue
            for name in os.listdir(leaf):
                marker = fov + "_Camera - " + channel + "_Z"
                if name.startswith(marker) and name.endswith(".tif"):
                    planes[name[len(marker):-4]] = os.path.join(leaf, name)
    return planes


def stretch(image, low=50.0, high=99.99):
    # Cell bodies sit only ~1.2x above background while puncta reach several
    # thousand counts above it, so the upper cut has to sit far out in the tail
    # or every spot clips to white and cells render as flat discs.
    lo, hi = np.percentile(image, [low, high])
    if hi <= lo:
        hi = lo + 1.0
    return np.clip((image.astype(np.float32) - lo) / (hi - lo), 0, 1)


def draw_panel(axis, background, mask, boundary, label_color, title):
    axis.imshow(background, cmap="gray", vmin=0, vmax=1, interpolation="nearest")

    colors = sorted(set(label_color.values()))
    color_index = {color: i + 1 for i, color in enumerate(colors)}
    lut = np.zeros(int(mask.max()) + 1, dtype=np.uint8)
    for label_value, color in label_color.items():
        lut[label_value] = color_index[color]

    boundary_category = np.zeros(mask.shape, dtype=np.uint8)
    boundary_category[boundary] = lut[mask[boundary]]

    overlay = np.zeros(background.shape + (4,), dtype=np.float32)
    for color, index in color_index.items():
        selection = binary_dilation(boundary_category == index, disk(1))
        overlay[selection] = matplotlib.colors.to_rgba(color)

    axis.imshow(overlay, interpolation="nearest")
    axis.set_title(title, fontsize=11)
    axis.set_xticks([])
    axis.set_yticks([])


def render_fov(sample, fov, category, arguments, stats):
    mask_path = find_mask(arguments.raw_dir, sample, fov)
    if mask_path is None:
        print("SKIP (no mask): " + fov)
        return None
    mask_z = os.path.splitext(os.path.basename(mask_path))[0].split("_")[-1].lstrip("Z")

    planes = find_channel_planes(arguments.raw_dir, sample, fov, arguments.channel)
    if not planes:
        print("SKIP (no {} planes): {}".format(arguments.channel, fov))
        return None

    mask = tifffile.imread(mask_path)
    if arguments.projection == "max":
        stack = np.stack([tifffile.imread(planes[z]) for z in sorted(planes)])
        image = stack.max(axis=0)
        plane_label = "MIP of {} z".format(len(planes))
    else:
        if mask_z not in planes:
            print("SKIP (mask z {} absent): {}".format(mask_z, fov))
            return None
        image = tifffile.imread(planes[mask_z])
        plane_label = "Z" + mask_z

    if mask.shape != image.shape:
        print("SKIP (shape {} vs {}): {}".format(mask.shape, image.shape, fov))
        return None

    background = stretch(image)
    boundary = find_boundaries(mask, mode="inner")

    old_color, new_color = {}, {}
    n_old = n_new = n_cells = 0
    for cell_id in [int(v) for v in np.unique(mask) if v > 0]:
        record = stats.get((fov, cell_id))
        if record is None:
            continue
        n_US, n_SS, n_total = record
        n_cells += 1
        is_old = n_total > 1
        is_new = (n_US + n_SS) > 0
        n_old += int(is_old)
        n_new += int(is_new)
        old_color[cell_id] = COLOR_POSITIVE if is_old else COLOR_NEGATIVE
        new_color[cell_id] = (COLOR_POSITIVE if is_new
                              else COLOR_DROPPED if is_old else COLOR_NEGATIVE)

    if n_cells == 0:
        print("SKIP (no matching cells): " + fov)
        return None

    pct_old = 100.0 * n_old / n_cells
    pct_new = 100.0 * n_new / n_cells

    figure, axes = plt.subplots(1, 2, figsize=(15, 7.8))
    draw_panel(axes[0], background, mask, boundary, old_color,
               "old rule  n_total_spots > 1\n{:.1f}% positive  ({}/{} cells)".format(
                   pct_old, n_old, n_cells))
    draw_panel(axes[1], background, mask, boundary, new_color,
               "new rule  n_US + n_SS > 0\n{:.1f}% positive  ({}/{} cells)".format(
                   pct_new, n_new, n_cells))
    axes[1].legend(handles=[
        Patch(facecolor=COLOR_POSITIVE, label="positive under the rule shown"),
        Patch(facecolor=COLOR_DROPPED, label="dropped by new rule"),
        Patch(facecolor=COLOR_NEGATIVE, label="negative under the rule shown"),
    ], loc="upper right", fontsize=8, framealpha=0.85)

    figure.suptitle("{}   [{}]   {} {}   -  {} dropped by new rule".format(
        fov, category, arguments.channel, plane_label, n_old - n_new), fontsize=12)
    figure.tight_layout(rect=[0, 0, 1, 0.94])

    name = (category + "_" if category else "") + fov + ".png"
    output_path = os.path.join(arguments.output_dir, name)
    figure.savefig(output_path, dpi=160)
    plt.close(figure)

    print("wrote {}  ({} cells, old {:.1f}% -> new {:.1f}%)".format(
        name, n_cells, pct_old, pct_new))
    return (category, fov, n_cells, pct_old, pct_new, n_old - n_new)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--fov_list", required=True,
        help="TSV of sample<TAB>fov[<TAB>category], one FOV per line",
    )
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--raw_dir", default=DEFAULT_RAW_DIR)
    parser.add_argument("--table_dir", default=DEFAULT_TABLE_DIR)
    parser.add_argument("--channel", default="FarRed")
    parser.add_argument("--projection", choices=("maskz", "max"), default="maskz")
    arguments = parser.parse_args()

    os.makedirs(arguments.output_dir, exist_ok=True)

    selected = []
    with open(arguments.fov_list) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                parser.error("bad line in --fov_list: " + line.rstrip("\n"))
            selected.append((fields[0], fields[1],
                             fields[2] if len(fields) > 2 else ""))

    unknown = sorted({s for s, _, _ in selected} - set(ARM_OF_SAMPLE))
    if unknown:
        parser.error("no table arm known for sample(s): " + ", ".join(unknown))

    stats_cache = {}
    summary_rows = []
    for sample, fov, category in selected:
        arm = ARM_OF_SAMPLE[sample]
        if arm not in stats_cache:
            stats_cache[arm] = load_stats(arguments.table_dir, arm)
        row = render_fov(sample, fov, category, arguments, stats_cache[arm])
        if row is not None:
            summary_rows.append(row)

    summary_path = os.path.join(arguments.output_dir, "overlay_summary.tsv")
    with open(summary_path, "w") as handle:
        handle.write("category\tfov\tn_cells\tpct_old\tpct_new\tn_dropped\n")
        for row in summary_rows:
            handle.write("{}\t{}\t{}\t{:.1f}\t{:.1f}\t{}\n".format(*row))
    print("wrote " + summary_path)


if __name__ == "__main__":
    main()
