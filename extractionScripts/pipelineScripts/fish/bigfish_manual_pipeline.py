#!/usr/bin/env python3
"""
big-FISH Multi-FOV Processing Pipeline with Manual Segmentation

Calculates spliceforms per cell using manually segmented cell masks from Fiji.
Processes multiple fields of view (FOVs) with automatic z-slice discovery.
Uses pre-made cell segmentation masks (label images) - no nuclei used.

Usage:
    python bigfish_manual_pipeline.py --input_dir /path/to/input --output_dir /path/to/output \
    --mask_dir /path/to/masks --channels Red Gold70 FarRed --fov_pattern=r'(\d{8}_[A-Z0-9]+_[a-z]+_slide\d+_\d+(?:_XY\d+)?)'
"""

import os
import re
import sys
import csv
import argparse
import traceback
import numpy as np
from pathlib import Path

# scipy
from scipy import ndimage

# big-FISH
import bigfish
import bigfish.stack as stack
import bigfish.multistack as multistack
import bigfish.detection as detection

# matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# tifffile for reading label images
import tifffile


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def discover_z_slices(path_input, fov, channels):
    """
    Auto-discover all z-slices available for a given FOV.

    Parameters
    ----------
    path_input : str
        Directory containing TIF files
    fov : str
        FOV identifier
    channels : list
        List of channel names to validate

    Returns
    -------
    list or None
        Sorted list of z-slice identifiers, or None if channels missing
    """
    z_slices = set()
    z_pattern = re.compile(r'(Z\d+)')

    for filename in os.listdir(path_input):
        if fov in filename and filename.endswith('.tif'):
            z_match = z_pattern.search(filename)
            if z_match:
                z_slices.add(z_match.group(1))

    if not z_slices:
        return None

    return sorted(list(z_slices))


def validate_fov_channels(path_input, fov, channels, opt_prefix="Camera - "):
    """
    Validate that all required channels exist for a given FOV.

    Parameters
    ----------
    path_input : str
        Directory containing TIF files
    fov : str
        FOV identifier
    channels : list
        Required channel names
    opt_prefix : str
        Optical configuration prefix in filename

    Returns
    -------
    tuple
        (bool, list) - (all_present, missing_channels)
    """
    files_in_dir = set(os.listdir(path_input))
    missing_channels = []

    for channel in channels:
        # Check if at least one file exists for this FOV+channel combination
        pattern_str = "{}.*{}{}.*\\.tif$".format(
            re.escape(fov),
            re.escape(opt_prefix),
            re.escape(channel)
        )
        channel_pattern = re.compile(pattern_str)
        channel_found = any(channel_pattern.search(f) for f in files_in_dir)

        if not channel_found:
            missing_channels.append(channel)

    return len(missing_channels) == 0, missing_channels


def discover_fovs(path_input, fov_pattern=r'(\d{8}_[A-Z0-9]+_[a-z]+_slide\d+_\d+(?:_XY\d+)?)'):
    """
    Auto-discover all unique FOVs in the directory.

    Parameters
    ----------
    path_input : str
        Directory containing TIF files
    fov_pattern : str
        Regex pattern matching FOV identifier format

    Returns
    -------
    list
        Sorted list of unique FOV identifiers
    """
    fovs = set()
    pattern = re.compile(fov_pattern)

    for filename in os.listdir(path_input):
        if filename.endswith('.tif'):
            match = pattern.search(filename)
            if match:
                fovs.add(match.group(1))

    if not fovs:
        raise ValueError("No FOVs found in directory: {}".format(path_input))

    return sorted(list(fovs))


def find_mask_for_fov(mask_dir, fov, fov_pattern=r'(\d{8}_[A-Z0-9]+_[a-z]+_slide\d+_\d+(?:_XY\d+)?)'):
    """
    Find the corresponding mask file for a given FOV.

    Mask files are expected to be named like:
    ROIs2Label_20251119_EGS023_bingo_slide6_001_XY1_Camera - FarRed_Z09.tif

    Parameters
    ----------
    mask_dir : str
        Directory containing mask TIF files
    fov : str
        FOV identifier to match
    fov_pattern : str
        Regex pattern for FOV extraction

    Returns
    -------
    str or None
        Path to mask file, or None if not found
    """
    pattern = re.compile(fov_pattern)

    for filename in os.listdir(mask_dir):
        if filename.endswith('.tif') and filename.startswith('ROIs2Label_'):
            match = pattern.search(filename)
            if match and match.group(1) == fov:
                return os.path.join(mask_dir, filename)

    return None


def load_cell_mask(mask_path):
    """
    Load a cell segmentation mask (label image) from a TIF file.

    Parameters
    ----------
    mask_path : str
        Path to the mask TIF file

    Returns
    -------
    np.ndarray
        Label image where each cell has a unique integer value, background is 0
    """
    mask = tifffile.imread(mask_path)

    # Ensure it's 2D
    if mask.ndim > 2:
        # Take the first slice if it's a stack
        mask = mask[0] if mask.ndim == 3 else mask[0, 0]
        print("  Warning: Mask was >2D, using first slice")

    # Ensure integer type
    if not np.issubdtype(mask.dtype, np.integer):
        mask = mask.astype(np.int32)

    return mask


def spots_to_structured(spots):
    """
    Convert spots array to structured array for exact matching.

    Parameters
    ----------
    spots : np.ndarray
        Spots array with shape (n_spots, n_coords)

    Returns
    -------
    np.recarray
        Structured array for efficient set operations
    """
    if spots.shape[1] == 2:
        dtype = [('y', spots.dtype), ('x', spots.dtype)]
    elif spots.shape[1] == 3:
        dtype = [('y', spots.dtype), ('x', spots.dtype), ('z', spots.dtype)]
    else:
        raise ValueError("Unexpected number of columns: {}".format(spots.shape[1]))
    return np.core.records.fromarrays(spots.T, dtype=dtype)


# Maximum separation, in nanometres, at which a spot in one channel is taken to
# mark the same molecule as a spot in another. Roughly three times the 130 nm
# spot_radius used for detection.
#
# This is a fixed value rather than detect_spots_colocalization's automatic one.
# The automatic threshold is fitted per image by elbow detection on that image's
# own pair-distance distribution, so it moves with spot density: measured across
# 23 FOVs it ranged 630-5765 nm, and its median differed about threefold between
# EGS023 and EGS027, meaning the two datasets were being classified under
# different rules. The pair-distance distributions themselves agree closely
# between the datasets and peak at 114 nm, near spot_radius, but they decay
# smoothly with no plateau, so no elbow exists to find; the automatic values sat
# far out in that tail. See plotData/Figure2/diagnostic_colocalization_elbow.pdf.
COLOCALIZATION_DISTANCE_NM = 400


def safe_colocalization(spots_1, spots_2, voxel_size,
                        threshold=COLOCALIZATION_DISTANCE_NM):
    """
    Wrapper for detect_spots_colocalization that handles empty arrays.

    Pairs are assigned by optimal one-to-one matching on the nanometre distance
    matrix, then any pair separated by more than `threshold` is dropped.

    Parameters
    ----------
    spots_1 : np.ndarray
        First spots array (n_spots, 2)
    spots_2 : np.ndarray
        Second spots array (n_spots, 2)
    voxel_size : tuple
        Voxel dimensions for distance calculation
    threshold : float
        Maximum separation in nanometres for a pair to count as colocalized

    Returns
    -------
    tuple
        (coloc_1, coloc_2, distances) - colocalized spots from each array
    """
    empty = lambda s: np.empty((0, s.shape[1] if s.ndim > 1 else 2), dtype=s.dtype)

    # A fixed threshold needs no distribution to fit, so the only case that
    # cannot be handled is an empty array. The previous minimum of three spots
    # per channel existed because the automatic threshold was fitted from the
    # pair distances and needed enough of them; with an explicit threshold that
    # minimum would discard real pairs in sparse cells and push their nef spots
    # into MS.
    if spots_1.size == 0 or spots_2.size == 0:
        return empty(spots_1), empty(spots_2), np.array([])

    return multistack.detect_spots_colocalization(
        spots_1=spots_1,
        spots_2=spots_2,
        voxel_size=voxel_size,
        threshold=threshold,
        return_indices=False
    )


def safe_decompose_dense(image, spots, voxel_size, spot_radius, alpha, beta, gamma):
    """
    Wrapper for decompose_dense that falls back to the undecomposed spots when the
    reference spot cannot be modelled.

    decompose_dense returns early for an empty spot array and for an all-zero
    reference spot, but a very small non-empty array still reaches modelize_spot,
    where curve_fit raises RuntimeError instead of returning. A channel with almost
    no signal inside the cell masks lands in exactly that gap.

    Parameters
    ----------
    image : np.ndarray
        Maximum intensity projection of the channel
    spots : np.ndarray
        Detected spots, shape (n_spots, 2)
    voxel_size : tuple
        Voxel dimensions
    spot_radius : tuple
        Expected spot radius
    alpha, beta, gamma : float
        decompose_dense parameters

    Returns
    -------
    np.ndarray
        Spots after dense decomposition, or the input spots if it could not run
    """
    if spots.size == 0:
        return spots

    try:
        spots_post_decomposition, _, _ = detection.decompose_dense(
            image=image,
            spots=spots,
            voxel_size=voxel_size,
            spot_radius=spot_radius,
            alpha=alpha,
            beta=beta,
            gamma=gamma
        )
        return spots_post_decomposition
    except RuntimeError as e:
        print("  Warning: dense decomposition failed on {} spots ({}), "
              "using undecomposed spots".format(spots.shape[0], e))
        return spots


def robust_remove_transcription_site(spots, clusters, cell_label, image_mip, ndim=2,
                                     ts_radius=15):
    """
    Remove transcription sites using intensity-based detection (no nuclear segmentation).

    Parameters
    ----------
    spots : np.ndarray
        Spots with cluster assignments, shape (n_spots, 3) for 2D [y, x, cluster_id]
    clusters : np.ndarray
        Cluster information, shape (n_clusters, 4) for 2D [y, x, n_spots, cluster_id]
    cell_label : np.ndarray
        Cell segmentation mask
    image_mip : np.ndarray
        Maximum intensity projection of the channel (used for TS detection)
    ndim : int
        Number of spatial dimensions (2 or 3)
    ts_radius : float
        Radius in pixels for the circular TS region around brightest point (default: 15)

    Returns
    -------
    tuple
        (spots_no_ts, foci, ts_spots) - spots after TS removal, foci array, TS spots
    """
    # Handle empty inputs
    if spots.size == 0 or clusters.size == 0:
        empty_spots = np.empty((0, spots.shape[1] if spots.ndim > 1 else 3), dtype=spots.dtype)
        empty_foci = np.empty((0, 2), dtype=np.int64)
        return spots.copy() if spots.size > 0 else empty_spots, empty_foci, empty_spots

    # Find and mask the brightest point in the image (likely the transcription site)
    smoothed = ndimage.gaussian_filter(image_mip.astype(float), sigma=2)
    max_idx = np.unravel_index(np.argmax(smoothed), smoothed.shape)
    ts_center_y, ts_center_x = max_idx

    y_coords, x_coords = np.ogrid[:image_mip.shape[0], :image_mip.shape[1]]
    distance_from_center = np.sqrt((y_coords - ts_center_y) ** 2 + (x_coords - ts_center_x) ** 2)
    ts_mask = (distance_from_center <= ts_radius).astype(np.uint8)

    # Find clusters within the TS mask
    cluster_coords = clusters[:, :2]
    clusters_in_ts, clusters_outside_ts = multistack.identify_objects_in_region(
        ts_mask, cluster_coords, ndim=2
    )

    # Get cluster IDs in the TS region
    ts_cluster_ids = []
    for ts_cluster_coord in clusters_in_ts:
        for i in range(clusters.shape[0]):
            if np.allclose(clusters[i, :2], ts_cluster_coord):
                cluster_id = int(clusters[i, -1]) if clusters.shape[1] > 3 else i
                ts_cluster_ids.append(cluster_id)
                break

    cluster_col = 2 if ndim == 2 else 3

    # Separate spots into TS and non-TS based on cluster membership
    if len(ts_cluster_ids) > 0:
        ts_spot_mask = np.isin(spots[:, cluster_col], ts_cluster_ids)
        spots_no_ts = spots[~ts_spot_mask]
        ts_spots = spots[ts_spot_mask]
        foci = np.array([[ts_center_y, ts_center_x]], dtype=np.int64)
    else:
        spots_no_ts = spots.copy()
        ts_spots = np.empty((0, spots.shape[1]), dtype=spots.dtype)
        foci = np.empty((0, 2), dtype=np.int64)

    print("  Intensity-based TS removal: TS center at ({}, {}), radius={}, {} clusters in TS, {} spots removed".format(
        ts_center_y, ts_center_x, ts_radius, len(ts_cluster_ids), ts_spots.shape[0]))

    return spots_no_ts, foci, ts_spots


def load_and_stack_fov(path_input, fov, z_slices, channels, opt_prefix="Camera - "):
    """
    Load and stack images for a single FOV.

    Parameters
    ----------
    path_input : str
        Directory containing TIF files
    fov : str
        FOV identifier
    z_slices : list
        List of z-slice identifiers (e.g., ['Z01', 'Z02', ...])
    channels : list
        List of channel names
    opt_prefix : str
        Optical configuration prefix in filename

    Returns
    -------
    np.ndarray
        Stacked image array with shape [n_channels, n_z, height, width]
    """
    stacks = []

    for channel in channels:
        channel_stack = []
        for z in z_slices:
            pattern_str = "{}.*{}{}.*{}.*\\.tif$".format(
                re.escape(fov),
                re.escape(opt_prefix),
                re.escape(channel),
                re.escape(z)
            )
            pattern = re.compile(pattern_str)

            matching_files = [f for f in os.listdir(path_input) if pattern.search(f)]

            if len(matching_files) == 0:
                raise FileNotFoundError("No file found for {} {} {} in {}".format(
                    fov, channel, z, path_input))
            elif len(matching_files) > 1:
                print("Warning: Multiple files match {} {} {}, using first".format(fov, channel, z))

            filepath = os.path.join(path_input, matching_files[0])
            img = stack.read_image(filepath)
            channel_stack.append(img)

        # Stack z-slices for this channel
        channel_array = np.stack(channel_stack, axis=0)
        stacks.append(channel_array)

    # Stack all channels
    return np.stack(stacks, axis=0)


# ============================================================================
# DIAGNOSTIC VISUALIZATION FUNCTIONS
# ============================================================================

def save_spot_diagnostic_images(path_output, pol_mip, env_mip, nef_mip,
                                pol_no_ts, env_no_ts, nef_no_ts,
                                cell_label, fov, ts_radius=15):
    """
    Save diagnostic images showing detected spots overlaid on MIPs.
    """
    # Create figure with 3 subplots (one per channel)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    fig.suptitle('Detected Spots (after TS removal) - {} [Manual Segmentation]'.format(fov), fontsize=14)

    # Channel names and colors
    channel_info = [
        ('pol (Red)', pol_mip, pol_no_ts, 'cyan'),
        ('env (Gold70)', env_mip, env_no_ts, 'magenta'),
        ('nef (FarRed)', nef_mip, nef_no_ts, 'yellow')
    ]

    for ax, (name, mip, spots, color) in zip(axes, channel_info):
        # Display MIP with auto-contrast
        vmin, vmax = np.percentile(mip, [1, 99.5])
        ax.imshow(mip, cmap='gray', vmin=vmin, vmax=vmax)

        # Overlay spots
        if spots.size > 0:
            # spots are in (y, x) format
            ax.scatter(spots[:, 1], spots[:, 0], c=color, s=8, alpha=0.7,
                       marker='o', linewidths=0.5, edgecolors='white')

        # Add cell boundary contours if cells exist
        cell_mask = cell_label > 0
        if np.any(cell_mask):
            ax.contour(cell_mask, colors='lime', linewidths=0.5, alpha=0.5)

        # Show circular TS region around brightest point
        smoothed = ndimage.gaussian_filter(mip.astype(float), sigma=2)
        max_idx = np.unravel_index(np.argmax(smoothed), smoothed.shape)
        ts_center_y, ts_center_x = max_idx

        # Draw circle around TS center
        circle = plt.Circle((ts_center_x, ts_center_y), ts_radius,
                            fill=False, color='red', linewidth=2, linestyle='--')
        ax.add_patch(circle)
        # Mark the center
        ax.plot(ts_center_x, ts_center_y, 'r+', markersize=10, markeredgewidth=2)

        ax.set_title('{}\n{} spots'.format(name, spots.shape[0] if spots.size > 0 else 0))
        ax.axis('off')

    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor='none', edgecolor='lime', label='Cell boundary (manual)'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                   markersize=8, label='Detected spots'),
        plt.Line2D([0], [0], marker='o', color='red', markerfacecolor='none',
                   markersize=10, linestyle='--', label='TS region (r={})'.format(ts_radius))
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3, fontsize=10)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.12)

    # Save figure
    output_path = os.path.join(path_output, 'spot_detection_diagnostic.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)


# ============================================================================
# PROCESSING FUNCTIONS
# ============================================================================

# False-discovery rate, in percent, for the env-intensity rescue below. The cut is
# set at the (100 - this)th percentile of the env signal measured at random
# positions inside the cell masks, so this many percent of arbitrary in-cell
# locations pass by chance, per FOV, by construction.
ENV_RESCUE_FP_PERCENT = 1.0


def signal_over_background(image, points, radius_px):
    """
    Peak signal within radius_px of each point, in units of the image background MAD.

    Parameters
    ----------
    image : np.ndarray
        Channel MIP
    points : np.ndarray
        Coordinates, shape (n, 2), as (y, x)
    radius_px : int
        Search radius in pixels

    Returns
    -------
    np.ndarray
        (peak - background median) / background MAD, one value per point
    """
    if len(points) == 0:
        return np.array([])
    image_float = image.astype(np.float64)
    background = np.median(image_float)
    mad = max(np.median(np.abs(image_float - background)), 1e-9)
    height, width = image.shape
    peaks = np.empty(len(points))
    for i, (y, x) in enumerate(points):
        y0 = max(0, int(y) - radius_px)
        y1 = min(height, int(y) + radius_px + 1)
        x0 = max(0, int(x) - radius_px)
        x1 = min(width, int(x) + radius_px + 1)
        peaks[i] = image_float[y0:y1, x0:x1].max()
    return (peaks - background) / mad


def env_rescue_cut(env_mip, cell_label, radius_px, fp_percent=ENV_RESCUE_FP_PERCENT,
                   n_sample=3000, seed=0):
    """
    Intensity above which env is taken to be present at a coordinate despite
    detect_spots not returning a spot there.

    Calibrated against a null drawn from the FOV itself: sample random positions
    inside the cell masks, measure the same peak-over-background statistic, and
    take the (100 - fp_percent)th percentile. Calibrating instead against the
    weakest detected env spot was tried and rejected -- it produced cuts ranging
    2 to 64 across FOVs and let 12-43% of random in-cell positions through,
    because the reference moves with however deep env detection happened to go.

    Returns
    -------
    float or None
        Threshold, or None if the mask is empty
    """
    ys, xs = np.where(cell_label > 0)
    if len(ys) == 0:
        return None
    rng = np.random.RandomState(seed)
    pick = rng.choice(len(ys), min(n_sample, len(ys)), replace=False)
    null_signal = signal_over_background(
        env_mip, np.stack([ys[pick], xs[pick]], axis=1), radius_px)
    return float(np.percentile(null_signal, 100.0 - fp_percent))


def process_fov(fov, path_input, path_output_base, mask_dir, channels, pipeline_params, fov_pattern):
    """
    Process a single FOV through the complete big-FISH pipeline with manual segmentation.

    Workflow:
    1. Validate all channels present
    2. Find and load manual cell segmentation mask
    3. Auto-discover z-slices
    4. Load and stack images
    5. Detect spots in each channel
    6. Cluster spots into foci
    7. Classify splicing patterns (US/SS/MS)
    8. Save results

    Parameters
    ----------
    fov : str
        FOV identifier
    path_input : str
        Input directory containing TIF files
    path_output_base : str
        Base output directory
    mask_dir : str
        Directory containing manual segmentation masks
    channels : list
        List of channel names (3 channels: Red, Gold70, FarRed)
    pipeline_params : dict
        Pipeline parameters
    fov_pattern : str
        Regex pattern for FOV matching

    Returns
    -------
    dict or None
        Processing summary, or None if FOV was skipped
    """
    print("\n" + "=" * 70)
    print("Processing FOV: {}".format(fov))
    print('=' * 70)

    # Validate channels
    all_present, missing = validate_fov_channels(path_input, fov, channels)
    if not all_present:
        print("WARNING: Skipping FOV {} - missing channels: {}".format(fov, missing))
        return None

    # Find and load manual segmentation mask
    mask_path = find_mask_for_fov(mask_dir, fov, fov_pattern)
    if mask_path is None:
        print("WARNING: Skipping FOV {} - no mask found in {}".format(fov, mask_dir))
        return None

    print("Loading manual segmentation mask: {}".format(os.path.basename(mask_path)))
    cell_label = load_cell_mask(mask_path)

    n_cells = len(np.unique(cell_label)) - 1  # subtract 1 for background
    print("Loaded mask with {} cells".format(n_cells))

    if n_cells == 0:
        print("WARNING: Skipping FOV {} - mask contains no cells".format(fov))
        return None

    # Discover z-slices
    z_slices = discover_z_slices(path_input, fov, channels)
    if z_slices is None or len(z_slices) == 0:
        print("WARNING: No z-slices found for FOV {}".format(fov))
        return None

    print("Found {} z-slices: {}".format(len(z_slices), z_slices))

    # Load images
    print("Loading images...")
    image_stack = load_and_stack_fov(path_input, fov, z_slices, channels)
    print("Image stack shape: {} [channels, z, y, x]".format(image_stack.shape))

    # Separate channels
    pol = image_stack[0]  # Red
    env = image_stack[1]  # Gold70
    nef = image_stack[2]  # FarRed

    # Create output directory
    path_output = os.path.join(path_output_base, fov)
    os.makedirs(path_output, exist_ok=True)

    # Create MIPs
    nef_mip = stack.maximum_projection(nef)
    pol_mip = stack.maximum_projection(pol)
    env_mip = stack.maximum_projection(env)

    # Verify mask dimensions match image dimensions
    if cell_label.shape != nef_mip.shape:
        print("ERROR: Mask shape {} does not match image shape {}".format(
            cell_label.shape, nef_mip.shape))
        print("Attempting to continue, but results may be incorrect...")

    # Spot detection
    print("Detecting spots...")
    spot_params = pipeline_params['spot_detection']

    spots_pol, threshold_pol = detection.detect_spots(
        images=pol_mip,
        return_threshold=True,
        voxel_size=spot_params['voxel_size'],
        spot_radius=spot_params['spot_radius']
    )
    print("pol: {} spots (threshold: {})".format(spots_pol.shape[0], threshold_pol))

    spots_env, threshold_env = detection.detect_spots(
        images=env_mip,
        return_threshold=True,
        voxel_size=spot_params['voxel_size'],
        spot_radius=spot_params['spot_radius']
    )
    print("env: {} spots (threshold: {})".format(spots_env.shape[0], threshold_env))

    spots_nef, threshold_nef = detection.detect_spots(
        images=nef_mip,
        return_threshold=True,
        voxel_size=spot_params['voxel_size'],
        spot_radius=spot_params['spot_radius']
    )
    print("nef: {} spots (threshold: {})".format(spots_nef.shape[0], threshold_nef))

    # Mask spots by cell segmentation (only keep spots within cells).
    # All three channels are restricted, not pol alone. An unrestricted channel
    # carries out-of-cell background into the shared stages below: decompose_dense
    # builds its reference spot from whatever spots it is given, detect_clusters can
    # bridge the cell boundary and robust_remove_transcription_site then drops spots
    # by cluster id rather than position, and detect_spots_colocalization fits its
    # distance threshold to the spots it is passed and applies that threshold to
    # every pair. All three feed the per-cell US/SS/MS counts, which
    # count_spots_in_cell filters by position alone.
    if cell_label.max() > 0:
        spots_pol, _ = multistack.identify_objects_in_region(cell_label, spots_pol, ndim=2)
        spots_env, _ = multistack.identify_objects_in_region(cell_label, spots_env, ndim=2)
        spots_nef, _ = multistack.identify_objects_in_region(cell_label, spots_nef, ndim=2)
        print("Masked to cells: pol={}, env={}, nef={}".format(
            spots_pol.shape[0], spots_env.shape[0], spots_nef.shape[0]))

    # Cluster detection
    print("Decomposing dense regions and clustering foci...")
    cluster_params = pipeline_params['clustering']

    # pol channel
    pol_post_decomposition = safe_decompose_dense(
        image=pol_mip,
        spots=spots_pol,
        voxel_size=spot_params['voxel_size'],
        spot_radius=spot_params['spot_radius'],
        alpha=cluster_params['alpha'],
        beta=cluster_params['beta'],
        gamma=cluster_params['gamma']
    )
    pol_post_clustering, clusters_pol = detection.detect_clusters(
        spots=pol_post_decomposition,
        voxel_size=spot_params['voxel_size'],
        radius=cluster_params['cluster_radius'],
        nb_min_spots=cluster_params['nb_min_spots']
    )

    # env channel
    env_post_decomposition = safe_decompose_dense(
        image=env_mip,
        spots=spots_env,
        voxel_size=spot_params['voxel_size'],
        spot_radius=spot_params['spot_radius'],
        alpha=cluster_params['alpha'],
        beta=cluster_params['beta'],
        gamma=cluster_params['gamma']
    )
    env_post_clustering, clusters_env = detection.detect_clusters(
        spots=env_post_decomposition,
        voxel_size=spot_params['voxel_size'],
        radius=cluster_params['cluster_radius'],
        nb_min_spots=cluster_params['nb_min_spots']
    )

    # nef channel
    nef_post_decomposition = safe_decompose_dense(
        image=nef_mip,
        spots=spots_nef,
        voxel_size=spot_params['voxel_size'],
        spot_radius=spot_params['spot_radius'],
        alpha=cluster_params['alpha'],
        beta=cluster_params['beta'],
        gamma=cluster_params['gamma']
    )
    nef_post_clustering, clusters_nef = detection.detect_clusters(
        spots=nef_post_decomposition,
        voxel_size=spot_params['voxel_size'],
        radius=cluster_params['cluster_radius'],
        nb_min_spots=cluster_params['nb_min_spots']
    )

    print("pol: {} clusters".format(clusters_pol.shape[0] if clusters_pol.size > 0 else 0))
    print("env: {} clusters".format(clusters_env.shape[0] if clusters_env.size > 0 else 0))
    print("nef: {} clusters".format(clusters_nef.shape[0] if clusters_nef.size > 0 else 0))

    # Transcription site removal
    print("Removing transcription sites...")
    ts_params = pipeline_params.get('transcription_site', {'ts_radius': 15})
    ts_radius = ts_params.get('ts_radius', 15)

    pol_no_ts, pol_foci, pol_ts = robust_remove_transcription_site(
        pol_post_clustering, clusters_pol, cell_label, pol_mip, ndim=2, ts_radius=ts_radius
    )
    env_no_ts, env_foci, env_ts = robust_remove_transcription_site(
        env_post_clustering, clusters_env, cell_label, env_mip, ndim=2, ts_radius=ts_radius
    )
    nef_no_ts, nef_foci, nef_ts = robust_remove_transcription_site(
        nef_post_clustering, clusters_nef, cell_label, nef_mip, ndim=2, ts_radius=ts_radius
    )

    print("After TS removal: pol={}, env={}, nef={}".format(
        pol_no_ts.shape[0] if pol_no_ts.size > 0 else 0,
        env_no_ts.shape[0] if env_no_ts.size > 0 else 0,
        nef_no_ts.shape[0] if nef_no_ts.size > 0 else 0
    ))

    # Colocalization analysis
    print("Performing colocalization analysis...")
    coloc_params = pipeline_params['colocalization']
    voxel_size = coloc_params['voxel_size']
    # search radius for the env rescue, the colocalization distance in pixels
    coloc_radius_px = int(round(COLOCALIZATION_DISTANCE_NM / float(voxel_size[0])))

    # Extract just coordinates (first 2 columns) for colocalization
    pol_coords = pol_no_ts[:, :2] if pol_no_ts.size > 0 else np.empty((0, 2), dtype=np.int64)
    env_coords = env_no_ts[:, :2] if env_no_ts.size > 0 else np.empty((0, 2), dtype=np.int64)
    nef_coords = nef_no_ts[:, :2] if nef_no_ts.size > 0 else np.empty((0, 2), dtype=np.int64)

    # Classify spots: US (unspliced), SS (singly spliced), MS (multiply spliced)
    #   US (unspliced): nef + env + pol (full-length transcript, all three probes bind)
    #   SS (singly spliced): nef + env only (pol intron removed)
    #   MS (multiply spliced): nef only (both introns removed)
    print("Classifying splicing patterns...")

    # Step 1: Find nef spots colocalized with env (these are US + SS combined)
    USSS_nef, USSS_env, _ = safe_colocalization(nef_coords, env_coords, voxel_size)

    # Step 2: From nef+env spots, find those also colocalized with pol (these are US)
    US_nef, US_pol, _ = safe_colocalization(USSS_nef, pol_coords, voxel_size)

    print("Colocalization: nef+env={}, nef+env+pol={}".format(
        USSS_nef.shape[0], US_nef.shape[0]))

    # Step 3: SS = nef+env spots that are NOT also colocalized with pol
    USSS_struct = spots_to_structured(USSS_nef) if USSS_nef.size > 0 else spots_to_structured(
        np.empty((0, 2), dtype=np.int64))
    US_struct = spots_to_structured(US_nef) if US_nef.size > 0 else spots_to_structured(
        np.empty((0, 2), dtype=np.int64))

    SS_mask = ~np.isin(USSS_struct, US_struct)
    SS_spots = USSS_nef[SS_mask] if USSS_nef.size > 0 else np.empty((0, 2), dtype=np.int64)

    # Step 4: MS = nef spots that are NOT colocalized with env
    nef_struct = spots_to_structured(nef_coords) if nef_coords.size > 0 else spots_to_structured(
        np.empty((0, 2), dtype=np.int64))

    MS_mask = ~np.isin(nef_struct, USSS_struct)
    MS_spots = nef_coords[MS_mask] if nef_coords.size > 0 else np.empty((0, 2), dtype=np.int64)

    # US spots are the nef coordinates from the triple colocalization
    US_spots = US_nef if US_nef.size > 0 else np.empty((0, 2), dtype=np.int64)

    print("Classification before rescue: US={}, SS={}, MS={}".format(
        US_spots.shape[0], SS_spots.shape[0], MS_spots.shape[0]))

    # Rescue MS spots where env is present but was not detected.
    #
    # Steps 1-4 gate everything on a detected env spot, so a nef spot whose env
    # signal fell below the detection threshold lands in MS no matter how bright
    # env actually is there, and no matter whether pol is present. That calls a
    # true US an MS, skipping SS entirely.
    #
    # The test is applied to every MS spot, not only those with pol: the evidence
    # is the same in both cases -- env is bright but undetected -- so testing only
    # the pol-bearing ones would hold SS fixed while US rose, an artefact of where
    # the test was applied rather than a property of the data.
    rescue_cut = env_rescue_cut(env_mip, cell_label, coloc_radius_px)
    n_rescued_US = 0
    n_rescued_SS = 0
    US_rescued_spots = np.empty((0, 2), dtype=np.int64)
    SS_rescued_spots = np.empty((0, 2), dtype=np.int64)
    if rescue_cut is not None and MS_spots.size > 0:
        env_signal = signal_over_background(env_mip, MS_spots, coloc_radius_px)
        bright = env_signal >= rescue_cut

        has_pol = np.zeros(MS_spots.shape[0], dtype=bool)
        if pol_coords.size > 0:
            pol_hits, _, _ = safe_colocalization(MS_spots, pol_coords, voxel_size)
            if pol_hits.size > 0:
                pol_struct = spots_to_structured(pol_hits)
                has_pol = np.isin(spots_to_structured(MS_spots), pol_struct)

        to_US = bright & has_pol
        to_SS = bright & ~has_pol
        n_rescued_US = int(to_US.sum())
        n_rescued_SS = int(to_SS.sum())

        if n_rescued_US or n_rescued_SS:
            US_rescued_spots = MS_spots[to_US]
            SS_rescued_spots = MS_spots[to_SS]
            US_spots = np.vstack([US_spots, US_rescued_spots]) if n_rescued_US else US_spots
            SS_spots = np.vstack([SS_spots, SS_rescued_spots]) if n_rescued_SS else SS_spots
            MS_spots = MS_spots[~bright]

    print("env rescue: cut {} at {}% FDR, {} -> US, {} -> SS".format(
        "{:.1f}".format(rescue_cut) if rescue_cut is not None else "n/a",
        ENV_RESCUE_FP_PERCENT, n_rescued_US, n_rescued_SS))
    print("Classification: US={}, SS={}, MS={}".format(
        US_spots.shape[0], SS_spots.shape[0], MS_spots.shape[0]))

    # Save diagnostic images
    save_spot_diagnostic_images(
        path_output, pol_mip, env_mip, nef_mip,
        pol_coords, env_coords, nef_coords,
        cell_label, fov, ts_radius=ts_radius
    )

    # Per-cell statistics
    cell_results = []
    if cell_label.max() > 0:
        cell_ids = [i for i in np.unique(cell_label) if i > 0]
        for cell_id in cell_ids:
            cell_mask = cell_label == cell_id
            cell_area = np.sum(cell_mask)

            # Count spots in this cell
            def count_spots_in_cell(spots, mask):
                if spots.size == 0:
                    return 0
                in_cell = mask[spots[:, 0].astype(int), spots[:, 1].astype(int)]
                return np.sum(in_cell)

            n_US = count_spots_in_cell(US_spots, cell_mask)
            n_US_rescued = count_spots_in_cell(US_rescued_spots, cell_mask)
            n_SS_rescued = count_spots_in_cell(SS_rescued_spots, cell_mask)
            n_SS = count_spots_in_cell(SS_spots, cell_mask)
            n_MS = count_spots_in_cell(MS_spots, cell_mask)
            n_pol = count_spots_in_cell(pol_coords, cell_mask)
            n_env = count_spots_in_cell(env_coords, cell_mask)
            n_nef = count_spots_in_cell(nef_coords, cell_mask)

            cell_results.append({
                'fov': fov,
                'cell_id': cell_id,
                'cell_area': cell_area,
                'n_US': n_US,
                'n_SS': n_SS,
                'n_MS': n_MS,
                'n_US_rescued': n_US_rescued,
                'n_SS_rescued': n_SS_rescued,
                'n_pol': n_pol,
                'n_env': n_env,
                'n_nef': n_nef,
                'n_total_spots': n_US + n_SS + n_MS
            })

    # Save per-FOV CSV
    if cell_results:
        csv_path = os.path.join(path_output, 'cell_stats.csv')
        fieldnames = ['fov', 'cell_id', 'cell_area',
                      'n_US', 'n_SS', 'n_MS', 'n_US_rescued', 'n_SS_rescued',
                      'n_pol', 'n_env', 'n_nef', 'n_total_spots']
        with open(csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(cell_results)

    # Return summary
    return {
        'fov': fov,
        'n_cells': n_cells,
        'cell_results': cell_results,
        'fov_totals': {
            'n_US': US_spots.shape[0],
            'n_SS': SS_spots.shape[0],
            'n_MS': MS_spots.shape[0]
        }
    }


def run_pipeline(path_input, path_output_base, mask_dir, channels, pipeline_params, fov_list=None, fov_pattern=None):
    """
    Run the pipeline on all FOVs.

    Parameters
    ----------
    path_input : str
        Input directory
    path_output_base : str
        Output directory
    mask_dir : str
        Directory containing manual segmentation masks
    channels : list
        Channel names (3 channels: Red, Gold70, FarRed)
    pipeline_params : dict
        Pipeline parameters
    fov_list : list, optional
        Specific FOVs to process
    fov_pattern : str, optional
        Regex pattern for FOV discovery

    Returns
    -------
    tuple
        (results, errors, skipped)
    """
    # Validate channel count (should be 3 for manual pipeline without DAPI)
    if len(channels) != 3:
        raise ValueError("Expected 3 channels (Red, Gold70, FarRed), got {}: {}".format(len(channels), channels))

    # Set default FOV pattern
    if fov_pattern is None:
        fov_pattern = r'(\d{8}_[A-Z0-9]+_[a-z]+_slide\d+_\d+(?:_XY\d+)?)'

    # Discover FOVs if not provided
    if fov_list is None:
        fov_list = discover_fovs(path_input, fov_pattern)

    print("\n" + "=" * 70)
    print("big-FISH Multi-FOV Pipeline (Manual Segmentation)")
    print('=' * 70)
    print("Input directory: {}".format(path_input))
    print("Mask directory: {}".format(mask_dir))
    print("Output directory: {}".format(path_output_base))
    print("Channels: {}".format(channels))
    print("FOV pattern: {}".format(fov_pattern))
    print("Found {} FOVs to process".format(len(fov_list)))
    print('=' * 70)

    results = []
    errors = []
    skipped = []

    for i, fov in enumerate(fov_list, 1):
        print("\n\nProcessing FOV {}/{}".format(i, len(fov_list)))
        try:
            result = process_fov(fov, path_input, path_output_base, mask_dir, channels, pipeline_params, fov_pattern)
            if result is None:
                skipped.append(fov)
            else:
                results.append(result)
        except Exception as e:
            error_msg = "Error processing FOV {}: {}".format(fov, str(e))
            print("\nERROR: {}".format(error_msg))
            traceback.print_exc()
            errors.append({'fov': fov, 'error': str(e)})
            continue

    # Print summary
    print("\n" + "=" * 70)
    print("PIPELINE SUMMARY")
    print("=" * 70)
    print("Total FOVs found: {}".format(len(fov_list)))
    print("Successfully processed: {}".format(len(results)))
    print("Skipped (missing channels/mask): {}".format(len(skipped)))
    print("Errors: {}".format(len(errors)))

    if results:
        # Collect all cell results across FOVs
        all_cell_results = []
        for r in results:
            all_cell_results.extend(r.get('cell_results', []))

        # Save combined cell statistics CSV (append if exists, create if not)
        if all_cell_results:
            combined_csv_path = os.path.join(path_output_base, 'all_cells_stats.csv')
            fieldnames = ['fov', 'cell_id', 'cell_area',
                          'n_US', 'n_SS', 'n_MS', 'n_US_rescued', 'n_SS_rescued',
                          'n_pol', 'n_env', 'n_nef', 'n_total_spots']

            # Check if file exists to determine write mode and header behavior
            file_exists = os.path.isfile(combined_csv_path)

            with open(combined_csv_path, 'a', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                # Only write header if file is new
                if not file_exists:
                    writer.writeheader()
                writer.writerows(all_cell_results)

            if file_exists:
                print("\nAppended {} rows to: {}".format(len(all_cell_results), combined_csv_path))
            else:
                print("\nCreated combined cell statistics: {}".format(combined_csv_path))

        # Print per-cell summary
        print("\n" + "-" * 70)
        print("PER-CELL RESULTS")
        print("-" * 70)
        print("{:<40} {:>6} {:>6} {:>6} {:>6} {:>8}".format(
            "FOV / Cell", "US", "SS", "MS", "Total", "Area"))
        print("-" * 70)

        for r in results:
            fov = r['fov']
            cell_results = r.get('cell_results', [])

            if cell_results:
                for cell in cell_results:
                    cell_name = "{}/cell_{}".format(fov, cell['cell_id'])
                    print("{:<40} {:>6} {:>6} {:>6} {:>6} {:>8}".format(
                        cell_name[:40],
                        cell['n_US'],
                        cell['n_SS'],
                        cell['n_MS'],
                        cell['n_total_spots'],
                        cell['cell_area']
                    ))
            else:
                print("{:<40} {:>6}".format(fov, "no cells"))

        # Print FOV-level summary
        print("\n" + "-" * 70)
        print("FOV SUMMARY")
        print("-" * 70)
        total_cells = sum(r['n_cells'] for r in results)
        total_US = sum(r['fov_totals']['n_US'] for r in results)
        total_SS = sum(r['fov_totals']['n_SS'] for r in results)
        total_MS = sum(r['fov_totals']['n_MS'] for r in results)

        print("Total cells across all FOVs: {}".format(total_cells))
        print("Total spots - US: {}, SS: {}, MS: {}".format(total_US, total_SS, total_MS))

        for r in results:
            fov_totals = r.get('fov_totals', {})
            print("\n{}:".format(r['fov']))
            print("  Cells: {}".format(r['n_cells']))
            print("  FOV totals - US: {}, SS: {}, MS: {}".format(
                fov_totals.get('n_US', 0),
                fov_totals.get('n_SS', 0),
                fov_totals.get('n_MS', 0)))

    if skipped:
        print("\nSkipped FOVs ({}):".format(len(skipped)))
        for fov in skipped:
            print("  {}".format(fov))

    if errors:
        print("\nErrors encountered ({}):".format(len(errors)))
        for e in errors:
            print("  {}: {}".format(e['fov'], e['error']))

    print("=" * 70)

    return results, errors, skipped


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main entry point for command-line execution"""
    parser = argparse.ArgumentParser(
        description='big-FISH Multi-FOV Processing Pipeline with Manual Segmentation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python bigfish_manual_pipeline.py \\
    --input_dir /path/to/TIFF \\
    --output_dir /path/to/output \\
    --mask_dir /path/to/masks \\
    --channels Red Gold70 FarRed \\
    --fov_pattern r'(\d{8}_[A-Z0-9]+_[a-z]+_slide\d+_\d+(?:_XY\d+)?)'

Mask files should be named like:
  ROIs2Label_20251119_EGS023_bingo_slide6_001_XY1_Camera - FarRed_Z09.tif

The FOV pattern will be extracted from both input images and mask files to match them.
        """
    )

    parser.add_argument(
        '--input_dir',
        type=str,
        required=True,
        help='Input directory containing TIF files'
    )

    parser.add_argument(
        '--output_dir',
        type=str,
        required=True,
        help='Output directory for processed results'
    )

    parser.add_argument(
        '--mask_dir',
        type=str,
        required=True,
        help='Directory containing manual segmentation masks (ROIs2Label_*.tif files from Fiji)'
    )

    parser.add_argument(
        '--channels',
        type=str,
        nargs='+',
        default=["Red", "Gold70", "FarRed"],
        help='Channel names in order (default: Red Gold70 FarRed). Do not include Blue/DAPI.'
    )

    parser.add_argument(
        '--fov_pattern',
        type=str,
        default=r'(\d{8}_[A-Z0-9]+_[a-z]+_slide\d+_\d+(?:_XY\d+)?)',
        help='Regex pattern for FOV discovery (default matches format like 20251119_EGS023_bingo_slide6_001_XY1)'
    )

    args = parser.parse_args()

    # Validate directories
    if not os.path.isdir(args.input_dir):
        print("ERROR: Input directory does not exist: {}".format(args.input_dir))
        sys.exit(1)

    if not os.path.isdir(args.mask_dir):
        print("ERROR: Mask directory does not exist: {}".format(args.mask_dir))
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    # Pipeline parameters
    pipeline_params = {
        'spot_detection': {
            'voxel_size': (107, 107),
            'spot_radius': (130, 130)
        },
        'clustering': {
            'alpha': 0.7,
            'beta': 1,
            'gamma': 5,
            'cluster_radius': 600,
            'nb_min_spots': 6
        },
        'colocalization': {
            'voxel_size': (107, 107)
        },
        'transcription_site': {
            'ts_radius': 15
        }
    }

    # Run pipeline
    results, errors, skipped = run_pipeline(
        path_input=args.input_dir,
        path_output_base=args.output_dir,
        mask_dir=args.mask_dir,
        channels=args.channels,
        pipeline_params=pipeline_params,
        fov_pattern=args.fov_pattern
    )

    # Exit with appropriate code
    if errors:
        sys.exit(1)
    else:
        sys.exit(0)


if __name__ == "__main__":
    main()