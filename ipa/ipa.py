import argparse
import os
import math
from tqdm import tqdm
import warnings

import bbi
import bioframe
import cooler
import cooltools
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import numpy as np
import pandas as pd

from ipa.lib import mask_out_diagonals, fetch_cis_matrix, calculate_observed_over_expected_matrix, warning_chromnames, create_stackup_plot, filter_regions


def ipa_track(clr_path, output_dir, expected=False, clr_weight_name='weight', min_dist=40_000, max_dist=100_000, nproc=4):
    """
    Calculate the Interaction Pattern Aggregation track (IPA) from a .cool file and save it to a .bw file.

    Args:
        clr_path: Path to the .cool file. Keep in mind that chromosome names in the .cool file and in all the files that will be used in the ipa_plot() function (e.g. TSS-TES sites, ATAC-Seq signal .bw file) should match each other.
        output_dir: Path to create the output directory which will store the output .bw file.
        expected: If True, generates an IPA track based on the observed over expected matrix (default: False).
        clr_weight_name: The name of the column in the .cool file that contains the balancing weights (default: 'weight').
        min_dist: Minimum distance (in bp) between two loci to consider for the IPA calculation, e.g. minimum loop size in bp. If None, restriction on minimum distance is not applied (default: 40_000).
        max_dist: Maximum distance (in bp) between two loci to consider for the IPA calculation, e.g. maximum loop size in bp. If None, restriction on maximum distance is not applied (default: 100_000).
        nproc: Number of processes to use for the calculation of expected (default: 4).
    """
    # Create output directory
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    else:
        warnings.warn(f"Directory {output_dir} already exists. The content of the directory could be overwritten.")

    # Read cool file
    clr = cooler.Cooler(clr_path)
    resolution, chromnames, chromsizes = clr.binsize, clr.chromnames, clr.chromsizes
    view_df = bioframe.make_viewframe(chromsizes)
    bins = clr.bins()[:][['chrom', 'start', 'end']]

    # Warning if any of chromosome names do not start with 'chr'
    warning_chromnames(chromnames, clr_path)

    # Get min and max diagonals
    get_min_diag = lambda dist: math.floor(dist / resolution) if dist is not None else 0
    get_max_diag = lambda dist: math.ceil(dist / resolution) if dist is not None else None
    min_diag, max_diag = get_min_diag(min_dist), get_max_diag(max_dist)
    
    # Create ipa track for each individual chromosome
    for chrom in tqdm(chromnames):
        # Fetch cis matrix
        cis_matrix = fetch_cis_matrix(clr, chrom, clr_weight_name)

        # Mask out diagonals in contact matrix based on min and max distance
        mask_out_diagonals(cis_matrix, min_diag, max_diag)

        # Expected calculation (optional)
        if expected:
            cis_matrix = calculate_observed_over_expected_matrix(cis_matrix, clr, chrom, view_df, min_diag, clr_weight_name, nproc)

        # Calculate average statistics
        ipa_track = np.nansum(cis_matrix, axis=1)
        ipa_track[ipa_track == 0.] = np.nan
        del cis_matrix

        # Final dataframe arrangement
        idx = bins.index[bins['chrom'] == chrom]
        bins.loc[idx, 'ipa'] = ipa_track
    
    # Save the ipa track to a `output_bw_file` file
    bioframe.to_bigwig(bins, chromsizes, os.path.join(output_dir, "ipa_track.bw"), value_field="ipa")

def ipa_plot(bw_file, roi_file, output_dir, extra_bw_file=None, roi_start_name=None, roi_end_name=None, flank=100_000, nbins=50, min_roi_size=None, max_roi_size=None):
    """
    Create an Interaction Pattern Aggregation (IPA) plot for a given region of interest (ROI) using up to two bigWig files.

    Args:
        bw_file: Path to the bigWig file. Keep in mind that chromosome names in the .bw file and in all the files that will be used in the ipa_plot() function (e.g. TSS-TES sites, ATAC-Seq signal .bw file) should match each other.
        roi_file: Path to the annotation file with the regions of interest (e.g. TSS-TES sites). The file should be in a BED format.
        output_dir: Path to the output directory which will store the output plot file.
        extra_bw_file: (optional) Path to the second bigWig file (default: None).
        roi_start_name: Alias for the start of the region of interest, e.g. TSS or loop start (default: None).
        roi_end_name: Alias for the end of the region of interest, e.g. TES or loop end (default: None).
        flank: Flank size in bp (default: 100_000).
        nbins: Number of bins to split the ROI into (default: 50).
        min_roi_size: Minimum size of the region of interest in bp to filter out small regions in the roi file (default: None).
        max_roi_size: Maximum size of the region of interest in bp to filter out large regions in the roi file (default: None).
    """
    # Check if the output directory exists
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    else:
        warnings.warn(f"Directory {output_dir} already exists. The content of the directory could be overwritten.")

    # Read the annotation file with the regions of interest (e.g. TSS-TES sites)
    roi_df = bioframe.read_table(roi_file, schema='bed')

    # Check if the annotation file contains the required columns
    assert all(column in roi_df.columns for column in ['chrom', 'start', 'end']), f"File {roi_file} must contain at least these three columns: {'chrom', 'start', 'end'}"
    
    # Warning if any of chromosome names do not start with 'chr'
    warning_chromnames(roi_df['chrom'], roi_file)

    # Region filtering based on the size
    roi_df = filter_regions(roi_df, min_roi_size, max_roi_size)

    # Create a stackup plot
    stackup_concat = create_stackup_plot(bw_file, roi_df, flank=flank, nbins=nbins)

    # IPA plot
    f, ax1 = plt.subplots(figsize=[15, 5])

    # Plot the first dataset
    line1, = ax1.plot(np.nanmean(stackup_concat, axis=0), label=os.path.basename(bw_file))
    ax1.set_xlabel(f'Distance from {roi_start_name}/{roi_end_name}, kbp')
    ax1.set_ylabel(os.path.basename(bw_file))
    ax1.set_title(os.path.basename(os.path.normpath(output_dir)))

    # Make x-axis ticks
    x_ticks = list(np.arange(0, nbins + 10, 10)) + list(np.arange(2 * nbins - 1, 3 * nbins, 10))
    x_labels = [int(x) for x in list(np.linspace(-flank // 1000, 0, 6))[:-1]] + [roi_start_name, roi_end_name] + [int(x) for x in list(np.linspace(0, flank // 1000, 6))[1:]]

    ax1.set_xticks(x_ticks)
    ax1.set_xticklabels(x_labels)

    # Create a second y-axis (optional)
    if extra_bw_file is not None:
        stackup_concat_2 = create_stackup_plot(extra_bw_file, roi_df, flank=flank, nbins=nbins)
        ax2 = ax1.twinx()
        line2, = ax2.plot(np.nanmean(stackup_concat_2, axis=0), color='r', label=os.path.basename(extra_bw_file))
        ax2.set_ylabel(os.path.basename(extra_bw_file))
        lines = [line1, line2]
        output_plot_filename = os.path.join(output_dir, f"{os.path.basename(bw_file).split('.')[0]}_{os.path.basename(extra_bw_file).split('.')[0]}.png")
    else:
        lines = [line1]
        output_plot_filename = os.path.join(output_dir, f"{os.path.basename(bw_file).split('.')[0]}.png")
    
    # Combine legends
    labels = [line.get_label() for line in lines]
    ax1.legend(lines, labels, loc='best')

    # Save the plot
    f.savefig(output_plot_filename, dpi=300, bbox_inches='tight')
    plt.close(f)

def ipa(clr_path, roi_file, output_dir, bw_dir=None, expected=False, clr_weight_name='weight', min_dist=40_000, max_dist=100_000, nproc=4, roi_start_name=None, roi_end_name=None, flank=100_000, nbins=50, min_roi_size=None, max_roi_size=None):
    """
    Run the Interaction Pattern Aggregation analysis (IPA).
    It consists of two steps:
    1. Create a .bw file with Interaction Pattern Aggregation track based on the .cool file.
    2. Create a stackup plot based on the .bw files.
    Args:
        clr_path: Path to the .cool file. Keep in mind that chromosome names in the .cool file and in all the files that will be used in the ipa_plot() function (e.g. TSS-TES sites, ATAC-Seq signal .bw file) should match each other.
        roi_file: Path to the annotation file with the regions of interest (e.g. TSS-TES sites). The file should be in a BED format.
        output_dir: Path to the output directory which will store the output ipa track and plot files.
        bw_dir: Path to the directory with the .bw files (ATAC-Seq, ChIP-Seq etc.) that will be used in the ipa_plot() function (default: None).
        expected: If True, generates an IPA track based on the observed over expected matrix (default: False).
        clr_weight_name: The name of the column in the .cool file that contains the balancing weights (default: 'weight').
        min_dist: Minimum distance (in bp) between two loci to consider for the IPA calculation, e.g. minimum loop size in bp. If None, restriction on minimum distance is not applied (default: 40_000).
        max_dist: Maximum distance (in bp) between two loci to consider for the IPA calculation, e.g. maximum loop size in bp. If None, restriction on maximum distance is not applied (default: 100_000).
        nproc: Number of processes to use for the calculation of expected (default: 4).
        roi_start_name: Alias for the start of the region of interest, e.g. TSS or loop start (default: None).
        roi_end_name: Alias for the end of the region of interest, e.g. TES or loop end (default: None).
        flank: Flank size in bp (default: 100_000).
        nbins: Number of bins to split the ROI into (default: 50).
        min_roi_size: Minimum size of the region of interest in bp to filter out small regions in the roi file (default: None).
        max_roi_size: Maximum size of the region of interest in bp to filter out small regions in the roi file (default: None).
    """
    # Step 1: Create a .bw file from .cool file
    ipa_track(clr_path, output_dir, expected=expected, clr_weight_name=clr_weight_name, min_dist=min_dist, max_dist=max_dist, nproc=nproc)

    # Step 2: Create a stackup plot from .bw files
    if bw_dir is None:
        ipa_plot(os.path.join(output_dir, "ipa_track.bw"), roi_file, os.path.join(output_dir, "ipa_track.png"), extra_bw_file=None, roi_start_name=roi_start_name, roi_end_name=roi_end_name, flank=flank, nbins=nbins, min_roi_size=min_roi_size, max_roi_size=max_roi_size)
    else:
        bw_files = [os.path.join(bw_dir, f) for f in os.listdir(bw_dir) if f.lower().endswith('.bw') or f.lower().endswith('.bigwig')]
        for bw_file in bw_files:
            ipa_plot(os.path.join(output_dir, "ipa_track.bw"), roi_file, output_dir, extra_bw_file=bw_file, roi_start_name=roi_start_name, roi_end_name=roi_end_name, flank=flank, nbins=nbins, min_roi_size=min_roi_size, max_roi_size=max_roi_size)
