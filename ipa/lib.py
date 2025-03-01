import warnings

import bbi
import cooler
import cooltools
import numpy as np
import pandas as pd


def mask_out_diagonals(matrix, min_diag, max_diag):
	"""
	Mask out first `min_diag` diagonals and diagonals starting 
	from `max_diag` diagonal until the end of the matrix `matrix`.
	
	Args:
		matrix: NumPy 2D array (matrix) to modify.
		min_diag: number of first diagonals to mask out.
		max_diag: number of last diagonals to mask out.
	"""
	# Define functions to mask out diagonals above and below the main diagonal
	mask_diagonal_above = lambda matrix, k: np.fill_diagonal(matrix[:, k:], np.nan) if matrix.shape[1] - k > 0 else None
	mask_diagonal_below = lambda matrix, k: np.fill_diagonal(matrix[k:, :], np.nan) if matrix.shape[0] - k > 0 else None

	# Mask out diagonals lower the `min_diag` diagonal
	for k in range(min_diag):
		# Mask out the main diagonal (k = 0)
		if k == 0:
			np.fill_diagonal(matrix, np.nan)
		else:
			# Mask out the diagonal k
			mask_diagonal_above(matrix, k)
			mask_diagonal_below(matrix, k)

	# Mask out diagonals higher the `max_diag` diagonal 
	if max_diag is not None:
		for k in range(max_diag + 1, matrix.shape[0]):
			# Mask out the diagonal k
			mask_diagonal_above(matrix, k)
			mask_diagonal_below(matrix, k)

def create_expected_matrix(expected_arr):
	"""
	Build a symmetric matrix based on expected values calculated by `cooltools.expected_cis` function.
	In the resulted expected matrix the main diagonal is filled with expected_arr[0],
	the 1st diagonal and -1st diagonal are filled with expected_arr[1], and so on.
	
	Args:
		expected_arr: 1D NumPy array of expected values calculated by `cooltools.expected_cis` function.
		
	Returns:
		A square NumPy matrix with diagonals filled symmetrically.
	"""
	n = len(expected_arr)  # The size of the square matrix (nxn)
	matrix = np.zeros((n, n), dtype=expected_arr.dtype)  # Initialize an empty square matrix

	# Fill the main diagonal (0th diagonal)
	np.fill_diagonal(matrix, expected_arr[0])
	
	# Fill the upper and lower diagonals symmetrically
	for k in range(1, n):
		# Fill k-th diagonal above the main diagonal (upper diagonal)
		if k < n:
			np.fill_diagonal(matrix[:, k:], expected_arr[k])  # Fill diagonal above the main diagonal
			np.fill_diagonal(matrix[k:, :], expected_arr[k])  # Fill corresponding lower diagonal

	return matrix

def fetch_cis_matrix(clr, chrom, clr_weight_name):
	"""
	Fetch the cis contact matrix for a given chromosome from a Cooler object.
	
	Args:
		clr: Cooler object.
		chrom: Chromosome name.
		clr_weight_name: The name of the column in the .cool file that contains the balancing weights (default: 'weight').
	
	Returns:
		A NumPy 2D array with the cis contact matrix.
	"""
	# Fetching the contact matrix based on the clr_weight_name
	if clr_weight_name is None:
		cis_matrix = clr.matrix(balance=False, sparse=True).fetch(chrom)
	else:
		cis_matrix = clr.matrix(balance=True, sparse=True).fetch(chrom)
	
	# Convert the sparse matrix to a dense NumPy array -- maybe to skip this step and put sparse=False in clr.matrix()?
	cis_matrix_np = cis_matrix.toarray()
	del cis_matrix

	# Convert the matrix to float type -- is there a way in Cooler to specify the type of the matrix?
	cis_matrix_np = cis_matrix_np.astype(float)

	return cis_matrix_np

def calculate_observed_over_expected_matrix(cis_matrix, clr, chrom, view_df, min_diag, clr_weight_name, nproc):
	"""
	Calculate the observed over expected matrix for a given chromosome.
	
	Args:
		cis_matrix: NumPy 2D array with the cis contact matrix.
		clr: Cooler object.
		chrom: Chromosome name.
		view_df: ViewFrame with the chromosome sizes.
		min_diag: Minimum diagonal to consider for the calculation.
		clr_weight_name: The name of the column in the .cool file that contains the balancing weights (default: 'weight').
		nproc: Number of processes to use for the calculation of expected (default: 4).
	
	Returns:
		A NumPy 2D array with the observed over expected matrix
	"""
	# Calculate expected
	expected_df = cooltools.expected_cis(clr, view_df=view_df[view_df['chrom'] == chrom], 
											ignore_diags=min_diag, nproc=nproc, chunksize=1_000_000, 
											clr_weight_name=clr_weight_name)
	
	# Extract the expected values from the DataFrame
	if clr_weight_name is None:
		expected_colname = "count.avg"
	else:
		expected_colname = "balanced.avg"
	expected_arr = expected_df[expected_colname]
	del expected_df

	# Obtaining expected matrix
	expected_matrix = create_expected_matrix(expected_arr)
	del expected_arr

	# Calculate observed over expected matrix
	cis_matrix = cis_matrix / expected_matrix
	del expected_matrix

	return cis_matrix

def warning_chromnames(chromnames, file_path):
	"""
	Check if all chromosome names start with 'chr' and issue a warning if not.
	
	Args:
		chromnames: List of chromosome names.
		file_path: Path to the file.
	"""
	if not all(str(chrom).startswith('chr') for chrom in chromnames):
		warnings.warn(f"Some values in the 'chrom' column of {file_path} do not start with 'chr' prefix. Keep in mind that chromosome names in the .cool file and in all the files that will be used in the ipa_plot() function (e.g. TSS-TES sites, ATAC-Seq signal .bw file) should match each other.")

def create_stackup_plot(bw_file, roi_df, flank=100_000, nbins=50):
	"""
	Create a stackup plot for a given region of interest (ROI) using a bigWig file.

	Args:
		bw_file: Path to the bigWig file.
		roi_df: DataFrame with the region of interest (ROI) information. It should contain at least 'chrom', 'start', and 'end'.
		flank: Size of the flanking regions in bp (default: 100_000).
		nbins: Number of bins to split the ROI into (default: 50).
	
	Returns:
		A NumPy 2D array with the stackup plot.
	"""
	# Fetch the stackup signal for the region of interest and flanking regions
	stackup_region = bbi.stackup(bw_file, roi_df.chrom, roi_df.start, roi_df.end, bins=nbins)
	stackup_left_flank = bbi.stackup(bw_file, roi_df.chrom, roi_df.start - flank, roi_df.start, bins=nbins)
	stackup_right_flank = bbi.stackup(bw_file, roi_df.chrom, roi_df.end, roi_df.end + flank, bins=nbins)

	# Flip the stackup signal if the strand is negative
	if 'strand' in roi_df.columns:
		stackup_region = np.asarray(list(map(lambda stackup, strand: stackup[::-1] if strand is '-' else stackup, stackup_region, roi_df['strand'])))
		stackup_left_flank_flipped = []
		stackup_right_flank_flipped = []
		for i, (s_left, s_right) in enumerate(zip(stackup_left_flank, stackup_right_flank)):
			if roi_df['strand'].iloc[i] == '-':
				stackup_left_flank_flipped.append(s_right[::-1].tolist())
				stackup_right_flank_flipped.append(s_left[::-1].tolist())
			else:
				stackup_left_flank_flipped.append(s_left.tolist())
				stackup_right_flank_flipped.append(s_right.tolist())
		
		stackup_left_flank = np.asarray(stackup_left_flank_flipped)
		stackup_right_flank = np.asarray(stackup_right_flank_flipped)

	# Concatenate the stackup signals for the left flank, region of interest, and right flank
	stackup_concat = np.hstack((stackup_left_flank, stackup_region, stackup_right_flank))

	return stackup_concat

def filter_regions(roi_df, min_roi_size=None, max_roi_size=None):
	"""
	Filter regions of interest in the roi file based on their size.

	Args:
		roi_df: DataFrame with the region of interest (ROI) information. It should contain at least 'chrom', 'start' and 'end' columns.
		min_roi_size: Minimum size of the region of interest in bp to filter out small regions in the roi file (default: None).
		max_roi_size: Maximum size of the region of interest in bp to filter out large regions in the roi file (default: None).

	Returns:
		Filtered DataFrame with the region of interest information.
	"""
	if min_roi_size is not None:
		roi_df = roi_df[roi_df['end'] - roi_df['start'] >= min_roi_size]
	if max_roi_size is not None:
		roi_df = roi_df[roi_df['end'] - roi_df['start'] <= max_roi_size]
	roi_df.reset_index(drop=True, inplace=True)

	return roi_df
