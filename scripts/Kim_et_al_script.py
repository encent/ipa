"""
A Python script to reproduce results of Kim et al. (2025) related to the IPA analysis (see the figure in a README file).
"""
import argparse
import os
import glob
from tqdm import tqdm
from ipa import ipa

# Define the parameters for each species (given in bp)
species_params = {
    "Cowc": {"resolution": 400, "min_dist": 4000, "max_dist": 100_000, "flank": 10_000},
    "Tadh": {"resolution": 400, "min_dist": 4000, "max_dist": 100_000, "flank": 10_000},
    "Mlei": {"resolution": 800, "min_dist": 5000, "max_dist": 150_000, "flank": 10_000},
    "Nvec": {"resolution": 500, "min_dist": 10_000, "max_dist": 360_000, "flank": 10_000},
    "Dmel": {"resolution": 400, "min_dist": 5000, "max_dist": 250_000, "flank": 10_000},
    "Hsap": {"resolution": 5000, "min_dist": 50_000, "max_dist": 1_060_000, "flank": 100_000},
    "Sarc": {"resolution": 2800, "min_dist": 50_000, "max_dist": 4_000_000, "flank": 30_000}
}
nbins = 50

# Set up argument parser
parser = argparse.ArgumentParser(description="Reproduce results of Kim et al. (2025) related to the IPA analysis.")
parser.add_argument("--input-dir", type=str, required=True, help="Downloaded directory containing input data files (data_Kim_et_al).")
parser.add_argument("--output-dir", type=str, required=True, help="Directory to save the output results.")
args = parser.parse_args()

# Define the input and output directories
input_dir = args.input_dir
output_dir = args.output_dir

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Perform IPA analysis for each species
for species, params in tqdm(species_params.items()):
    # Get cooler filename
    clr_files = glob.glob(os.path.join(input_dir, "coolers", f"{species}*.mcool"))
    if not clr_files:
        print(f"No cooler file found for species {species}")
        continue
    clr_path = f"{clr_files[0]}::resolutions/{params['resolution']}"
    
	# Get the ROI filename
    roi_files = glob.glob(os.path.join(input_dir, "roi", f"{species}*.bed"))
    if not roi_files:
        print(f"No cooler file found for species {species}")
        continue
    roi_file = roi_files[0]

	# Get the species-specific bigwigs directory
    bw_dir = os.path.join(input_dir, "bigwigs", species)
    
    species_output_dir = os.path.join(output_dir, f"{species}_res_{params['resolution']}bp_min_dist_{params['min_dist']}bp_max_dist_{params['max_dist']}bp")
    # Ensure the output directory exists
    os.makedirs(species_output_dir, exist_ok=True)

    # Call the ipa function with the specified parameters
    ipa(
        clr_path=clr_path,
        roi_file=roi_file,
        output_dir=species_output_dir,
        bw_dir=bw_dir,
        expected=False,
        clr_weight_name='weight',
        min_dist=params["min_dist"],
        max_dist=params["max_dist"],
        roi_start_name="TSS" if species != "Sarc" else "RSS",
        roi_end_name="TES" if species != "Sarc" else "RES",
        flank=params["flank"],
        nbins=nbins
    )

    print(f"IPA analysis completed for {species}. Results saved to {species_output_dir}.")
