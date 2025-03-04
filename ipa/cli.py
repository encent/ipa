import argparse
from ipa import ipa, ipa_track, ipa_plot

def main():
    parser = argparse.ArgumentParser(prog="ipa", description="Interaction Pattern Aggregation (IPA)")
    parser.add_argument("--version", "-v", action="version", version="%(prog)s 0.1.0")

    subparsers = parser.add_subparsers(dest="command")

    # IPA main command
    parser.add_argument("--cool-path", "--cool_path", "-c", help="Path to the .cool file. Keep in mind that chromosome names in the .cool file and in all the files that will be used in the ipa_plot() function (e.g. TSS-TES sites, ATAC-Seq signal .bw file) should match each other.")
    parser.add_argument("--roi-path", "--roi_path", "-roi", help="Path to the annotation file with the regions of interest (e.g. TSS-TES sites). The file should be in a BED format.")
    parser.add_argument("--output-dir", "--output_dir", "-o", help="Path to the output directory which will store the output .bw file and the output plot files.")
    parser.add_argument("--bw-dir", "--bw_dir", "-bw", default=None, help="Path to the directory with the .bw files (ATAC-Seq, ChIP-Seq etc.) that will be used in the ipa_plot() function (default: None).")
    parser.add_argument("--expected", "-e", action="store_true", default=False, required=False, help="If True, generates an IPA track based on the observed over expected matrix (default: False).")
    parser.add_argument("--clr-weight-name", "--clr_weight_name", "-b", default="weight", required=False, help="The name of the column in the .cool file that contains the balancing weights (default: 'weight').")
    parser.add_argument("--min-dist", "--min_dist", type=int, default=40_000, required=False, help="Minimum distance (in bp) between two loci to consider for the IPA calculation, e.g. minimum loop size in bp. If None, restriction on minimum distance is not applied (default: 40_000).")
    parser.add_argument("--max-dist", "--max_dist", type=int, default=100_000, required=False, help="Maximum distance (in bp) between two loci to consider for the IPA calculation, e.g. maximum loop size in bp. If None, restriction on maximum distance is not applied (default: 100_000).")
    parser.add_argument("--nproc", "-np", type=int, default=4, required=False, help="Number of processes to use for the calculation of expected (default: 4).")
    parser.add_argument("--roi-start-name", "--roi_start_name", default=None, required=False, help="Alias for the start of the region of interest, e.g. TSS or loop start (default: None).")
    parser.add_argument("--roi-end-name", "--roi_end_name", default=None, required=False, help="Alias for the end of the region of interest, e.g. TES or loop end (default: None).")
    parser.add_argument("--flank", type=int, default=100_000, required=False, help="Size of the flanking regions in bp (default: 100_000).")
    parser.add_argument("--nbins", type=int, default=50, required=False, help="Number of bins for the stackup plot (default: 50).")
    parser.add_argument("--min-roi-size", "--min_roi_size", type=int, default=None, required=False, help="Minimum size of the region of interest (ROI) in bp to filter out small regions in the roi file (default: None).")
    parser.add_argument("--max-roi-size", "--max_roi_size", type=int, default=None, required=False, help="Maximum size of the region of interest (ROI) in bp to filter out large regions in the roi file (default: None).")

    # IPA track arguments
    parser_track = subparsers.add_parser("track", help="Calculate the IPA track from a .cool file")
    parser_track.add_argument("--cool-path", "--cool_path", "-c", required=True, help="Path to the .cool file. Keep in mind that chromosome names in the .cool file and in all the files that will be used in the ipa_plot() function (e.g. TSS-TES sites, ATAC-Seq signal .bw file) should match each other.")
    parser_track.add_argument("--output-dir", "--output_dir", "-o", required=True, help="Path to create the output directory which will store the output .bw file.")
    parser_track.add_argument("--expected", "-e", action="store_true", default=False, required=False, help="If True, generates an IPA track based on the observed over expected matrix (default: False).")
    parser_track.add_argument("--clr-weight-name", "--clr_weight_name", "-b", default="weight", required=False, help="The name of the column in the .cool file that contains the balancing weights (default: 'weight').")
    parser_track.add_argument("--min-dist", "--min_dist", type=int, default=40_000, required=False, help="Minimum distance (in bp) between two loci to consider for the IPA calculation, e.g. minimum loop size in bp. If None, restriction on minimum distance is not applied (default: 40_000).")
    parser_track.add_argument("--max-dist", "--max_dist", type=int, default=100_000, required=False, help="Maximum distance (in bp) between two loci to consider for the IPA calculation, e.g. maximum loop size in bp. If None, restriction on maximum distance is not applied (default: 100_000).")
    parser_track.add_argument("--nproc", "-np", type=int, default=4, required=False, help="Number of processes to use for the calculation of expected (default: 4).")

    # IPA plot arguments
    parser_plot = subparsers.add_parser("plot", help="Create an Interaction Pattern Aggregation (IPA) plot for a given region of interest (ROI) using up to two bigWig files.")
    parser_plot.add_argument("--bw-path", "--bw_path", "-bw", required=True, help="Path to the bigWig file. Keep in mind that chromosome names in the .bw files and in all the files that will be used in the ipa_plot() function (e.g. TSS-TES sites, ATAC-Seq signal .bw file) should match each other.")
    parser_plot.add_argument("--roi-path", "--roi_path", "-roi", required=True, help="Path to the annotation file with the regions of interest (e.g. TSS-TES sites). The file should be in a BED format.")
    parser_plot.add_argument("--output-dir", "--output_dir", "-o", required=True, help="Path to the output directory which will store the output plot file.")
    parser_plot.add_argument("--extra-bw-path", "--extra_bw_path", "-extra", default=None, required=False, help="Path to the second bigWig file (default: None).")
    parser_plot.add_argument("--roi-start-name", "--roi_start_name", default=None, required=False, help="Alias for the start of the region of interest, e.g. TSS or loop start (default: None).")
    parser_plot.add_argument("--roi-end-name", "--roi_end_name", default=None, required=False, help="Alias for the end of the region of interest, e.g. TES or loop end (default: None).")
    parser_plot.add_argument("--flank", type=int, default=100_000, required=False, help="Size of the flanking regions in bp (default: 100_000).")
    parser_plot.add_argument("--nbins", type=int, default=50, required=False, help="Number of bins for the stackup plot (default: 50).")
    parser_plot.add_argument("--min-roi-size", "--min_roi_size", type=int, default=None, required=False, help="Minimum size of the region of interest (ROI) in bp to filter out small regions in the roi file (default: None).")
    parser_plot.add_argument("--max-roi-size", "--max_roi_size", type=int, default=None, required=False, help="Maximum size of the region of interest (ROI) in bp to filter out large regions in the roi file (default: None).")

    args = parser.parse_args()

    # Handle the different command cases
    if args.command == "track":
        ipa_track(args.cool_path, args.output_dir, args.expected, 
                 args.clr_weight_name, args.min_dist, args.max_dist, args.nproc)
    elif args.command == "plot":
        ipa_plot(args.bw_path, args.roi_path, args.output_dir, 
                args.extra_bw_path, args.roi_start_name, args.roi_end_name,
                args.flank, args.nbins, args.min_roi_size, args.max_roi_size)
    else:
        # This is the main 'ipa' command without subcommands
        # Check that required arguments are present
        missing_args = []
        if not args.cool_path:
            missing_args.append("--cool-path/--cool_path/-c")
        if not args.roi_path:
            missing_args.append("--roi-path/--roi_path/-roi")
        if not args.output_dir:
            missing_args.append("--output-dir/--output_dir/-o")
            
        if missing_args:
            parser.error(f"the following arguments are required: {', '.join(missing_args)}")
        
        ipa(args.cool_path, args.roi_path, args.output_dir, args.bw_dir, 
           args.expected, args.clr_weight_name, args.min_dist, args.max_dist, 
           args.nproc, args.roi_start_name, args.roi_end_name, args.flank, 
           args.nbins, args.min_roi_size, args.max_roi_size)

if __name__ == "__main__":
    main()
