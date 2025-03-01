import argparse
from ipa import ipa, ipa_track, ipa_plot

def main():
    parser = argparse.ArgumentParser(prog="ipa", description="Interaction Pattern Aggregation (IPA)")
    parser.add_argument("--version", "-v", action="version", version="%(prog)s 0.1.0")

    subparsers = parser.add_subparsers(dest="command")

    # IPA arguments
    parser.add_argument("--cool_path", "-c", required=True, help="Path to the .cool file. Keep in mind that chromosome names in the .cool file and in all the files that will be used in the ipa_plot() function (e.g. TSS-TES sites, ATAC-Seq signal .bw file) should match each other.")
    parser.add_argument("--roi_path", "-roi", required=True, help="Path to the annotation file with the regions of interest (e.g. TSS-TES sites). The file should contain at least three columns: 'chrom', 'start', 'end' without header.")
    parser.add_argument("--output_dir", "-o", required=True, help="Path to the output directory which will store the output .bw file and the output plot files.")
    parser.add_argument("--bw_dir", "-bw", required=True, help="Path to the directory with the .bw files (ATAC-Seq, ChIP-Seq etc.) that will be used in the ipa_plot() function (default: None).")
    parser.add_argument("--expected", "-e", action="store_true", default=False, required=False, help="If True, generates an IPA track based on the observed over expected matrix (default: False).")
    parser.add_argument("--clr_weight_name", "-b", default="weight", required=False, help="The name of the column in the .cool file that contains the balancing weights (default: 'weight').")
    parser.add_argument("--min_dist", type=int, default=40_000, required=False, help="Minimum distance (in bp) between two loci to consider for the IPA calculation, e.g. minimum loop size in bp. If None, restriction on minimum distance is not applied (default: 40_000).")
    parser.add_argument("--max_dist", type=int, default=100_000, required=False, help="Maximum distance (in bp) between two loci to consider for the IPA calculation, e.g. maximum loop size in bp. If None, restriction on maximum distance is not applied (default: 100_000).")
    parser.add_argument("--nproc", "-np", type=int, default=4, required=False, help="Number of processes to use for the calculation of expected (default: 4).")
    parser.add_argument("--roi_start_name", default=None, required=False, help="Alias for the start of the region of interest, e.g. TSS or loop start (default: None).")
    parser.add_argument("--roi_end_name", default=None, required=False, help="Alias for the end of the region of interest, e.g. TES or loop end (default: None).")
    parser.add_argument("--flank", type=int, default=100_000, required=False, help="Size of the flanking regions in bp (default: 100_000).")
    parser.add_argument("--nbins",type=int, default=50, required=False, help="Number of bins for the stackup plot (default: 50).")
    parser.add_argument("--min_roi_size", type=int, default=None, required=False, help="Minimum size of the region of interest (ROI) in bp to filter out small regions in the roi file (default: None).")
    parser.add_argument("--max_roi_size", type=int, default=None, required=False, help="Maximum size of the region of interest (ROI) in bp to filter out large regions in the roi file (default: None).")
    parser.set_defaults(func=ipa)

    # IPA track arguments
    parser_ipa_track = subparsers.add_parser("track", help="Calculate the Interaction Pattern Aggregation (IPA) track from a .cool file and save it to a .bw file.")
    parser_ipa_track.add_argument("--cool_path", "-c", required=True, help="Path to the .cool file. Keep in mind that chromosome names in the .cool file and in all the files that will be used in the ipa_plot() function (e.g. TSS-TES sites, ATAC-Seq signal .bw file) should match each other.")
    parser_ipa_track.add_argument("--output_dir", "-o", required=True, help="Path to create the output directory which will store the output .bw file.")
    parser_ipa_track.add_argument("--expected", "-e", action="store_true", default=False, required=False, help="If True, generates an IPA track based on the observed over expected matrix (default: False).")
    parser_ipa_track.add_argument("--clr_weight_name", "-b", default="weight", required=False, help="The name of the column in the .cool file that contains the balancing weights (default: 'weight').")
    parser_ipa_track.add_argument("--min_dist", type=int, default=40_000, required=False, help="Minimum distance (in bp) between two loci to consider for the IPA calculation, e.g. minimum loop size in bp. If None, restriction on minimum distance is not applied (default: 40_000).")
    parser_ipa_track.add_argument("--max_dist", type=int, default=100_000, required=False, help="Maximum distance (in bp) between two loci to consider for the IPA calculation, e.g. maximum loop size in bp. If None, restriction on maximum distance is not applied (default: 100_000).")
    parser_ipa_track.add_argument("--nproc", "-np", type=int, default=4, required=False, help="Number of processes to use for the calculation of expected (default: 4).")
    parser_ipa_track.set_defaults(func=ipa_track)

    # IPA plot arguments
    parser_ipa_plot = subparsers.add_parser("plot", help="Create an Interaction Pattern Aggregation (IPA) plot for a given region of interest (ROI) using up to two bigWig files.")
    parser_ipa_plot.add_argument("--bw_path", "-bw", required=True, help="Path to the bigWig file. Keep in mind that chromosome names in the .bw file and in all the files that will be used in the ipa_plot() function (e.g. TSS-TES sites, ATAC-Seq signal .bw file) should match each other.")
    parser_ipa_plot.add_argument("--roi_path", "-roi", required=True, help="Path to the annotation file with the regions of interest (e.g. TSS-TES sites). The file should contain at least three columns: 'chrom', 'start', 'end' without header.")
    parser_ipa_plot.add_argument("--output_dir", "-o", required=True, help="Path to the output directory which will store the output plot file.")
    parser_ipa_plot.add_argument("--extra_bw_path", "-extra", default=None, required=False, help="Path to the second bigWig file (default: None).")
    parser_ipa_plot.add_argument("--roi_start_name", default=None, required=False, help="Alias for the start of the region of interest, e.g. TSS or loop start (default: None).")
    parser_ipa_plot.add_argument("--roi_end_name", default=None, required=False, help="Alias for the end of the region of interest, e.g. TES or loop end (default: None).")
    parser_ipa_plot.add_argument("--flank", type=int, default=100_000, required=False, help="Size of the flanking regions in bp (default: 100_000).")
    parser_ipa_plot.add_argument("--nbins",type=int, default=50, required=False, help="Number of bins for the stackup plot (default: 50).")
    parser_ipa_plot.add_argument("--min_roi_size", type=int, default=None, required=False, help="Minimum size of the region of interest (ROI) in bp to filter out small regions in the roi file (default: None).")
    parser_ipa_plot.add_argument("--max_roi_size", type=int, default=None, required=False, help="Maximum size of the region of interest (ROI) in bp to filter out large regions in the roi file (default: None).")
    parser_ipa_plot.set_defaults(func=ipa_plot)

    args = parser.parse_args()

    if hasattr(args, "func") and args.func == ipa:
        ipa(args.cool_path, args.roi_path, args.output_dir, args.bw_dir, args.expected, args.clr_weight_name, args.min_dist, args.max_dist, args.nproc, args.roi_start_name, args.roi_end_name, args.flank, args.nbins, args.min_roi_size, args.max_roi_size)
    elif hasattr(args, "func") and args.func == ipa_track:
        ipa_track(args.cool_path, args.output_dir, args.expected, args.clr_weight_name, args.min_dist, args.max_dist, args.nproc)
    elif hasattr(args, "func") and args.func == ipa_plot:
        ipa_plot(args.bw_path, args.roi_path, args.output_dir, args.extra_bw_path, args.roi_start_name, args.roi_end_name, args.flank, args.nbins, args.min_roi_size, args.max_roi_size)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
