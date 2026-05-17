import argparse
import sys
import os


def build_parser():
    parser = argparse.ArgumentParser(
        description="Proxy Species Edit Burden Calculator: computes edit burden to transform proxy genome toward target."
    )
    parser.add_argument("--demo", action="store_true",
                        help="Run with synthetic genomes (no input files needed).")
    parser.add_argument("--proxy", type=str, default=None,
                        help="Path to proxy species FASTA genome.")
    parser.add_argument("--target", type=str, default=None,
                        help="Path to target (extinct) species FASTA genome.")
    parser.add_argument("--gff", type=str, default=None,
                        help="Optional GFF3 annotation file for impact annotation.")
    parser.add_argument("--output-dir", dest="output_dir", type=str, default="output_calculator",
                        help="Directory for output files.")
    parser.add_argument("--no-plot", dest="no_plot", action="store_true",
                        help="Skip generating plots.")
    parser.add_argument("--plot-format", dest="plot_format", choices=["png", "svg"],
                        default="png", help="Output format for plots.")
    parser.add_argument("--n-samples", dest="n_samples", type=int, default=1000,
                        help="Number of 500bp windows to sample for alignment.")
    parser.add_argument("--window-size", dest="window_size", type=int, default=500,
                        help="Window size for pairwise alignment.")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility.")
    return parser


def validate_args(args, parser):
    if not args.demo and (args.proxy is None or args.target is None):
        parser.error("Provide --proxy and --target, or use --demo.")
    if args.proxy and not os.path.isfile(args.proxy):
        parser.error(f"Proxy FASTA not found: {args.proxy}")
    if args.target and not os.path.isfile(args.target):
        parser.error(f"Target FASTA not found: {args.target}")
    if args.gff and not os.path.isfile(args.gff):
        parser.error(f"GFF3 file not found: {args.gff}")
    if args.n_samples < 1:
        parser.error("--n-samples must be >= 1")


def main():
    parser = build_parser()
    args = parser.parse_args()
    validate_args(args, parser)

    from src.pipeline import run_pipeline
    result = run_pipeline(args)
    b = result["burden"]
    print(f"Total variants called: {b['total_edits']}")
    print(f"Weighted edit burden: {b['weighted_burden']}")
    print(f"Normalized burden: {b['normalized_burden_per_mb']:.2f} per Mb")
    for vc, count in sorted(b["class_counts"].items()):
        print(f"  {vc}: {count}")
    print(f"Outputs written to: {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
