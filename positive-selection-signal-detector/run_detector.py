import argparse
import sys
import os


def build_parser():
    parser = argparse.ArgumentParser(
        description="Positive Selection Signal Detector: dN/dS analysis using Nei-Gojobori 1986 method."
    )
    parser.add_argument("--demo", action="store_true",
                        help="Run with synthetic data (no input files needed).")
    parser.add_argument("--input", type=str, default=None,
                        help="Path to a directory of per-gene FASTA alignments or a single multi-gene FASTA.")
    parser.add_argument("--output-dir", dest="output_dir", type=str, default="output_detector",
                        help="Directory for output files.")
    parser.add_argument("--no-plot", dest="no_plot", action="store_true",
                        help="Skip generating plots.")
    parser.add_argument("--plot-format", dest="plot_format", choices=["png", "svg"],
                        default="png", help="Output format for plots.")
    parser.add_argument("--alpha", type=float, default=0.05,
                        help="FDR significance threshold.")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for demo mode.")
    return parser


def validate_args(args, parser):
    if not args.demo and args.input is None:
        parser.error("Provide --input or use --demo.")
    if args.input and not os.path.exists(args.input):
        parser.error(f"Input path not found: {args.input}")
    if not (0 < args.alpha < 1):
        parser.error("--alpha must be between 0 and 1.")


def main():
    parser = build_parser()
    args = parser.parse_args()
    validate_args(args, parser)

    from src.pipeline import run_pipeline
    result = run_pipeline(args)
    n_sig = sum(1 for g in result["tested_genes"] if g.get("significant", False))
    print(f"Analyzed {len(result['tested_genes'])} genes")
    print(f"Significant positive selection signals: {n_sig}")
    if n_sig:
        sig_genes = [g["gene"] for g in result["tested_genes"] if g.get("significant", False)]
        print(f"Significant genes: {', '.join(sig_genes)}")
    print(f"Outputs written to: {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
