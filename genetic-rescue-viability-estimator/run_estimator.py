import argparse
import sys
import os


def build_parser():
    parser = argparse.ArgumentParser(
        description="Genetic Rescue Viability Estimator: Wright-Fisher simulation of genetic rescue interventions."
    )
    parser.add_argument("--demo", action="store_true",
                        help="Run with synthetic data (no input files needed).")
    parser.add_argument("--vcf", type=str, default=None,
                        help="Path to input VCF file with population genotypes.")
    parser.add_argument("--output-dir", dest="output_dir", type=str, default="output_estimator",
                        help="Directory for output files.")
    parser.add_argument("--no-plot", dest="no_plot", action="store_true",
                        help="Skip generating plots.")
    parser.add_argument("--plot-format", dest="plot_format", choices=["png", "svg"],
                        default="png", help="Output format for plots.")
    parser.add_argument("--N-current", dest="N_current", type=int, default=50,
                        help="Current effective population size.")
    parser.add_argument("--N-target", dest="N_target", type=int, default=200,
                        help="Target effective population size after recovery.")
    parser.add_argument("--n-rescue", dest="n_rescue", type=int, default=5,
                        help="Number of rescue individuals to introduce.")
    parser.add_argument("--n-generations", dest="n_generations", type=int, default=50,
                        help="Number of generations to simulate.")
    parser.add_argument("--n-replicates", dest="n_replicates", type=int, default=100,
                        help="Number of simulation replicates.")
    parser.add_argument("--selection-coeff", dest="selection_coeff", type=float, default=0.1,
                        help="Selection coefficient against inbred individuals.")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility.")
    return parser


def validate_args(args, parser):
    if not args.demo and args.vcf is None:
        parser.error("Provide --vcf or use --demo.")
    if args.vcf and not os.path.isfile(args.vcf):
        parser.error(f"VCF file not found: {args.vcf}")
    if args.N_current < 1:
        parser.error("--N-current must be >= 1")
    if args.n_rescue < 0:
        parser.error("--n-rescue must be >= 0")
    if args.n_generations < 1:
        parser.error("--n-generations must be >= 1")
    if args.n_replicates < 1:
        parser.error("--n-replicates must be >= 1")


def main():
    parser = build_parser()
    args = parser.parse_args()
    validate_args(args, parser)

    from src.pipeline import run_pipeline
    result = run_pipeline(args)
    sr = result["score_result"]
    print(f"Viability Score: {sr['score']:.1f}/100  Grade: {sr['grade']}")
    print(f"Interpretation: {sr['interpretation']}")
    print(f"Initial mean F: {sr['initial_F']:.4f}  Final mean F: {sr['final_mean_F']:.4f}")
    print(f"Outputs written to: {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
