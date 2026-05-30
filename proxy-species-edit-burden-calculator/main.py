"""Entry-point for the Proxy Species Edit Burden Calculator.

Usage examples
--------------
Run with synthetic demo genomes:
    python main.py --demo --output-dir results/

Run with real FASTA files:
    python main.py --proxy proxy.fa --target target.fa --output-dir results/

With GFF3 annotation:
    python main.py --proxy proxy.fa --target target.fa --gff genes.gff3 --output-dir results/
"""

import sys
import os

# Ensure the package root is on the path when invoked directly.
_here = os.path.dirname(os.path.abspath(__file__))
if _here not in sys.path:
    sys.path.insert(0, _here)


def _build_parser():
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Proxy Species Edit Burden Calculator: compute the number and nature "
            "of genome edits required to convert a proxy species genome toward a "
            "reconstructed extinct-species target sequence."
        )
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run with synthetic demo genomes (no input files needed).",
    )
    parser.add_argument(
        "--proxy",
        type=str,
        default=None,
        help="Path to proxy species FASTA genome.",
    )
    parser.add_argument(
        "--target",
        type=str,
        default=None,
        help="Path to target (extinct) species FASTA genome.",
    )
    parser.add_argument(
        "--gff",
        type=str,
        default=None,
        help="Optional GFF3 annotation file for impact annotation.",
    )
    parser.add_argument(
        "--output-dir",
        dest="output_dir",
        type=str,
        default="output_calculator",
        help="Directory for output files (default: output_calculator).",
    )
    parser.add_argument(
        "--no-plot",
        dest="no_plot",
        action="store_true",
        help="Skip generating plots.",
    )
    parser.add_argument(
        "--plot-format",
        dest="plot_format",
        choices=["png", "svg"],
        default="png",
        help="Output format for plots (default: png).",
    )
    parser.add_argument(
        "--n-samples",
        dest="n_samples",
        type=int,
        default=1000,
        help="Number of 500 bp windows to sample for alignment (default: 1000).",
    )
    parser.add_argument(
        "--window-size",
        dest="window_size",
        type=int,
        default=500,
        help="Window size for pairwise alignment (default: 500).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42).",
    )
    return parser


def _validate(args, parser):
    if not args.demo and (args.proxy is None or args.target is None):
        parser.error("Provide --proxy and --target, or use --demo.")
    if args.proxy and not os.path.isfile(args.proxy):
        parser.error(f"Proxy FASTA not found: {args.proxy}")
    if args.target and not os.path.isfile(args.target):
        parser.error(f"Target FASTA not found: {args.target}")
    if args.gff and not os.path.isfile(args.gff):
        parser.error(f"GFF3 file not found: {args.gff}")
    if args.n_samples < 1:
        parser.error("--n-samples must be >= 1.")


def main() -> int:
    parser = _build_parser()
    args = parser.parse_args()
    _validate(args, parser)

    from src.pipeline import run_pipeline

    result = run_pipeline(args)
    b = result["burden"]

    print(f"Total variants called : {b['total_edits']}")
    print(f"Weighted edit burden  : {b['weighted_burden']}")
    print(f"Normalized burden     : {b['normalized_burden_per_mb']:.2f} per Mb")

    crispr_cycles = b.get(
        "estimated_crispr_cycles",
        -(-b["total_edits"] // 10),   # ceiling division fallback
    )
    print(f"Estimated CRISPR cycles (10 edits/cycle): {crispr_cycles}")

    print("\nBreakdown by variant class:")
    for vc, count in sorted(b["class_counts"].items()):
        print(f"  {vc:20s}: {count}")

    out = result["outputs"]
    print(f"\nOutputs written to: {args.output_dir}")
    print(f"  VCF             : {out['vcf']}")
    print(f"  Burden CSV      : {out['burden_csv']}")
    print(f"  Prioritized CSV : {out['prioritized_csv']}")
    print(f"  Burden JSON     : {out['burden_json']}")
    for p in out.get("plots", []):
        print(f"  Plot            : {p}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
