"""CLI entry point for erv-risk-mapper."""

import argparse
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_mapper",
        description="Identify ERV elements and risk-score them for reactivation potential.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--genome", type=Path, metavar="FASTA",
                   help="Input genome FASTA file")
    p.add_argument("--output-dir", type=Path, default=Path("results"),
                   metavar="PATH", help="Output directory")
    p.add_argument("--no-plot", action="store_true", help="Skip chromosome density plot")
    p.add_argument("--plot-format", choices=["png", "svg"], default="png",
                   help="Output plot format")
    p.add_argument("--window-size", type=int, default=5000,
                   help="Genome scanning window size (bp)")
    p.add_argument("--demo", action="store_true",
                   help="Run demo: synthetic 100kb genome (5 chroms x 20kb) with 8 planted ERVs")
    return p


def validate_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    if not args.demo:
        if args.genome is None:
            parser.error("--genome is required unless --demo is set")
        if not args.genome.exists():
            parser.error(f"Genome file not found: {args.genome}")
    if args.window_size < 1000:
        parser.error("--window-size must be at least 1000 bp")


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    validate_args(args, parser)

    from src.fasta_parser import parse_fasta, make_demo_genome
    from src.pipeline import run_pipeline

    if args.demo:
        genome = make_demo_genome()
        total_bp = sum(len(v) for v in genome.values())
        print(f"Demo mode: {len(genome)} chromosomes | {total_bp:,} bp | 8 planted ERV elements")
    else:
        genome = parse_fasta(args.genome)
        total_bp = sum(len(v) for v in genome.values())
        print(f"Genome: {len(genome)} sequences | {total_bp:,} bp")

    result = run_pipeline(
        genome=genome,
        output_dir=args.output_dir,
        no_plot=args.no_plot,
        plot_fmt=args.plot_format,
    )

    print(f"ERV elements detected: {result['n_erv_hits']}")
    print(f"  High risk:     {result['tier_high']}")
    print(f"  Moderate risk: {result['tier_moderate']}")
    print(f"  Low risk:      {result['tier_low']}")
    print(f"  BED file:      {result['bed']}")
    print(f"  CSV report:    {result['csv']}")
    if result["plot"]:
        print(f"  Density plot:  {result['plot']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
