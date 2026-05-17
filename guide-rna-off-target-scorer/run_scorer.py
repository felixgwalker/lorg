"""CLI entry point for guide-rna-off-target-scorer."""

import argparse
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_scorer",
        description="Predict CRISPR off-target cleavage sites using CFD scoring.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--guides", type=Path, metavar="FILE",
                   help="Guide RNA sequences (FASTA or CSV, 20nt protospacers)")
    p.add_argument("--genome", type=Path, metavar="FASTA",
                   help="Genome FASTA to search for off-targets")
    p.add_argument("--output-dir", type=Path, default=Path("results"),
                   metavar="PATH", help="Directory for output files")
    p.add_argument("--no-plot", action="store_true", help="Skip plot generation")
    p.add_argument("--plot-format", choices=["png", "svg"], default="png",
                   help="Plot file format")
    p.add_argument("--max-mismatches", type=int, default=3, choices=[1, 2, 3],
                   help="Maximum mismatches to consider as off-target")
    p.add_argument("--demo", action="store_true",
                   help="Run demo with synthetic 5000bp genome and 3 guides")
    return p


def validate_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    if not args.demo:
        if args.guides is None:
            parser.error("--guides is required unless --demo is set")
        if args.genome is None:
            parser.error("--genome is required unless --demo is set")
        if not args.guides.exists():
            parser.error(f"Guide file not found: {args.guides}")
        if not args.genome.exists():
            parser.error(f"Genome file not found: {args.genome}")


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    validate_args(args, parser)

    from src.guide_parser import parse_guides, validate_guide, make_demo_guides
    from src.genome_indexer import load_fasta_genome, make_demo_genome
    from src.pipeline import run_pipeline

    if args.demo:
        guides = make_demo_guides()
        genome = make_demo_genome()
        print("Demo mode: synthetic 5000bp genome, 3 guide RNAs")
    else:
        raw_guides = parse_guides(args.guides)
        guides = [(name, validate_guide(seq)) for name, seq in raw_guides]
        genome = load_fasta_genome(args.genome)
        total_bp = sum(len(v) for v in genome.values())
        print(f"Loaded {len(guides)} guides | genome {total_bp:,} bp")

    result = run_pipeline(
        guides=guides,
        genome=genome,
        output_dir=args.output_dir,
        no_plot=args.no_plot,
        plot_fmt=args.plot_format,
    )

    print(f"Guides: {result['n_guides']}  |  Hits (<={args.max_mismatches} mm): {result['n_hits']}")
    print(f"  Off-target table:    {result['off_target_table']}")
    print(f"  Specificity summary: {result['specificity_summary']}")
    if result["plot"]:
        print(f"  Manhattan plot:      {result['plot']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
