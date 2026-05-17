#!/usr/bin/env python3
"""CLI for the Synteny Block Visualiser.

Usage:
    python run_visualiser.py --fasta1 genome1.fa --fasta2 genome2.fa --output-dir results/
    python run_visualiser.py --demo --output-dir results/
"""

import argparse
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_visualiser.py",
        description=(
            "Synteny Block Visualiser — identifies conserved synteny blocks between "
            "two genome assemblies and visualises chromosomal rearrangements."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--fasta1", metavar="FASTA", default=None,
                        help="First genome assembly FASTA.")
    parser.add_argument("--fasta2", metavar="FASTA", default=None,
                        help="Second genome assembly FASTA.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        dest="output_dir",
                        help="Output directory. Created if absent.")
    parser.add_argument("--kmer-size", type=int, default=15, dest="kmer_size",
                        help="K-mer size for seeding.")
    parser.add_argument("--min-chain-score", type=int, default=3, dest="min_chain_score",
                        help="Minimum number of seeds per synteny block.")
    parser.add_argument("--min-block-length", type=int, default=1000, dest="min_block_length",
                        help="Minimum synteny block length in bp.")
    parser.add_argument("--demo", action="store_true",
                        help="Run with synthetic demo genomes (no input files required).")
    parser.add_argument("--no-plot", action="store_true", dest="no_plot",
                        help="Skip plot generation.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.demo:
        return
    if not args.fasta1:
        print("ERROR: --fasta1 is required unless --demo is used.", file=sys.stderr)
        sys.exit(1)
    if not args.fasta2:
        print("ERROR: --fasta2 is required unless --demo is used.", file=sys.stderr)
        sys.exit(1)
    if not Path(args.fasta1).exists():
        print(f"ERROR: FASTA file not found: {args.fasta1}", file=sys.stderr)
        sys.exit(1)
    if not Path(args.fasta2).exists():
        print(f"ERROR: FASTA file not found: {args.fasta2}", file=sys.stderr)
        sys.exit(1)


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args)

    from src.pipeline import run_pipeline
    result = run_pipeline(args)
    return 0 if result.get("status") == "ok" else 1


if __name__ == "__main__":
    sys.exit(main())
