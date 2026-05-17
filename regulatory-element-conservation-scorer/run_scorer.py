import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

from src.pipeline import run_pipeline


def build_parser():
    parser = argparse.ArgumentParser(
        description="Score conservation of non-coding regulatory elements across species.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bed", help="BED file of regulatory elements")
    parser.add_argument("--fastas", nargs="+",
                        help="FASTA files, one per species (first is treated as reference)")
    parser.add_argument("--output-dir", default="output", dest="output_dir",
                        help="Directory for output files")
    parser.add_argument("--demo", action="store_true",
                        help="Run demo with synthetic data (no input files needed)")
    parser.add_argument("--no-plot", action="store_true", dest="no_plot",
                        help="Skip plot generation")
    return parser


def validate_args(args):
    if args.demo:
        return True
    if not args.bed:
        print("ERROR: --bed is required (or use --demo)", file=sys.stderr)
        return False
    if not os.path.isfile(args.bed):
        print(f"ERROR: BED file not found: {args.bed}", file=sys.stderr)
        return False
    if not args.fastas or len(args.fastas) < 2:
        print("ERROR: --fastas requires at least 2 FASTA files (or use --demo)", file=sys.stderr)
        return False
    for p in args.fastas:
        if not os.path.isfile(p):
            print(f"ERROR: FASTA file not found: {p}", file=sys.stderr)
            return False
    return True


def main():
    parser = build_parser()
    args = parser.parse_args()
    if not validate_args(args):
        sys.exit(1)
    result = run_pipeline(args)
    return 0 if result.get("status") == "ok" else 1


if __name__ == "__main__":
    sys.exit(main())
