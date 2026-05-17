import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

from src.pipeline import run_pipeline


def build_parser():
    parser = argparse.ArgumentParser(
        description="Predict genomic locus tolerance to large DNA inserts.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bed", help="BED file of target loci")
    parser.add_argument("--fasta", help="Reference FASTA file")
    parser.add_argument("--gff3", help="Optional GFF3 annotation file")
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
    if not args.fasta:
        print("ERROR: --fasta is required (or use --demo)", file=sys.stderr)
        return False
    if not os.path.isfile(args.fasta):
        print(f"ERROR: FASTA file not found: {args.fasta}", file=sys.stderr)
        return False
    if args.gff3 and not os.path.isfile(args.gff3):
        print(f"ERROR: GFF3 file not found: {args.gff3}", file=sys.stderr)
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
