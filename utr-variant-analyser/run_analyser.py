#!/usr/bin/env python3
"""CLI for the UTR Variant Analyser.

Usage:
    python run_analyser.py variants.vcf --fasta genome.fa --annotation utrs.bed --output-dir results/
    python run_analyser.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_analyser.py",
        description=(
            "UTR Variant Analyser — identifies uORF creation/disruption events, "
            "Kozak context changes, and polyadenylation signal disruptions caused by "
            "variants in 5' and 3' UTR regions."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="VCF of UTR variants. Omit with --demo.")
    parser.add_argument("--fasta", metavar="FASTA", default=None,
                        help="Reference genome FASTA.")
    parser.add_argument("--annotation", metavar="BED/GTF", default=None,
                        help="UTR annotation BED or GTF.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--utr-type", choices=["5prime", "3prime", "both"], default="both",
                        help="Which UTR regions to analyse.")
    parser.add_argument("--kozak-threshold", type=float, default=0.5, metavar="FLOAT",
                        help="Kozak score delta to flag as significant.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.vcf is None:
        _fatal("Provide a VCF file or use --demo.")
    if not args.demo and args.vcf and not Path(args.vcf).exists():
        _fatal(f"VCF not found: {args.vcf}")
    if not args.demo and args.fasta and not Path(args.fasta).exists():
        _fatal(f"FASTA not found: {args.fasta}")
    if not args.demo and args.annotation and not Path(args.annotation).exists():
        _fatal(f"Annotation file not found: {args.annotation}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            vcf_path=Path(args.vcf) if args.vcf else None,
            fasta_path=Path(args.fasta) if args.fasta else None,
            annotation_path=Path(args.annotation) if args.annotation else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            utr_type=args.utr_type,
            kozak_threshold=args.kozak_threshold,
            no_plot=args.no_plot,
        )
    except NotImplementedError:
        print("ERROR: This tool is not yet implemented.", file=sys.stderr)
        return 1
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        logging.exception("Unhandled exception")
        return 1
    _print_summary(result)
    return 0


def _print_summary(result: dict) -> None:
    results = result.get("results", [])
    print(f"\nUTR Variant Analyser v{result['pipeline_version']}")
    print(f"Variants analysed: {len(results)}")
    out = result.get("output_files", {})
    print("Output files:")
    for k, v in out.items():
        print(f"  {k:20s}: {v}")
    print()


def _fatal(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
