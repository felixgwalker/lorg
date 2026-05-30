#!/usr/bin/env python3
"""CLI for the Allele Frequency Comparator.

Usage:
    python run_comparator.py variants.vcf --af-table gnomad_af.tsv --output-dir results/
    python run_comparator.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_comparator.py",
        description=(
            "Allele Frequency Comparator — compares allele frequencies across gnomAD "
            "populations, flagging population-specific variants and computing pairwise "
            "fold changes and Fst estimates."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="VCF or variant list to query. Omit with --demo.")
    parser.add_argument("--af-table", metavar="TSV", default=None,
                        help="gnomAD-format population AF TSV.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--populations", nargs="*", metavar="POP",
                        default=None,
                        help="gnomAD population codes to compare (e.g. afr eas nfe).")
    parser.add_argument("--fold-change-threshold", type=float, default=5.0, metavar="FLOAT",
                        help="Minimum fold change to flag a differential variant.")
    parser.add_argument("--min-af", type=float, default=1e-5, metavar="FLOAT",
                        help="Minimum AF in at least one population to include a variant.")
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
    if not args.demo and args.af_table and not Path(args.af_table).exists():
        _fatal(f"AF table not found: {args.af_table}")


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
            af_table=Path(args.af_table) if args.af_table else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            populations=args.populations,
            fold_change_threshold=args.fold_change_threshold,
            min_af=args.min_af,
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
    comparisons = result.get("comparisons", [])
    print(f"\nAllele Frequency Comparator v{result['pipeline_version']}")
    print(f"Variants compared: {len(comparisons)}")
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
