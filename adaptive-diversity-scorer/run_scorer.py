#!/usr/bin/env python3
"""CLI for the Adaptive Diversity Scorer.

Usage:
    python run_scorer.py variants.vcf --environment env.tsv --output-dir results/
    python run_scorer.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_scorer.py",
        description=(
            "Adaptive Diversity Scorer — identifies putatively adaptive loci via Fst "
            "outlier detection and environment–genotype associations, then scores the "
            "population's adaptive diversity relative to neutral baseline."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="VCF of the population. Omit with --demo.")
    parser.add_argument("--environment", metavar="TSV", default=None,
                        help="TSV of environmental variables per sample.")
    parser.add_argument("--adaptive-loci", metavar="BED", default=None,
                        help="Pre-identified adaptive loci BED.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.vcf is None:
        _fatal("Provide a VCF or use --demo.")
    for p in [args.vcf, args.environment, args.adaptive_loci]:
        if p and not Path(p).exists():
            _fatal(f"File not found: {p}")


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
            environment_data=Path(args.environment) if args.environment else None,
            adaptive_loci_bed=Path(args.adaptive_loci) if args.adaptive_loci else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
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
    score = result.get("score")
    print(f"\nAdaptive Diversity Scorer v{result['pipeline_version']}")
    if score:
        print(f"Diversity class: {score.get('diversity_class', 'N/A')}")
    print(f"Adaptive loci identified: {len(result.get('adaptive_loci', []))}")
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
