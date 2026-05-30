#!/usr/bin/env python3
"""CLI for the Population Viability Genomics Estimator.

Usage:
    python run_estimator.py variants.vcf --census-size 500 --output-dir results/
    python run_estimator.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_estimator.py",
        description=(
            "Population Viability Genomics Estimator — integrates inbreeding, Ne, "
            "genetic load, and adaptive diversity to project extinction probability "
            "over 50/100/200 year horizons and estimate genomic MVP size."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="VCF of the population. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--census-size", type=int, default=None, metavar="INT",
                        help="Current census population size.")
    parser.add_argument("--generation-time", type=float, default=5.0, metavar="YEARS",
                        help="Generation time in years.")
    parser.add_argument("--time-horizons", nargs="*", type=int, metavar="YEARS",
                        default=None, help="Time horizons for viability projection (years).")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.vcf is None:
        _fatal("Provide a VCF or use --demo.")
    if args.vcf and not Path(args.vcf).exists():
        _fatal(f"VCF not found: {args.vcf}")


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
            output_dir=Path(args.output_dir),
            demo=args.demo,
            census_size=args.census_size,
            generation_time=args.generation_time,
            time_horizons=args.time_horizons,
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
    projs = result.get("projections", [])
    print(f"\nPopulation Viability Genomics Estimator v{result['pipeline_version']}")
    print(f"Genomic MVP estimate: {result.get('minimum_viable_population_genomic', 'N/A')}")
    print(f"Time horizons projected: {len(projs)}")
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
