#!/usr/bin/env python3
"""CLI for the Inbreeding Risk Forecaster.

Usage:
    python run_forecaster.py variants.vcf --output-dir results/
    python run_forecaster.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_forecaster.py",
        description=(
            "Inbreeding Risk Forecaster — estimates current inbreeding from runs of "
            "homozygosity and Fis, infers Ne, and projects inbreeding accumulation over "
            "future generations to forecast conservation risk."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="VCF of the population to assess. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--min-roh-kb", type=int, default=500, metavar="INT",
                        help="Minimum ROH length in kilobases for FROH calculation.")
    parser.add_argument("--generation-time", type=float, default=5.0, metavar="YEARS",
                        help="Generation time in years.")
    parser.add_argument("--n-generations", type=int, default=50, metavar="INT",
                        help="Number of generations to project inbreeding.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.vcf is None:
        _fatal("Provide a VCF file or use --demo.")
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
            min_roh_kb=args.min_roh_kb,
            generation_time=args.generation_time,
            n_generations_forecast=args.n_generations,
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
    fc = result.get("forecast")
    print(f"\nInbreeding Risk Forecaster v{result['pipeline_version']}")
    if fc:
        print(f"Risk level: {fc.get('risk_level', 'N/A')}")
        print(f"Current F: {fc.get('current_f', 'N/A'):.3f}")
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
