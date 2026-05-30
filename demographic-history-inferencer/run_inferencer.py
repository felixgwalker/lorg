#!/usr/bin/env python3
"""CLI for the Demographic History Inferencer.

Usage:
    python run_inferencer.py variants.vcf --output-dir results/
    python run_inferencer.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_inferencer.py",
        description=(
            "Demographic History Inferencer — fits parametric demographic models "
            "(constant Ne, exponential growth, two-epoch, three-epoch) to the observed "
            "site frequency spectrum and selects the best model by AIC."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="VCF of SNPs for one population. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--models", nargs="*",
                        choices=["constant", "exponential_growth", "two_epoch",
                                 "three_epoch", "isolation_with_migration"],
                        default=None, help="Demographic models to fit (default: all).")
    parser.add_argument("--generation-time", type=float, default=30.0, metavar="YEARS",
                        help="Assumed generation time in years.")
    parser.add_argument("--mutation-rate", type=float, default=1.25e-8, metavar="FLOAT",
                        help="Per-base per-generation mutation rate.")
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
    if args.generation_time <= 0:
        _fatal("--generation-time must be positive.")
    if args.mutation_rate <= 0:
        _fatal("--mutation-rate must be positive.")


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
            models=args.models,
            generation_time=args.generation_time,
            mutation_rate=args.mutation_rate,
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
    best = result.get("best_model")
    all_fits = result.get("all_fits", [])
    print(f"\nDemographic History Inferencer v{result['pipeline_version']}")
    if best:
        print(f"Best model: {best.get('model', 'N/A')}  AIC: {best.get('aic', 'N/A'):.2f}")
    print(f"Models fitted: {len(all_fits)}")
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
