#!/usr/bin/env python3
"""CLI for the Genomic Diversity Index.

Usage:
    python run_index.py variants.vcf --output-dir results/
    python run_index.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_index.py",
        description=(
            "Genomic Diversity Index — computes θW, θπ, Tajima's D, Ho, He, and Fis "
            "in sliding windows across the genome and produces a genome-wide diversity "
            "summary for one population."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="VCF of SNPs. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--window-size", type=int, default=50000, metavar="INT",
                        help="Sliding window size in base pairs.")
    parser.add_argument("--step-size", type=int, default=10000, metavar="INT",
                        help="Step size in base pairs.")
    parser.add_argument("--metrics", nargs="*",
                        choices=["theta_w", "theta_pi", "tajimas_d", "Ho", "He", "Fis"],
                        default=None, help="Metrics to compute (default: all).")
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
    if args.step_size > args.window_size:
        _fatal("--step-size must be ≤ --window-size.")


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
            window_size=args.window_size,
            step_size=args.step_size,
            metrics=args.metrics,
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
    idx = result.get("population_index")
    print(f"\nGenomic Diversity Index v{result['pipeline_version']}")
    if idx:
        print(f"θW (genome-wide): {idx.get('genome_wide_theta_w', 'N/A'):.6f}")
    print(f"Windows computed: {len(result.get('windows', []))}")
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
