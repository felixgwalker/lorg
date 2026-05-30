#!/usr/bin/env python3
"""CLI for the Population Differentiation Scorer.

Usage:
    python run_scorer.py variants.vcf --pop-map populations.tsv --output-dir results/
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
            "Population Differentiation Scorer — computes Weir-Cockerham Fst, Gst, "
            "and Jost's D between all population pairs genome-wide and in sliding windows, "
            "highlighting outlier regions consistent with divergent selection."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="Multi-population VCF. Omit with --demo.")
    parser.add_argument("--pop-map", metavar="TSV", default=None,
                        help="TSV mapping sample IDs to population labels.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--metrics", nargs="*",
                        choices=["Fst", "Gst", "Jost_D", "Phi_st"],
                        default=None, help="Differentiation metrics to compute (default: all).")
    parser.add_argument("--window-size", type=int, default=50000, metavar="INT",
                        help="Sliding window size in base pairs.")
    parser.add_argument("--outlier-percentile", type=float, default=99.0, metavar="FLOAT",
                        help="Percentile threshold to flag Fst outlier windows.")
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
    if not args.demo and args.pop_map and not Path(args.pop_map).exists():
        _fatal(f"Population map not found: {args.pop_map}")


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
            population_map=Path(args.pop_map) if args.pop_map else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            metrics=args.metrics,
            window_size=args.window_size,
            outlier_percentile=args.outlier_percentile,
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
    gw = result.get("genome_wide", [])
    windows = result.get("windows", [])
    print(f"\nPopulation Differentiation Scorer v{result['pipeline_version']}")
    print(f"Population pairs scored: {len(gw)}")
    print(f"Windows computed: {len(windows)}")
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
