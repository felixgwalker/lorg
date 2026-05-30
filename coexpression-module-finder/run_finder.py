#!/usr/bin/env python3
"""CLI for the Coexpression Module Finder.

Usage:
    python run_finder.py --expression expression.tsv --output-dir results/
    python run_finder.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_finder.py",
        description=(
            "Coexpression Module Finder — identifies groups of co-expressed genes using "
            "WGCNA, clique-based clustering, or k-means, and correlates module eigengenes "
            "with sample traits."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--expression", metavar="TSV", default=None,
                        help="Normalised expression matrix (genes × samples).")
    parser.add_argument("--trait-data", metavar="TSV", default=None,
                        help="Sample trait data TSV (samples × traits).")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--method", choices=["WGCNA", "clique_based", "kmeans"],
                        default="WGCNA", help="Module detection method.")
    parser.add_argument("--min-module-size", type=int, default=30, metavar="INT",
                        help="Minimum number of genes per module.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.expression is None:
        _fatal("Provide --expression or use --demo.")
    for p in [args.expression, args.trait_data]:
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
            expression_matrix=Path(args.expression) if args.expression else None,
            trait_data=Path(args.trait_data) if args.trait_data else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            method=args.method,
            min_module_size=args.min_module_size,
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
    modules = result.get("modules", [])
    print(f"\nCoexpression Module Finder v{result['pipeline_version']}")
    print(f"Modules detected: {len(modules)}")
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
