#!/usr/bin/env python3
"""CLI for the Off-Target Cluster Detector.

Usage:
    python run_detector.py offtargets.bed --output-dir results/
    python run_detector.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_detector.py",
        description=(
            "Off-Target Cluster Detector — identifies genomic hotspots of CRISPR "
            "off-target sites using sliding-window density scoring and "
            "DBSCAN-style single-linkage clustering."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("offtargets", metavar="BED", nargs="?", default=None,
                        help="BED/TSV of off-target sites with score column. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Use synthetic off-target site list.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--window-size", type=int, default=100_000, metavar="INT",
                        help="Sliding window size in bp for density scan.")
    parser.add_argument("--min-cluster-size", type=int, default=3, metavar="INT",
                        help="Minimum off-target sites to call a cluster.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip Manhattan plot generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.offtargets is None:
        _fatal("Provide an off-target BED/TSV or use --demo.")
    if not args.demo and args.offtargets and not Path(args.offtargets).exists():
        _fatal(f"Off-targets file not found: {args.offtargets}")
    if args.min_cluster_size < 2:
        _fatal("--min-cluster-size must be >= 2.")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            offtargets_path=Path(args.offtargets) if args.offtargets else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            window_size=args.window_size,
            min_cluster_size=args.min_cluster_size,
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
    clusters = result.get("clusters", [])
    print(f"\nOff-Target Cluster Detector v{result['pipeline_version']}")
    print(f"Clusters detected: {len(clusters)}")
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
