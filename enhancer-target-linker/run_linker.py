#!/usr/bin/env python3
"""CLI for the Enhancer Target Linker.

Usage:
    python run_linker.py --enhancers enhancers.bed --gene-tss tss.bed --output-dir results/
    python run_linker.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_linker.py",
        description=(
            "Enhancer Target Linker — links enhancer elements to candidate target genes "
            "using the activity-by-contact (ABC) model, expression correlation, or "
            "distance-based scoring within a configurable window."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--enhancers", metavar="BED", default=None,
                        help="BED of enhancer elements with activity scores.")
    parser.add_argument("--gene-tss", metavar="BED", default=None,
                        help="BED of gene TSS positions.")
    parser.add_argument("--activity-matrix", metavar="TSV", default=None,
                        help="Enhancer activity matrix across samples (for correlation).")
    parser.add_argument("--hic-matrix", metavar="COOL/HIC", default=None,
                        help="Hi-C contact matrix (for ABC model).")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--method",
                        choices=["activity_by_contact", "correlation", "distance", "hi_c"],
                        default="activity_by_contact",
                        help="Enhancer-gene linking method.")
    parser.add_argument("--max-distance-bp", type=int, default=500000, metavar="INT",
                        help="Maximum enhancer-TSS distance to consider.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.enhancers is None:
        _fatal("Provide --enhancers and --gene-tss, or use --demo.")
    for p in [args.enhancers, args.gene_tss, args.activity_matrix, args.hic_matrix]:
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
            enhancers=Path(args.enhancers) if args.enhancers else None,
            gene_tss=Path(args.gene_tss) if args.gene_tss else None,
            activity_matrix=Path(args.activity_matrix) if args.activity_matrix else None,
            hic_matrix=Path(args.hic_matrix) if args.hic_matrix else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            method=args.method,
            max_distance_bp=args.max_distance_bp,
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
    links = result.get("links", [])
    print(f"\nEnhancer Target Linker v{result['pipeline_version']}")
    print(f"Enhancers: {result.get('n_enhancers', 'N/A')}")
    print(f"Enhancer-gene links: {len(links)}")
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
