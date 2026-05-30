#!/usr/bin/env python3
"""CLI for the Palaeogenomic Coverage Assessor.

Usage:
    python run_assessor.py sample.bam --output-dir results/
    python run_assessor.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_assessor.py",
        description=(
            "Palaeogenomic Coverage Assessor — computes mapping rate, mean depth, "
            "breadth of coverage at multiple thresholds, endogenous fraction, and "
            "duplication rate for an ancient DNA BAM."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("bam", metavar="BAM", nargs="?", default=None,
                        help="BAM file of aligned ancient DNA reads. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--min-mapq", type=int, default=25, metavar="INT",
                        help="Minimum mapping quality for depth calculation.")
    parser.add_argument("--keep-duplicates", action="store_true",
                        help="Include duplicate reads in coverage calculation.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.bam is None:
        _fatal("Provide a BAM file or use --demo.")
    if args.bam and not Path(args.bam).exists():
        _fatal(f"BAM not found: {args.bam}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            bam_path=Path(args.bam) if args.bam else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            min_mapq=args.min_mapq,
            remove_duplicates=not args.keep_duplicates,
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
    a = result.get("assessment")
    print(f"\nPalaeogenomic Coverage Assessor v{result['pipeline_version']}")
    if a:
        print(f"Coverage class: {a.get('coverage_class', 'N/A')}")
        print(f"Mean depth: {a.get('mean_depth', 'N/A'):.3f}×")
        print(f"1× breadth: {a.get('fraction_covered_1x', 'N/A'):.1%}")
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
