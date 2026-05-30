#!/usr/bin/env python3
"""CLI for the Chromatin Accessibility Scorer.

Usage:
    python run_scorer.py sample.bam --peaks peaks.bed --output-dir results/
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
            "Chromatin Accessibility Scorer — counts ATAC-seq fragment insertions within "
            "peak regions, normalises by library size, and classifies regions as open, "
            "intermediate, or closed."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("bam", metavar="BAM", nargs="?", default=None,
                        help="ATAC-seq BAM file. Omit with --demo.")
    parser.add_argument("--peaks", metavar="BED", default=None,
                        help="BED of peak regions to score.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--normalisation", choices=["RPM", "RPKM", "TMM"],
                        default="RPM", help="Signal normalisation method.")
    parser.add_argument("--q-value-threshold", type=float, default=0.05, metavar="FLOAT",
                        help="FDR threshold for open chromatin classification.")
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
    if not args.demo and args.peaks and not Path(args.peaks).exists():
        _fatal(f"Peaks BED not found: {args.peaks}")


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
            peaks_bed=Path(args.peaks) if args.peaks else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            normalisation=args.normalisation,
            q_value_threshold=args.q_value_threshold,
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
    print(f"\nChromatin Accessibility Scorer v{result['pipeline_version']}")
    print(f"Open regions: {result.get('n_open_regions', 'N/A')}")
    print(f"Regions scored: {len(result.get('scores', []))}")
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
