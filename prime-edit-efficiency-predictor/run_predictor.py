#!/usr/bin/env python3
"""CLI for the Prime Edit Efficiency Predictor.

Usage:
    python run_predictor.py designs.tsv --target target.fa --output-dir results/
    python run_predictor.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_predictor.py",
        description=(
            "Prime Edit Efficiency Predictor — scores pegRNA designs using "
            "DeepPrime-style features (PBS GC, RT length, MFE approximation, "
            "nick distance) and returns efficiency predictions."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("pegrna", metavar="PEGRNA", nargs="?", default=None,
                        help="TSV/JSON of pegRNA designs. Omit with --demo.")
    parser.add_argument("--target", metavar="FASTA", default=None,
                        help="Target locus FASTA for context features.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic designs.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.pegrna is None:
        _fatal("Provide a pegRNA TSV/JSON or use --demo.")
    if not args.demo and args.pegrna and not Path(args.pegrna).exists():
        _fatal(f"pegRNA file not found: {args.pegrna}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            pegrna_path=Path(args.pegrna) if args.pegrna else None,
            target_fasta=Path(args.target) if args.target else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
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
    preds = result.get("predictions", [])
    print(f"\nPrime Edit Efficiency Predictor v{result['pipeline_version']}")
    print(f"Predictions: {len(preds)}")
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
