#!/usr/bin/env python3
"""CLI for the Base Edit Outcome Predictor.

Usage:
    python run_predictor.py targets.tsv --editor ABE8e --output-dir results/
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
            "Base Edit Outcome Predictor — predicts per-base editing probabilities "
            "for CBE and ABE editors using BE-Hive-derived position and "
            "trinucleotide context weights."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("targets", metavar="TARGETS", nargs="?", default=None,
                        help="TSV with id, spacer, target_sequence columns. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic spacer/target pairs.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--editor", default="CBE3",
                        choices=["CBE3", "BE4max", "ABE8e", "ABEmax", "CBE4-max"],
                        help="Base editor type.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip heatmap generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.targets is None:
        _fatal("Provide a targets TSV or use --demo.")
    if not args.demo and args.targets and not Path(args.targets).exists():
        _fatal(f"Targets file not found: {args.targets}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            targets_path=Path(args.targets) if args.targets else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            editor=args.editor,
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
    outcomes = result.get("outcomes", [])
    print(f"\nBase Edit Outcome Predictor v{result['pipeline_version']}")
    print(f"Targets predicted: {len(outcomes)}")
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
