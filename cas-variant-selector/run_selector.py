#!/usr/bin/env python3
"""CLI for the Cas Variant Selector.

Usage:
    python run_selector.py locus.fa --goal knockout --output-dir results/
    python run_selector.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_selector.py",
        description=(
            "Cas Variant Selector — ranks Cas nucleases and editors for a target "
            "locus and editing goal based on PAM density, delivery constraints, "
            "and editing window compatibility."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("locus", metavar="LOCUS", nargs="?", default=None,
                        help="Target locus FASTA. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Use a synthetic 300 bp locus.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--goal", default="knockout",
                        choices=["knockout", "base-edit", "prime-edit", "activation", "repression"],
                        help="Desired editing goal.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.locus is None:
        _fatal("Provide a locus FASTA or use --demo.")
    if not args.demo and args.locus and not Path(args.locus).exists():
        _fatal(f"Locus FASTA not found: {args.locus}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            locus_fasta=Path(args.locus) if args.locus else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            editing_goal=args.goal,
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
    rankings = result.get("rankings", [])
    print(f"\nCas Variant Selector v{result['pipeline_version']}")
    print(f"Variants ranked: {len(rankings)}")
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
