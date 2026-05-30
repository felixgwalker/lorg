#!/usr/bin/env python3
"""CLI for the PAM Flexibility Predictor.

Usage:
    python run_predictor.py target.fa --output-dir results/
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
            "PAM Flexibility Predictor — scores PAM site availability for a "
            "panel of Cas variants at a target locus using IUPAC position "
            "weight matrices."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("fasta", metavar="FASTA", nargs="?", default=None,
                        help="Target locus FASTA. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Use a synthetic 500 bp target.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--cas-variants", nargs="+", default=None, metavar="VARIANT",
                        help="Specific Cas variants to score (default: all).")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.fasta is None:
        _fatal("Provide a FASTA file or use --demo.")
    if not args.demo and args.fasta and not Path(args.fasta).exists():
        _fatal(f"FASTA not found: {args.fasta}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            fasta_path=Path(args.fasta) if args.fasta else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            cas_variants=args.cas_variants,
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
    print(f"\nPAM Flexibility Predictor v{result['pipeline_version']}")
    scores = result.get("pam_scores", {})
    for variant, data in scores.items():
        print(f"  {variant}: {data}")
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
