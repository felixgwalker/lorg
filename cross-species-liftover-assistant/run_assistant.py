#!/usr/bin/env python3
"""CLI for the Cross-Species Liftover Assistant.

Usage:
    python run_assistant.py intervals.bed --chain human_to_mouse.chain --output-dir results/
    python run_assistant.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_assistant.py",
        description=(
            "Cross-Species Liftover Assistant — maps genomic BED intervals from a source "
            "species to a target species using a UCSC chain file, with synteny-based "
            "fallback for intervals lacking chain coverage."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("bed", metavar="BED", nargs="?", default=None,
                        help="BED of intervals to lift over. Omit with --demo.")
    parser.add_argument("--chain", metavar="CHAIN", default=None,
                        help="UCSC-format chain file (source→target).")
    parser.add_argument("--source-fasta", metavar="FASTA", default=None,
                        help="Source species FASTA (for BLAST fallback).")
    parser.add_argument("--target-fasta", metavar="FASTA", default=None,
                        help="Target species FASTA (for BLAST fallback).")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--min-identity", type=float, default=70.0, metavar="FLOAT",
                        help="Minimum alignment identity (%%) to accept a liftover.")
    parser.add_argument("--no-synteny-fallback", action="store_true",
                        help="Disable synteny-based fallback.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.bed is None:
        _fatal("Provide a BED file or use --demo.")
    if args.bed and not Path(args.bed).exists():
        _fatal(f"BED not found: {args.bed}")
    if args.chain and not Path(args.chain).exists():
        _fatal(f"Chain file not found: {args.chain}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            input_bed=Path(args.bed) if args.bed else None,
            chain_file=Path(args.chain) if args.chain else None,
            source_fasta=Path(args.source_fasta) if args.source_fasta else None,
            target_fasta=Path(args.target_fasta) if args.target_fasta else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            min_identity=args.min_identity,
            synteny_fallback=not args.no_synteny_fallback,
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
    print(f"\nCross-Species Liftover Assistant v{result['pipeline_version']}")
    print(f"Intervals input:   {result.get('n_input', 'N/A')}")
    print(f"Successfully mapped: {result.get('n_success', 'N/A')}")
    print(f"Failed:            {result.get('n_failed', 'N/A')}")
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
