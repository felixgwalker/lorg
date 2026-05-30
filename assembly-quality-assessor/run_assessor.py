#!/usr/bin/env python3
"""CLI for the Assembly Quality Assessor.

Usage:
    python run_assessor.py assembly.fa --output-dir results/
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
            "Assembly Quality Assessor — computes N50, N90, L50, L90, GC content, "
            "gap statistics, and ambiguous base count from a genome assembly FASTA, "
            "classifying the assembly quality level."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("assembly", metavar="FASTA", nargs="?", default=None,
                        help="Assembly FASTA file. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.assembly is None:
        _fatal("Provide an assembly FASTA or use --demo.")
    if args.assembly and not Path(args.assembly).exists():
        _fatal(f"Assembly FASTA not found: {args.assembly}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            assembly_fasta=Path(args.assembly) if args.assembly else None,
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
    stats = result.get("stats")
    print(f"\nAssembly Quality Assessor v{result['pipeline_version']}")
    if stats:
        print(f"Quality class: {stats.get('quality_class', 'N/A')}")
        print(f"N50: {stats.get('n50', 'N/A'):,} bp")
        print(f"Total length: {stats.get('total_length_bp', 'N/A'):,} bp")
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
