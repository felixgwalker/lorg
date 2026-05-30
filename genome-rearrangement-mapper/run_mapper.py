#!/usr/bin/env python3
"""CLI for the Genome Rearrangement Mapper.

Usage:
    python run_mapper.py --synteny-blocks blocks.tsv --output-dir results/
    python run_mapper.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_mapper.py",
        description=(
            "Genome Rearrangement Mapper — identifies inversions, translocations, "
            "chromosome fusions, and fissions between two genomes from synteny block "
            "coordinates."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--synteny-blocks", metavar="TSV", default=None,
                        help="TSV of synteny blocks between species A and B.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--rearrangement-types", nargs="*",
                        choices=["inversion", "translocation", "fusion", "fission",
                                 "transposition"],
                        default=None, help="Rearrangement types to detect (default: all).")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.synteny_blocks is None:
        _fatal("Provide --synteny-blocks or use --demo.")
    if args.synteny_blocks and not Path(args.synteny_blocks).exists():
        _fatal(f"Synteny blocks file not found: {args.synteny_blocks}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            synteny_blocks=Path(args.synteny_blocks) if args.synteny_blocks else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            rearrangement_types=args.rearrangement_types,
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
    rearrangements = result.get("rearrangements", [])
    print(f"\nGenome Rearrangement Mapper v{result['pipeline_version']}")
    print(f"Rearrangements detected: {len(rearrangements)}")
    print(f"Breakpoints identified: {result.get('n_breakpoints', 'N/A')}")
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
