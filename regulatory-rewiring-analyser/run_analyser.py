#!/usr/bin/env python3
"""CLI for the Regulatory Rewiring Analyser.

Usage:
    python run_analyser.py --elements-a speciesA.bed --elements-b speciesB.bed --orthologs orthologs.tsv --output-dir results/
    python run_analyser.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_analyser.py",
        description=(
            "Regulatory Rewiring Analyser — compares regulatory elements (enhancers, "
            "promoters) near orthologous genes between two species to classify each "
            "element as gained, lost, conserved, or relocated."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--elements-a", metavar="BED", default=None,
                        help="BED of regulatory elements in species A.")
    parser.add_argument("--elements-b", metavar="BED", default=None,
                        help="BED of regulatory elements in species B.")
    parser.add_argument("--orthologs", metavar="TSV", default=None,
                        help="TSV of ortholog gene pairs.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--min-conservation", type=float, default=0.7, metavar="FLOAT",
                        help="Minimum sequence identity to call an element conserved.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.elements_a is None:
        _fatal("Provide --elements-a and --elements-b or use --demo.")
    for p in [args.elements_a, args.elements_b, args.orthologs]:
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
            elements_a=Path(args.elements_a) if args.elements_a else None,
            elements_b=Path(args.elements_b) if args.elements_b else None,
            ortholog_table=Path(args.orthologs) if args.orthologs else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            min_conservation=args.min_conservation,
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
    rewirings = result.get("rewirings", [])
    print(f"\nRegulatory Rewiring Analyser v{result['pipeline_version']}")
    print(f"Elements assessed: {result.get('n_elements_assessed', 'N/A')}")
    print(f"Rewiring events identified: {len(rewirings)}")
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
