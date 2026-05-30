#!/usr/bin/env python3
"""CLI for the Assembly Gap Analyser.

Usage:
    python run_analyser.py assembly.fa --annotation genes.gtf --output-dir results/
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
            "Assembly Gap Analyser — identifies and classifies all N-runs in a genome "
            "assembly, determines whether gaps interrupt gene models or fall between "
            "synteny-supported gene pairs."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("assembly", metavar="FASTA", nargs="?", default=None,
                        help="Assembly FASTA. Omit with --demo.")
    parser.add_argument("--annotation", metavar="GTF", default=None,
                        help="Gene annotation GTF for gap-gene overlap analysis.")
    parser.add_argument("--reference-annotation", metavar="GTF", default=None,
                        help="Reference species annotation for synteny-gap analysis.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--min-gap-length", type=int, default=10, metavar="INT",
                        help="Minimum N run length to report.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.assembly is None:
        _fatal("Provide an assembly FASTA or use --demo.")
    for p in [args.assembly, args.annotation, args.reference_annotation]:
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
            assembly_fasta=Path(args.assembly) if args.assembly else None,
            gene_annotation=Path(args.annotation) if args.annotation else None,
            reference_annotation=Path(args.reference_annotation) if args.reference_annotation else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            min_gap_length=args.min_gap_length,
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
    gaps = result.get("gaps", [])
    s = result.get("summary")
    print(f"\nAssembly Gap Analyser v{result['pipeline_version']}")
    print(f"Gaps identified: {len(gaps)}")
    if s:
        print(f"Total gap length: {s.get('total_gap_length_bp', 'N/A'):,} bp")
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
