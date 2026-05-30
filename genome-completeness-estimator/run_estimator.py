#!/usr/bin/env python3
"""CLI for the Genome Completeness Estimator.

Usage:
    python run_estimator.py assembly.fa --lineage vertebrata_odb10 --output-dir results/
    python run_estimator.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_estimator.py",
        description=(
            "Genome Completeness Estimator — searches BUSCO conserved orthologs from "
            "a specified lineage database in the assembly, reporting complete, "
            "fragmented, and missing gene fractions."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("assembly", metavar="FASTA", nargs="?", default=None,
                        help="Assembly FASTA. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--lineage", default="vertebrata_odb10", metavar="DATASET",
                        help="BUSCO lineage dataset.")
    parser.add_argument("--mode", choices=["genome", "proteins", "transcriptome"],
                        default="genome", help="BUSCO analysis mode.")
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
            lineage=args.lineage,
            mode=args.mode,
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
    print(f"\nGenome Completeness Estimator v{result['pipeline_version']}")
    if stats:
        print(f"Complete (S+D): {stats.get('complete_fraction', 0):.1%}")
        print(f"Fragmented: {stats.get('fragmented_fraction', 0):.1%}")
        print(f"Missing: {stats.get('missing_fraction', 0):.1%}")
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
