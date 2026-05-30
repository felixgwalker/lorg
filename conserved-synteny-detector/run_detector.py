#!/usr/bin/env python3
"""CLI for the Conserved Synteny Detector.

Usage:
    python run_detector.py --orthologs orthologs.tsv --positions-a speciesA.bed --positions-b speciesB.bed --output-dir results/
    python run_detector.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_detector.py",
        description=(
            "Conserved Synteny Detector — identifies synteny blocks between two genomes "
            "by chaining ortholog anchor pairs with a collinearity algorithm and reporting "
            "conserved regions with orientation."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--orthologs", metavar="TSV", default=None,
                        help="TSV of ortholog gene pairs between species A and B.")
    parser.add_argument("--positions-a", metavar="BED", default=None,
                        help="BED of gene positions in species A.")
    parser.add_argument("--positions-b", metavar="BED", default=None,
                        help="BED of gene positions in species B.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--min-anchors", type=int, default=5, metavar="INT",
                        help="Minimum anchor genes to form a synteny block.")
    parser.add_argument("--min-block-length-kb", type=int, default=100, metavar="INT",
                        help="Minimum block length in kilobases.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.orthologs is None:
        _fatal("Provide --orthologs or use --demo.")
    for p in [args.orthologs, args.positions_a, args.positions_b]:
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
            ortholog_table=Path(args.orthologs) if args.orthologs else None,
            gene_positions_a=Path(args.positions_a) if args.positions_a else None,
            gene_positions_b=Path(args.positions_b) if args.positions_b else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            min_anchors=args.min_anchors,
            min_block_length_kb=args.min_block_length_kb,
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
    blocks = result.get("synteny_blocks", [])
    print(f"\nConserved Synteny Detector v{result['pipeline_version']}")
    print(f"Synteny blocks detected: {len(blocks)}")
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
