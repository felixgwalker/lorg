#!/usr/bin/env python3
"""CLI for the Gene Family Expansion Detector.

Usage:
    python run_detector.py --family-table families.tsv --phylogeny tree.nwk --output-dir results/
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
            "Gene Family Expansion Detector — identifies gene families that have expanded "
            "significantly on specific branches using a birth-death model, comparing observed "
            "family sizes against a genome-wide rate."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--family-table", metavar="TSV", default=None,
                        help="TSV of gene family sizes per species.")
    parser.add_argument("--phylogeny", metavar="NWK", default=None,
                        help="Newick species tree with branch lengths.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--min-fold-expansion", type=float, default=2.0, metavar="FLOAT",
                        help="Minimum fold change to report as expanded.")
    parser.add_argument("--p-value-threshold", type=float, default=0.05, metavar="FLOAT",
                        help="FDR-corrected p-value threshold for significance.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.family_table is None:
        _fatal("Provide --family-table or use --demo.")
    if args.family_table and not Path(args.family_table).exists():
        _fatal(f"Family table not found: {args.family_table}")
    if args.phylogeny and not Path(args.phylogeny).exists():
        _fatal(f"Phylogeny not found: {args.phylogeny}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            family_table=Path(args.family_table) if args.family_table else None,
            phylogeny=Path(args.phylogeny) if args.phylogeny else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            min_fold_expansion=args.min_fold_expansion,
            p_value_threshold=args.p_value_threshold,
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
    expansions = result.get("expansions", [])
    print(f"\nGene Family Expansion Detector v{result['pipeline_version']}")
    print(f"Families assessed: {result.get('n_families_assessed', 'N/A')}")
    print(f"Significant expansions: {len(expansions)}")
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
