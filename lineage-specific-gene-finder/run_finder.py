#!/usr/bin/env python3
"""CLI for the Lineage-Specific Gene Finder.

Usage:
    python run_finder.py query.fa --outgroups out1.fa out2.fa --output-dir results/
    python run_finder.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_finder.py",
        description=(
            "Lineage-Specific Gene Finder — identifies orphan / taxonomically restricted "
            "genes in a query species by detecting proteins with no detectable homolog in "
            "a set of outgroup proteomes."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("query", metavar="FASTA", nargs="?", default=None,
                        help="Protein FASTA of the focal species. Omit with --demo.")
    parser.add_argument("--outgroups", nargs="*", metavar="FASTA", default=None,
                        help="Outgroup species protein FASTAs.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--min-outgroup-species", type=int, default=3, metavar="INT",
                        help="Minimum outgroup species that must lack a hit.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.query is None:
        _fatal("Provide a query FASTA or use --demo.")
    if args.query and not Path(args.query).exists():
        _fatal(f"Query FASTA not found: {args.query}")
    for t in (args.outgroups or []):
        if not Path(t).exists():
            _fatal(f"Outgroup FASTA not found: {t}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            query_proteome=Path(args.query) if args.query else None,
            outgroup_proteomes=[Path(t) for t in args.outgroups] if args.outgroups else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            min_outgroup_species=args.min_outgroup_species,
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
    genes = result.get("lineage_specific_genes", [])
    print(f"\nLineage-Specific Gene Finder v{result['pipeline_version']}")
    print(f"Genes screened: {result.get('n_genes_screened', 'N/A')}")
    print(f"Lineage-specific genes found: {len(genes)}")
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
