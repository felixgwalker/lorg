#!/usr/bin/env python3
"""CLI for the Ortholog Mapper.

Usage:
    python run_mapper.py query.fa --targets target1.fa target2.fa --output-dir results/
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
            "Ortholog Mapper — maps orthologous genes between a query and target species "
            "using reciprocal best BLAST hits or OMA, classifying 1:1, 1:N, N:1 and N:N "
            "relationships with optional synteny validation."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("query", metavar="FASTA", nargs="?", default=None,
                        help="Query species protein FASTA. Omit with --demo.")
    parser.add_argument("--targets", nargs="*", metavar="FASTA", default=None,
                        help="Target species protein FASTAs.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--method",
                        choices=["reciprocal_best_hits", "OMA", "inparanoid"],
                        default="reciprocal_best_hits",
                        help="Ortholog detection method.")
    parser.add_argument("--min-identity", type=float, default=30.0, metavar="FLOAT",
                        help="Minimum sequence identity (%%) to accept a BLAST hit.")
    parser.add_argument("--synteny-support", action="store_true",
                        help="Require synteny evidence to confirm orthologs.")
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
    for t in (args.targets or []):
        if not Path(t).exists():
            _fatal(f"Target FASTA not found: {t}")


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
            target_proteomes=[Path(t) for t in args.targets] if args.targets else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            method=args.method,
            min_identity=args.min_identity,
            synteny_support=args.synteny_support,
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
    groups = result.get("ortholog_groups", [])
    print(f"\nOrtholog Mapper v{result['pipeline_version']}")
    print(f"Query genes: {result.get('n_query_genes', 'N/A')}")
    print(f"Ortholog groups: {len(groups)}")
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
