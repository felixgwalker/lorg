#!/usr/bin/env python3
"""CLI for the Ancestral Gene Content Reconstructor.

Usage:
    python run_reconstructor.py --presence-table genes.tsv --phylogeny tree.nwk --output-dir results/
    python run_reconstructor.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_reconstructor.py",
        description=(
            "Ancestral Gene Content Reconstructor — reconstructs the set of genes "
            "present in ancestral genomes at internal phylogenetic nodes using "
            "Dollo parsimony or a Bayesian gain/loss model."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--presence-table", metavar="TSV", default=None,
                        help="TSV of gene presence/absence per species.")
    parser.add_argument("--phylogeny", metavar="NWK", default=None,
                        help="Newick species tree.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--method", choices=["parsimony", "bayesian"],
                        default="parsimony", help="Reconstruction method.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.presence_table is None:
        _fatal("Provide --presence-table or use --demo.")
    if args.presence_table and not Path(args.presence_table).exists():
        _fatal(f"Presence table not found: {args.presence_table}")
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
            presence_absence_table=Path(args.presence_table) if args.presence_table else None,
            phylogeny=Path(args.phylogeny) if args.phylogeny else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            method=args.method,
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
    genomes = result.get("ancestral_genomes", [])
    print(f"\nAncestral Gene Content Reconstructor v{result['pipeline_version']}")
    print(f"Genes assessed: {result.get('n_genes_assessed', 'N/A')}")
    print(f"Ancestral genomes reconstructed: {len(genomes)}")
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
