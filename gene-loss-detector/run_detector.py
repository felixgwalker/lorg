#!/usr/bin/env python3
"""CLI for the Gene Loss Detector.

Usage:
    python run_detector.py --ortholog-table orthologs.tsv --phylogeny tree.nwk --output-dir results/
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
            "Gene Loss Detector — identifies genes present in most species but absent "
            "or pseudogenised in specific lineages, mapping the inferred branch of loss "
            "onto a species phylogeny."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--ortholog-table", metavar="TSV", default=None,
                        help="TSV of ortholog presence/absence per species.")
    parser.add_argument("--phylogeny", metavar="NWK", default=None,
                        help="Newick species tree.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--min-species-with-gene", type=int, default=3, metavar="INT",
                        help="Minimum species that must have a gene to assess for loss.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.ortholog_table is None:
        _fatal("Provide --ortholog-table or use --demo.")
    if args.ortholog_table and not Path(args.ortholog_table).exists():
        _fatal(f"Ortholog table not found: {args.ortholog_table}")
    if args.phylogeny and not Path(args.phylogeny).exists():
        _fatal(f"Phylogeny file not found: {args.phylogeny}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            ortholog_table=Path(args.ortholog_table) if args.ortholog_table else None,
            species_proteomes=None,
            phylogeny=Path(args.phylogeny) if args.phylogeny else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            min_species_with_gene=args.min_species_with_gene,
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
    losses = result.get("gene_losses", [])
    print(f"\nGene Loss Detector v{result['pipeline_version']}")
    print(f"Genes assessed: {result.get('n_genes_assessed', 'N/A')}")
    print(f"Gene losses detected: {len(losses)}")
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
