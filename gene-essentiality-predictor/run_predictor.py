#!/usr/bin/env python3
"""CLI for the Gene Essentiality Predictor.

Usage:
    python run_predictor.py --gene-list genes.txt --depmap-scores depmap.tsv --output-dir results/
    python run_predictor.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_predictor.py",
        description=(
            "Gene Essentiality Predictor — aggregates DepMap CRISPR fitness scores, "
            "RNAi dependency data, and gnomAD constraint metrics (LOEUF, pLI) into a "
            "composite essentiality score with three-class classification."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--gene-list", metavar="FILE", default=None,
                        help="Text file of gene IDs or symbols to assess.")
    parser.add_argument("--depmap-scores", metavar="TSV", default=None,
                        help="DepMap CRISPR gene effect score TSV.")
    parser.add_argument("--constraint-data", metavar="TSV", default=None,
                        help="gnomAD constraint TSV (LOEUF, pLI).")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--cell-line", metavar="ID", default=None,
                        help="DepMap cell line ID for context-specific essentiality.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.gene_list is None:
        _fatal("Provide --gene-list or use --demo.")
    for p in [args.gene_list, args.depmap_scores, args.constraint_data]:
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
            gene_list=Path(args.gene_list) if args.gene_list else None,
            depmap_scores=Path(args.depmap_scores) if args.depmap_scores else None,
            constraint_data=Path(args.constraint_data) if args.constraint_data else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            cell_line_context=args.cell_line,
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
    preds = result.get("predictions", [])
    print(f"\nGene Essentiality Predictor v{result['pipeline_version']}")
    print(f"Genes assessed: {len(preds)}")
    print(f"Essential: {result.get('n_essential', 'N/A')}")
    print(f"Context-dependent: {result.get('n_context_dependent', 'N/A')}")
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
