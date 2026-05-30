#!/usr/bin/env python3
"""CLI for the Expression Divergence Scorer.

Usage:
    python run_scorer.py --expression-a speciesA.tsv --expression-b speciesB.tsv --orthologs orthologs.tsv --output-dir results/
    python run_scorer.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_scorer.py",
        description=(
            "Expression Divergence Scorer — computes expression divergence for orthologous "
            "gene pairs between two species using log fold change, JSI, tau, or Euclidean "
            "distance metrics across matched tissues."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--expression-a", metavar="TSV", default=None,
                        help="TPM/count matrix for species A (genes × tissues).")
    parser.add_argument("--expression-b", metavar="TSV", default=None,
                        help="TPM/count matrix for species B (genes × tissues).")
    parser.add_argument("--orthologs", metavar="TSV", default=None,
                        help="TSV of ortholog gene pairs.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--metric",
                        choices=["log_fold_change", "tau", "JSI", "euclidean"],
                        default="log_fold_change",
                        help="Divergence metric.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.expression_a is None:
        _fatal("Provide --expression-a and --expression-b, or use --demo.")
    for p in [args.expression_a, args.expression_b, args.orthologs]:
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
            expression_a=Path(args.expression_a) if args.expression_a else None,
            expression_b=Path(args.expression_b) if args.expression_b else None,
            ortholog_table=Path(args.orthologs) if args.orthologs else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            metric=args.metric,
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
    scores = result.get("divergence_scores", [])
    print(f"\nExpression Divergence Scorer v{result['pipeline_version']}")
    print(f"Ortholog pairs scored: {len(scores)}")
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
