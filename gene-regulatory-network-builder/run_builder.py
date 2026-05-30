#!/usr/bin/env python3
"""CLI for the Gene Regulatory Network Builder.

Usage:
    python run_builder.py --expression expression.tsv --tf-list tfs.txt --output-dir results/
    python run_builder.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_builder.py",
        description=(
            "Gene Regulatory Network Builder — infers regulatory edges between "
            "transcription factors and target genes from expression data using "
            "GENIE3, ARACNE, or Pearson correlation."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--expression", metavar="TSV", default=None,
                        help="Gene expression matrix (genes × samples).")
    parser.add_argument("--tf-list", metavar="FILE", default=None,
                        help="Text file of transcription factor gene IDs.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--method", choices=["GENIE3", "correlation", "ARACNE"],
                        default="GENIE3", help="Network inference method.")
    parser.add_argument("--min-edge-weight", type=float, default=0.01, metavar="FLOAT",
                        help="Minimum edge weight to include in output.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.expression is None:
        _fatal("Provide --expression or use --demo.")
    for p in [args.expression, args.tf_list]:
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
            expression_matrix=Path(args.expression) if args.expression else None,
            tf_list=Path(args.tf_list) if args.tf_list else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            method=args.method,
            min_edge_weight=args.min_edge_weight,
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
    edges = result.get("edges", [])
    print(f"\nGene Regulatory Network Builder v{result['pipeline_version']}")
    print(f"TF regulators: {result.get('n_tf_regulators', 'N/A')}")
    print(f"Regulatory edges: {len(edges)}")
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
