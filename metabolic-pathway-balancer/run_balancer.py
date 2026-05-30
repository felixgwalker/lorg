#!/usr/bin/env python3
"""CLI for the Metabolic Pathway Balancer.

Usage:
    python run_balancer.py --pathway pathway.json --output-dir results/
    python run_balancer.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_balancer.py",
        description=(
            "Metabolic Pathway Balancer — runs FBA on a heterologous pathway, "
            "identifies flux-limiting bottlenecks and cofactor imbalances, and "
            "recommends enzyme expression adjustments."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--pathway", metavar="JSON", default=None,
                        help="JSON of pathway reactions with stoichiometry and kinetics.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--objective",
                        choices=["maximise_yield", "maximise_productivity",
                                 "minimise_byproducts"],
                        default="maximise_yield", help="Optimisation objective.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.pathway is None:
        _fatal("Provide --pathway or use --demo.")
    if args.pathway and not Path(args.pathway).exists():
        _fatal(f"Pathway JSON not found: {args.pathway}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            pathway_json=Path(args.pathway) if args.pathway else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            objective=args.objective,
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
    bp = result.get("balanced_pathway")
    print(f"\nMetabolic Pathway Balancer v{result['pipeline_version']}")
    if bp:
        print(f"Bottlenecks identified: {len(bp.get('bottlenecks', []))}")
        print(f"Predicted yield: {bp.get('predicted_yield', 'N/A')}")
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
