#!/usr/bin/env python3
"""CLI for the Gene Circuit Stability Estimator.

Usage:
    python run_estimator.py --circuit circuit.json --output-dir results/
    python run_estimator.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_estimator.py",
        description=(
            "Gene Circuit Stability Estimator — simulates ODE dynamics for a synthetic "
            "gene circuit, classifies behaviour (stable/oscillating/bistable/unstable), "
            "and quantifies robustness under parameter perturbations."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--circuit", metavar="JSON", default=None,
                        help="JSON of circuit nodes and edges.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--simulation-time", type=float, default=100.0, metavar="HOURS",
                        help="Simulation duration in hours.")
    parser.add_argument("--n-perturbations", type=int, default=100, metavar="INT",
                        help="Number of parameter perturbations for robustness analysis.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.circuit is None:
        _fatal("Provide --circuit or use --demo.")
    if args.circuit and not Path(args.circuit).exists():
        _fatal(f"Circuit JSON not found: {args.circuit}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            circuit_json=Path(args.circuit) if args.circuit else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            simulation_time=args.simulation_time,
            n_perturbations=args.n_perturbations,
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
    a = result.get("analysis")
    print(f"\nGene Circuit Stability Estimator v{result['pipeline_version']}")
    if a:
        print(f"Behaviour: {a.get('behaviour', 'N/A')}")
        print(f"Robustness score: {a.get('robustness_score', 'N/A'):.3f}")
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
