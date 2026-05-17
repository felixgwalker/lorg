#!/usr/bin/env python3
"""CLI for the Phylogenetic Distance Estimator.

Usage:
    python run_estimator.py alignment.fasta --output-dir results/
    python run_estimator.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser."""
    parser = argparse.ArgumentParser(
        prog="run_estimator.py",
        description=(
            "Phylogenetic Distance Estimator — computes pairwise phylogenetic distances "
            "from a sequence alignment, builds a Neighbor-Joining tree, and renders a "
            "heatmap + dendrogram. Outputs: distance matrix CSV, Newick tree, ranked "
            "comparisons CSV, and PNG/SVG plots."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "fasta",
        metavar="FASTA",
        nargs="?",
        default=None,
        help="FASTA or CLUSTAL alignment file. Omit when using --demo.",
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Generate synthetic alignment and run full pipeline without real input.",
    )
    parser.add_argument(
        "-o", "--output-dir",
        metavar="DIR",
        default="./results",
        help="Directory for output files. Created if absent.",
    )
    parser.add_argument(
        "--model",
        choices=["JC69", "K2P"],
        default="K2P",
        help="Substitution model for distance calculation.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip matplotlib figure generation.",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging verbosity.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 1.0.0",
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    """Validate argument values; exits on violation."""
    if not args.demo and args.fasta is None:
        _fatal("Provide a FASTA alignment file or use --demo.")
    if not args.demo and args.fasta and not Path(args.fasta).exists():
        _fatal(f"FASTA file not found: {args.fasta}")


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s: %(message)s",
        stream=sys.stderr,
    )

    validate_args(args)

    from src.pipeline import run_pipeline

    try:
        result = run_pipeline(
            fasta_path=Path(args.fasta) if args.fasta else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            distance_model=args.model,
            no_plot=args.no_plot,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        logging.exception("Unhandled exception")
        return 1

    _print_summary(result)
    return 0


def _print_summary(result: dict) -> None:
    names = result["names"]
    matrix = result["distance_matrix"]
    print(f"\nPhylogenetic Distance Estimator v{result['pipeline_version']}")
    print(f"Sequences analyzed : {len(names)}")
    print(f"Species            : {', '.join(names)}")
    if len(names) >= 2:
        import numpy as np
        mask = ~np.eye(len(names), dtype=bool)
        print(f"Distance range     : {matrix[mask].min():.4f} – {matrix[mask].max():.4f}")
        closest_idx = np.unravel_index(np.argmin(matrix + np.eye(len(names)) * 999), matrix.shape)
        print(f"Closest pair       : {names[closest_idx[0]]} — {names[closest_idx[1]]} ({matrix[closest_idx]:.4f})")
    print(f"Newick tree        : {result['newick'][:80]}{'...' if len(result['newick']) > 80 else ''}")
    out = result.get("output_files", {})
    print("Output files:")
    for k, v in out.items():
        print(f"  {k:25s}: {v}")
    print()


def _fatal(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
