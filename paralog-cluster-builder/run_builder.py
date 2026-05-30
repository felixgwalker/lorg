#!/usr/bin/env python3
"""CLI for the Paralog Cluster Builder.

Usage:
    python run_builder.py proteome.fa --output-dir results/
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
            "Paralog Cluster Builder — performs all-vs-self BLAST on a proteome and "
            "clusters paralog pairs with MCL to produce gene family clusters, annotated "
            "by duplication type (tandem, segmental, dispersed, retroposed)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("proteome", metavar="FASTA", nargs="?", default=None,
                        help="Single-species protein FASTA. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--min-identity", type=float, default=30.0, metavar="FLOAT",
                        help="Minimum sequence identity (%%) for paralog pairs.")
    parser.add_argument("--inflation", type=float, default=2.0, metavar="FLOAT",
                        help="MCL inflation parameter (higher = more clusters).")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.proteome is None:
        _fatal("Provide a protein FASTA or use --demo.")
    if args.proteome and not Path(args.proteome).exists():
        _fatal(f"Proteome FASTA not found: {args.proteome}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            proteome=Path(args.proteome) if args.proteome else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            min_identity=args.min_identity,
            inflation=args.inflation,
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
    clusters = result.get("clusters", [])
    print(f"\nParalog Cluster Builder v{result['pipeline_version']}")
    print(f"Paralog clusters built: {len(clusters)}")
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
