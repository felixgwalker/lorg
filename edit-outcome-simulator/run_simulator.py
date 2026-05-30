#!/usr/bin/env python3
"""CLI for the Edit Outcome Simulator.

Usage:
    python run_simulator.py target.fa --guide SPACER --output-dir results/
    python run_simulator.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_simulator.py",
        description=(
            "Edit Outcome Simulator — simulates CRISPR indel distributions at a "
            "target site using an inDelphi-style approach: 1-bp templated insertions "
            "and microhomology-weighted deletions up to 30 bp."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("fasta", metavar="FASTA", nargs="?", default=None,
                        help="FASTA of ≥60 bp flanking the cut site. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Use a synthetic target sequence.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--guide", default=None, metavar="SEQ",
                        help="20 nt guide spacer to locate cut site in FASTA.")
    parser.add_argument("--cas", default="SpCas9", metavar="VARIANT",
                        choices=["SpCas9", "SaCas9", "Cas12a", "Cas12b"],
                        help="Cas variant for cut position offset.")
    parser.add_argument("--n-simulations", type=int, default=10000, metavar="INT",
                        help="Number of Monte Carlo draws.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip stacked bar chart generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.fasta is None:
        _fatal("Provide a FASTA file or use --demo.")
    if not args.demo and args.fasta and not Path(args.fasta).exists():
        _fatal(f"FASTA not found: {args.fasta}")
    if args.n_simulations < 100:
        _fatal("--n-simulations must be >= 100.")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            target_fasta=Path(args.fasta) if args.fasta else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            guide=args.guide,
            cas_variant=args.cas,
            n_simulations=args.n_simulations,
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
    print(f"\nEdit Outcome Simulator v{result['pipeline_version']}")
    print(f"Frameshift rate: {result.get('frameshift_rate', 'N/A')}")
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
