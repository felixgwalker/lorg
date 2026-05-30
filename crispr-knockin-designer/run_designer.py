#!/usr/bin/env python3
"""CLI for the CRISPR Knock-in Designer.

Usage:
    python run_designer.py locus.fa --insert insert.fa --output-dir results/
    python run_designer.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_designer.py",
        description=(
            "CRISPR Knock-in Designer — designs sgRNA and HDR donor template "
            "for precision knock-in, with configurable homology arm length and "
            "silent PAM mutation to prevent re-cutting."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("locus", metavar="LOCUS", nargs="?", default=None,
                        help="Target locus FASTA (≥500 bp). Omit with --demo.")
    parser.add_argument("--insert", metavar="FASTA", default=None,
                        help="Insert sequence FASTA.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic sequences.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--pam", default="NGG", metavar="PAM",
                        help="PAM motif.")
    parser.add_argument("--arm-length", type=int, default=100, metavar="INT",
                        help="Homology arm length in bp.")
    parser.add_argument("--max-cut-distance", type=int, default=30, metavar="INT",
                        help="Max bp from insertion point to cut site.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.locus is None:
        _fatal("Provide a locus FASTA or use --demo.")
    if not args.demo and args.locus and not Path(args.locus).exists():
        _fatal(f"Locus FASTA not found: {args.locus}")
    if args.insert and not Path(args.insert).exists():
        _fatal(f"Insert FASTA not found: {args.insert}")
    if args.arm_length < 50:
        _fatal("--arm-length must be >= 50 bp.")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            locus_fasta=Path(args.locus) if args.locus else None,
            insert_fasta=Path(args.insert) if args.insert else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            pam=args.pam,
            arm_length=args.arm_length,
            max_cut_distance=args.max_cut_distance,
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
    print(f"\nCRISPR Knock-in Designer v{result['pipeline_version']}")
    print(f"Guide: {result.get('guide', {})}")
    print(f"Donor length: {len(result.get('donor_sequence', ''))}")
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
