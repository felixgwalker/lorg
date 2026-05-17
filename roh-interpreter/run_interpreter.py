#!/usr/bin/env python3
"""CLI for the ROH Interpreter.

Usage:
    python run_interpreter.py genotypes.vcf --output-dir results/
    python run_interpreter.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser."""
    parser = argparse.ArgumentParser(
        prog="run_interpreter.py",
        description=(
            "ROH Interpreter — detects Runs of Homozygosity from VCF genotype data, "
            "computes FROH per individual, estimates effective population size (Ne), "
            "and outputs a ROH BED catalog, FROH CSV, and analysis plots."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "vcf",
        metavar="VCF",
        nargs="?",
        default=None,
        help="VCF file with GT genotype fields. Omit when using --demo.",
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Generate synthetic data and run full pipeline without real input.",
    )
    parser.add_argument(
        "-o", "--output-dir",
        metavar="DIR",
        default="./results",
        help="Directory for output files. Created if absent.",
    )
    parser.add_argument(
        "--window-size",
        metavar="INT",
        type=int,
        default=50,
        help="Sliding window size in SNPs for ROH detection.",
    )
    parser.add_argument(
        "--homo-threshold",
        metavar="FLOAT",
        type=float,
        default=0.95,
        help="Minimum homozygosity fraction in window to call ROH.",
    )
    parser.add_argument(
        "--genome-length",
        metavar="INT",
        type=int,
        default=2_700_000_000,
        help="Autosomal genome length in bp for FROH denominator.",
    )
    parser.add_argument(
        "--generation-time",
        metavar="FLOAT",
        type=float,
        default=6.0,
        help="Generation time in years for Ne estimation.",
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
    if not args.demo and args.vcf is None:
        _fatal("Provide a VCF file or use --demo.")
    if not args.demo and args.vcf and not Path(args.vcf).exists():
        _fatal(f"VCF file not found: {args.vcf}")
    if not (0.5 <= args.homo_threshold <= 1.0):
        _fatal(f"--homo-threshold must be in [0.5, 1.0], got {args.homo_threshold}.")
    if args.window_size < 5:
        _fatal(f"--window-size must be >= 5, got {args.window_size}.")


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
            vcf_path=Path(args.vcf) if args.vcf else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            window_size=args.window_size,
            homo_threshold=args.homo_threshold,
            genome_length=args.genome_length,
            generation_time=args.generation_time,
            no_plot=args.no_plot,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        logging.exception("Unhandled exception")
        return 1

    _print_summary(result)
    return 0


def _print_summary(result: dict) -> None:
    segs = result["segments"]
    frohs = result["froh_results"]
    print(f"\nROH Interpreter v{result['pipeline_version']}")
    print(f"ROH segments detected : {len(segs)}")
    if segs:
        classes = {"SHORT": 0, "MEDIUM": 0, "LONG": 0, "VERY_SHORT": 0}
        for s in segs:
            classes[s.length_class] = classes.get(s.length_class, 0) + 1
        for cls, cnt in sorted(classes.items()):
            print(f"  {cls:12s}: {cnt}")
    print(f"Individuals analyzed  : {len(frohs)}")
    for fr in frohs:
        print(f"  {fr.individual_id}: FROH={fr.froh:.4f}, n_ROH={fr.n_roh}")
    out = result.get("output_files", {})
    print("Output files:")
    for k, v in out.items():
        print(f"  {k:15s}: {v}")
    print()


def _fatal(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
