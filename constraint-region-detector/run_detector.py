#!/usr/bin/env python3
"""CLI for the Constraint Region Detector.

Usage:
    python run_detector.py variants.vcf --constraint-file gnomad_constraint.tsv --output-dir results/
    python run_detector.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_detector.py",
        description=(
            "Constraint Region Detector — intersects a VCF with gnomAD constraint "
            "metrics (LOEUF, pLI, missense Z-score) to flag variants in "
            "genomically constrained genes and regions."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="VCF of variants to evaluate. Omit with --demo.")
    parser.add_argument("--constraint-file", metavar="TSV", default=None,
                        help="gnomAD-format constraint TSV.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--loeuf-threshold", type=float, default=0.35, metavar="FLOAT",
                        help="LOEUF upper bound; genes below are flagged as constrained.")
    parser.add_argument("--z-score-threshold", type=float, default=3.09, metavar="FLOAT",
                        help="Missense Z-score above which a gene is constrained.")
    parser.add_argument("--primary-metric", choices=["LOEUF", "pLI", "z_score"],
                        default="LOEUF", help="Primary constraint metric.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.vcf is None:
        _fatal("Provide a VCF file or use --demo.")
    if not args.demo and args.vcf and not Path(args.vcf).exists():
        _fatal(f"VCF not found: {args.vcf}")
    if not args.demo and args.constraint_file and not Path(args.constraint_file).exists():
        _fatal(f"Constraint file not found: {args.constraint_file}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            vcf_path=Path(args.vcf) if args.vcf else None,
            constraint_file=Path(args.constraint_file) if args.constraint_file else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            loeuf_threshold=args.loeuf_threshold,
            z_score_threshold=args.z_score_threshold,
            primary_metric=args.primary_metric,
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
    overlaps = result.get("overlaps", [])
    constrained = sum(1 for o in overlaps if isinstance(o, dict) and o.get("is_constrained"))
    print(f"\nConstraint Region Detector v{result['pipeline_version']}")
    print(f"Variants evaluated: {len(overlaps)}")
    print(f"In constrained regions: {constrained}")
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
