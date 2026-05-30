#!/usr/bin/env python3
"""CLI for the Founder Effect Estimator.

Usage:
    python run_estimator.py --study-vcf study.vcf --reference-vcf ref.vcf --output-dir results/
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
            "Founder Effect Estimator — compares nucleotide diversity, private variant "
            "fraction, and haplotype block length between a study and reference population "
            "to infer whether and how strongly a founder event occurred."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--study-vcf", metavar="VCF", default=None,
                        help="VCF of the study (potentially founded) population.")
    parser.add_argument("--reference-vcf", metavar="VCF", default=None,
                        help="VCF of the ancestral/reference population.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--study-pop", default="study", metavar="LABEL",
                        help="Label for the study population.")
    parser.add_argument("--reference-pop", default="reference", metavar="LABEL",
                        help="Label for the reference population.")
    parser.add_argument("--generation-time", type=float, default=30.0, metavar="YEARS",
                        help="Assumed generation time in years.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and (args.study_vcf is None or args.reference_vcf is None):
        _fatal("Provide --study-vcf and --reference-vcf, or use --demo.")
    if args.study_vcf and not Path(args.study_vcf).exists():
        _fatal(f"Study VCF not found: {args.study_vcf}")
    if args.reference_vcf and not Path(args.reference_vcf).exists():
        _fatal(f"Reference VCF not found: {args.reference_vcf}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            study_vcf=Path(args.study_vcf) if args.study_vcf else None,
            reference_vcf=Path(args.reference_vcf) if args.reference_vcf else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            study_pop=args.study_pop,
            reference_pop=args.reference_pop,
            generation_time=args.generation_time,
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
    print(f"\nFounder Effect Estimator v{result['pipeline_version']}")
    est = result.get("estimate")
    if est:
        print(f"Founder signature: {est.get('signature', 'N/A')}")
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
