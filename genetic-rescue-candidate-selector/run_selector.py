#!/usr/bin/env python3
"""CLI for the Genetic Rescue Candidate Selector.

Usage:
    python run_selector.py --recipient recipient.vcf --donors donor1.vcf donor2.vcf --output-dir results/
    python run_selector.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_selector.py",
        description=(
            "Genetic Rescue Candidate Selector — ranks donor populations for genetic "
            "rescue of an inbred recipient by balancing heterozygosity gain, kinship "
            "distance, ecotype compatibility, and outbreeding depression risk."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--recipient", metavar="VCF", default=None,
                        help="VCF of the recipient population.")
    parser.add_argument("--donors", nargs="*", metavar="VCF", default=None,
                        help="VCFs of candidate donor populations.")
    parser.add_argument("--pop-map", metavar="TSV", default=None,
                        help="TSV mapping samples to populations.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--recipient-pop", default="recipient", metavar="LABEL",
                        help="Population label for the recipient.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.recipient is None:
        _fatal("Provide --recipient or use --demo.")
    for p in ([args.recipient] if args.recipient else []) + (args.donors or []):
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
            recipient_vcf=Path(args.recipient) if args.recipient else None,
            donor_vcfs=[Path(d) for d in args.donors] if args.donors else None,
            population_map=Path(args.pop_map) if args.pop_map else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            recipient_pop=args.recipient_pop,
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
    candidates = result.get("candidates", [])
    print(f"\nGenetic Rescue Candidate Selector v{result['pipeline_version']}")
    print(f"Donor candidates ranked: {len(candidates)}")
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
