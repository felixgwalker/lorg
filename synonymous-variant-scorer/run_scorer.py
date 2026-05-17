#!/usr/bin/env python3
"""CLI for the Synonymous Variant Scorer.

Usage:
    python run_scorer.py variants.vcf --output-dir results/
    python run_scorer.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser."""
    parser = argparse.ArgumentParser(
        prog="run_scorer.py",
        description=(
            "Synonymous Variant Scorer — scores synonymous variants across four "
            "functional mechanisms: splicing disruption, codon usage bias, "
            "mRNA stability, and cotranslational folding."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "vcf",
        metavar="VCF",
        nargs="?",
        default=None,
        help="VCF file containing synonymous variants. Omit when using --demo.",
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Generate synthetic variants and run the full pipeline without real input.",
    )
    parser.add_argument(
        "-o", "--output-dir",
        metavar="DIR",
        default="./results",
        help="Directory for output files. Created if absent.",
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
    """Validate argument values; exits on any violation."""
    if not args.demo and args.vcf is None:
        _fatal("Provide a VCF file or use --demo.")
    if not args.demo and args.vcf and not Path(args.vcf).exists():
        _fatal(f"VCF file not found: {args.vcf}")


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
            no_plot=args.no_plot,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        logging.exception("Unhandled exception")
        return 1

    _print_summary(result)
    return 0


def _print_summary(result: dict) -> None:
    scored = result["scored"]
    print(f"\nSynonymous Variant Scorer v{result['pipeline_version']}")
    print(f"Variants scored : {len(scored)}")
    if scored:
        tiers = {"HIGH": 0, "MODERATE": 0, "LOW": 0}
        for s in scored:
            tiers[s.impact_tier] = tiers.get(s.impact_tier, 0) + 1
        for tier, count in sorted(tiers.items()):
            print(f"  {tier:10s}: {count}")
        top = max(scored, key=lambda s: s.composite_score)
        print(f"Top variant     : {top.variant_id} ({top.gene}) — composite={top.composite_score:.3f} [{top.impact_tier}]")
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
