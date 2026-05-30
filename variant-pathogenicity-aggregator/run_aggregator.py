#!/usr/bin/env python3
"""CLI for the Variant Pathogenicity Aggregator.

Usage:
    python run_aggregator.py variants.vcf --output-dir results/
    python run_aggregator.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_aggregator.py",
        description=(
            "Variant Pathogenicity Aggregator — collects evidence from ClinVar and "
            "in silico tools, maps it to ACMG/AMP criteria, and produces a five-tier "
            "composite pathogenicity classification."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="Annotated VCF with ClinVar and tool score INFO fields. Omit with --demo.")
    parser.add_argument("--tool-scores", metavar="TSV", default=None,
                        help="Optional TSV of additional in silico tool scores per variant.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--clinvar-stars", type=int, default=1, metavar="INT",
                        help="Minimum ClinVar review star count to trust ClinVar evidence.")
    parser.add_argument("--classification-threshold", type=int, default=4, metavar="INT",
                        help="Minimum ACMG point total to call likely-pathogenic.")
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
    if args.tool_scores and not Path(args.tool_scores).exists():
        _fatal(f"Tool scores TSV not found: {args.tool_scores}")


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
            tool_scores=Path(args.tool_scores) if args.tool_scores else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            clinvar_stars=args.clinvar_stars,
            classification_threshold=args.classification_threshold,
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
    results = result.get("results", [])
    print(f"\nVariant Pathogenicity Aggregator v{result['pipeline_version']}")
    print(f"Variants classified: {len(results)}")
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
