#!/usr/bin/env python3
"""CLI for the Missense Impact Scorer.

Usage:
    python run_scorer.py variants.vcf --output-dir results/
    python run_scorer.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_scorer.py",
        description=(
            "Missense Impact Scorer — scores the functional impact of missense variants "
            "using conservation, BLOSUM62 substitution cost, and physicochemical property "
            "changes, producing a composite five-tier pathogenicity estimate."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="VCF of missense variants with HGVS_P annotations. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--conservation-tool",
                        choices=["phyloP", "GERP", "phastCons"], default="phyloP",
                        help="Conservation metric to use.")
    parser.add_argument("--benign-threshold", type=float, default=0.3, metavar="FLOAT",
                        help="Composite score below which variants are called benign.")
    parser.add_argument("--pathogenic-threshold", type=float, default=0.7, metavar="FLOAT",
                        help="Composite score above which variants are called pathogenic.")
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
    if args.benign_threshold >= args.pathogenic_threshold:
        _fatal("--benign-threshold must be less than --pathogenic-threshold.")


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
            output_dir=Path(args.output_dir),
            demo=args.demo,
            conservation_tool=args.conservation_tool,
            benign_threshold=args.benign_threshold,
            pathogenic_threshold=args.pathogenic_threshold,
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
    variants = result.get("scored_variants", [])
    print(f"\nMissense Impact Scorer v{result['pipeline_version']}")
    print(f"Variants scored: {len(variants)}")
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
