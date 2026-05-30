#!/usr/bin/env python3
"""CLI for the Admixture Signal Scanner.

Usage:
    python run_scanner.py variants.vcf --output-dir results/
    python run_scanner.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_scanner.py",
        description=(
            "Admixture Signal Scanner — fits parametric admixture models for K = k_min…k_max, "
            "selects best K by cross-validation, and reports per-sample ancestry proportions "
            "as a structure-style bar chart."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="VCF of LD-pruned SNPs. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--k-min", type=int, default=2, metavar="INT",
                        help="Minimum number of ancestry components to test.")
    parser.add_argument("--k-max", type=int, default=10, metavar="INT",
                        help="Maximum number of ancestry components to test.")
    parser.add_argument("--model", choices=["unsupervised", "supervised"],
                        default="unsupervised", help="Admixture model type.")
    parser.add_argument("--reference-panel", metavar="VCF", default=None,
                        help="Reference population panel VCF (supervised mode only).")
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
    if args.k_min >= args.k_max:
        _fatal("--k-min must be less than --k-max.")
    if args.model == "supervised" and args.reference_panel is None:
        _fatal("--reference-panel is required for supervised mode.")
    if args.reference_panel and not Path(args.reference_panel).exists():
        _fatal(f"Reference panel not found: {args.reference_panel}")


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
            k_min=args.k_min,
            k_max=args.k_max,
            model=args.model,
            reference_panel=Path(args.reference_panel) if args.reference_panel else None,
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
    print(f"\nAdmixture Signal Scanner v{result['pipeline_version']}")
    print(f"Best K: {result.get('best_k', 'N/A')}")
    print(f"Samples analysed: {len(result.get('sample_ancestries', []))}")
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
