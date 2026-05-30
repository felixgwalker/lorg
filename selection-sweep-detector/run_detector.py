#!/usr/bin/env python3
"""CLI for the Selection Sweep Detector.

Usage:
    python run_detector.py variants.vcf --pop-map populations.tsv --output-dir results/
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
            "Selection Sweep Detector — identifies genomic regions under recent positive "
            "selection using iHS, XP-EHH, the composite likelihood ratio (CLR), and "
            "Tajima's D in sliding windows."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="Phased multi-population VCF. Omit with --demo.")
    parser.add_argument("--pop-map", metavar="TSV", default=None,
                        help="TSV mapping sample IDs to population labels.")
    parser.add_argument("--gene-annotation", metavar="BED", default=None,
                        help="Gene annotation BED for candidate gene reporting.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--tests", nargs="*",
                        choices=["iHS", "XP-EHH", "CLR", "tajimas_d"],
                        default=None, help="Tests to run (default: all).")
    parser.add_argument("--outlier-percentile", type=float, default=99.0, metavar="FLOAT",
                        help="Percentile threshold for outlier sweep windows.")
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
    if not args.demo and args.pop_map and not Path(args.pop_map).exists():
        _fatal(f"Population map not found: {args.pop_map}")
    if args.gene_annotation and not Path(args.gene_annotation).exists():
        _fatal(f"Gene annotation not found: {args.gene_annotation}")


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
            population_map=Path(args.pop_map) if args.pop_map else None,
            gene_annotation=Path(args.gene_annotation) if args.gene_annotation else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            tests=args.tests,
            outlier_percentile=args.outlier_percentile,
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
    sweeps = result.get("sweep_regions", [])
    print(f"\nSelection Sweep Detector v{result['pipeline_version']}")
    print(f"Sweep regions detected: {len(sweeps)}")
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
