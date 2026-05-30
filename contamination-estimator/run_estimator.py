#!/usr/bin/env python3
"""CLI for the Contamination Estimator.

Usage:
    python run_estimator.py sample.bam --reference genome.fa --output-dir results/
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
            "Contamination Estimator — estimates the fraction of modern human "
            "contamination in an ancient DNA BAM using mitochondrial consensus deviation, "
            "X-chromosome heterozygosity, and ANGSD-based approaches."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("bam", metavar="BAM", nargs="?", default=None,
                        help="BAM file of ancient DNA reads. Omit with --demo.")
    parser.add_argument("--reference", metavar="FASTA", default=None,
                        help="Reference genome FASTA.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--methods", nargs="*",
                        choices=["mt_consensus", "nuclear_X", "ANGSD", "schmutzi"],
                        default=None, help="Contamination methods to use (default: all).")
    parser.add_argument("--contamination-threshold", type=float, default=0.03, metavar="FLOAT",
                        help="Contamination rate above which sample fails QC.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.bam is None:
        _fatal("Provide a BAM file or use --demo.")
    if args.bam and not Path(args.bam).exists():
        _fatal(f"BAM not found: {args.bam}")
    if not args.demo and args.reference and not Path(args.reference).exists():
        _fatal(f"Reference FASTA not found: {args.reference}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            bam_path=Path(args.bam) if args.bam else None,
            reference_fasta=Path(args.reference) if args.reference else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            methods=args.methods,
            contamination_threshold=args.contamination_threshold,
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
    print(f"\nContamination Estimator v{result['pipeline_version']}")
    print(f"Combined estimate: {result.get('combined_estimate', 'N/A')}")
    print(f"Passes QC: {result.get('passes_threshold', 'N/A')}")
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
