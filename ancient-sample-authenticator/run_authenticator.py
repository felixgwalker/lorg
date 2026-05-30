#!/usr/bin/env python3
"""CLI for the Ancient Sample Authenticator.

Usage:
    python run_authenticator.py sample.bam --reference genome.fa --output-dir results/
    python run_authenticator.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_authenticator.py",
        description=(
            "Ancient Sample Authenticator — evaluates fragment length, deamination "
            "damage, contamination, endogenous DNA fraction, and coverage to produce "
            "a composite authentication verdict for an ancient DNA BAM."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("bam", metavar="BAM", nargs="?", default=None,
                        help="BAM file of aligned ancient DNA reads. Omit with --demo.")
    parser.add_argument("--reference", metavar="FASTA", default=None,
                        help="Reference genome FASTA.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--sample-id", default="sample", metavar="ID",
                        help="Sample identifier for the report.")
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
            sample_id=args.sample_id,
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
    res = result.get("result")
    print(f"\nAncient Sample Authenticator v{result['pipeline_version']}")
    if res:
        print(f"Verdict: {res.get('verdict', 'N/A')}")
        print(f"Composite score: {res.get('composite_score', 'N/A'):.3f}")
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
