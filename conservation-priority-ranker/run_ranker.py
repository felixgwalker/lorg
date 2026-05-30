#!/usr/bin/env python3
"""CLI for the Conservation Priority Ranker.

Usage:
    python run_ranker.py --vcfs pop1.vcf pop2.vcf pop3.vcf --output-dir results/
    python run_ranker.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_ranker.py",
        description=(
            "Conservation Priority Ranker — ranks populations by conservation urgency "
            "using a composite genomic priority score combining inbreeding, Ne, adaptive "
            "diversity, unique alleles, and threat status."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--vcfs", nargs="*", metavar="VCF", default=None,
                        help="Population VCFs to rank.")
    parser.add_argument("--metadata", metavar="TSV", default=None,
                        help="Population metadata TSV (size, threat status, etc.).")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and not args.vcfs:
        _fatal("Provide --vcfs or use --demo.")
    for p in (args.vcfs or []):
        if not Path(p).exists():
            _fatal(f"VCF not found: {p}")
    if args.metadata and not Path(args.metadata).exists():
        _fatal(f"Metadata file not found: {args.metadata}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            vcfs=[Path(v) for v in args.vcfs] if args.vcfs else None,
            population_metadata=Path(args.metadata) if args.metadata else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
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
    ranked = result.get("ranked_populations", [])
    print(f"\nConservation Priority Ranker v{result['pipeline_version']}")
    print(f"Populations ranked: {len(ranked)}")
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
