#!/usr/bin/env python3
"""CLI for the Guide RNA Specificity Ranker.

Usage:
    python run_ranker.py guides.fa --genome genome.fa --output-dir results/
    python run_ranker.py guides.fa --offtargets sites.bed --output-dir results/
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
            "Guide RNA Specificity Ranker — ranks sgRNAs by predicted specificity "
            "using CFD (Cutting Frequency Determination) scoring against enumerated "
            "off-target sites with up to 3 mismatches."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("guides", metavar="GUIDES", nargs="?", default=None,
                        help="FASTA or TSV of 20 nt guides. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic guides.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--genome", default=None, metavar="FASTA",
                        help="Reference genome FASTA for off-target search.")
    parser.add_argument("--offtargets", default=None, metavar="BED",
                        help="Pre-computed off-target BED (alternative to --genome).")
    parser.add_argument("--max-mismatches", type=int, default=3, metavar="INT",
                        help="Maximum mismatches to enumerate.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.guides is None:
        _fatal("Provide a guides file or use --demo.")
    if not args.demo and args.guides and not Path(args.guides).exists():
        _fatal(f"Guides file not found: {args.guides}")
    if args.genome and not Path(args.genome).exists():
        _fatal(f"Genome FASTA not found: {args.genome}")
    if args.offtargets and not Path(args.offtargets).exists():
        _fatal(f"Off-targets BED not found: {args.offtargets}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            guides_path=Path(args.guides) if args.guides else None,
            genome_fasta=Path(args.genome) if args.genome else None,
            offtargets_bed=Path(args.offtargets) if args.offtargets else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            max_mismatches=args.max_mismatches,
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
    guides = result.get("ranked_guides", [])
    print(f"\nGuide RNA Specificity Ranker v{result['pipeline_version']}")
    print(f"Guides ranked: {len(guides)}")
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
