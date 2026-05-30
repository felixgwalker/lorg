#!/usr/bin/env python3
"""CLI for the Guide RNA GC Optimiser.

Usage:
    python run_optimiser.py guides.fa --output-dir results/
    python run_optimiser.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_optimiser.py",
        description=(
            "Guide RNA GC Optimiser — scores and ranks sgRNAs by GC content "
            "features: total GC (30–70 %), seed-region GC (40–60 %), "
            "homopolymer runs, and poly-T stretches."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("guides", metavar="GUIDES", nargs="?", default=None,
                        help="FASTA or TSV of guide sequences. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic guides.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--gc-min", type=float, default=0.30, metavar="FLOAT",
                        help="Minimum acceptable GC fraction.")
    parser.add_argument("--gc-max", type=float, default=0.70, metavar="FLOAT",
                        help="Maximum acceptable GC fraction.")
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
    if not (0.0 <= args.gc_min < args.gc_max <= 1.0):
        _fatal("--gc-min must be less than --gc-max and both in [0, 1].")


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
            output_dir=Path(args.output_dir),
            demo=args.demo,
            gc_min=args.gc_min,
            gc_max=args.gc_max,
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
    guides = result.get("scored_guides", [])
    print(f"\nGuide RNA GC Optimiser v{result['pipeline_version']}")
    print(f"Guides scored: {len(guides)}")
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
