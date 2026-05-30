#!/usr/bin/env python3
"""CLI for the Guide RNA Secondary Structure Analyser.

Usage:
    python run_analyser.py guides.fa --output-dir results/
    python run_analyser.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_analyser.py",
        description=(
            "Guide RNA Secondary Structure Analyser — predicts RNA secondary "
            "structure of sgRNA spacer+scaffold, reports MFE, seed-region "
            "accessibility, and guide-scaffold duplex formation flags."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("guides", metavar="GUIDES", nargs="?", default=None,
                        help="FASTA or TSV of 20 nt spacer sequences. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Use synthetic guide sequences.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--scaffold", default="SpCas9",
                        choices=["SpCas9", "SaCas9"],
                        help="sgRNA scaffold sequence identifier.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip dot-bracket structure plot.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.guides is None:
        _fatal("Provide a guides file or use --demo.")
    if not args.demo and args.guides and not Path(args.guides).exists():
        _fatal(f"Guides file not found: {args.guides}")


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
            scaffold=args.scaffold,
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
    structs = result.get("structures", [])
    print(f"\nGuide RNA Secondary Structure Analyser v{result['pipeline_version']}")
    print(f"Guides analysed: {len(structs)}")
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
