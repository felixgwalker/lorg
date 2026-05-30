#!/usr/bin/env python3
"""CLI for the Repeat Element Classifier.

Usage:
    python run_classifier.py assembly.fa --output-dir results/
    python run_classifier.py --repeat-masker-out assembly.fa.out --output-dir results/
    python run_classifier.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_classifier.py",
        description=(
            "Repeat Element Classifier — classifies repeat elements in a genome assembly "
            "by parsing RepeatMasker output or running de novo repeat identification, "
            "producing a BED annotation and repeat landscape plot."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("assembly", metavar="FASTA", nargs="?", default=None,
                        help="Assembly FASTA to annotate. Omit with --demo.")
    parser.add_argument("--repeat-masker-out", metavar="OUT", default=None,
                        help="Pre-existing RepeatMasker .out file.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--min-length", type=int, default=50, metavar="INT",
                        help="Minimum element length to include.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.assembly is None and args.repeat_masker_out is None:
        _fatal("Provide an assembly FASTA or --repeat-masker-out, or use --demo.")
    if args.assembly and not Path(args.assembly).exists():
        _fatal(f"Assembly FASTA not found: {args.assembly}")
    if args.repeat_masker_out and not Path(args.repeat_masker_out).exists():
        _fatal(f"RepeatMasker output not found: {args.repeat_masker_out}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            assembly_fasta=Path(args.assembly) if args.assembly else None,
            repeat_masker_out=Path(args.repeat_masker_out) if args.repeat_masker_out else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            min_length=args.min_length,
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
    elements = result.get("elements", [])
    print(f"\nRepeat Element Classifier v{result['pipeline_version']}")
    print(f"Repeat elements classified: {len(elements)}")
    print(f"Total repeat fraction: {result.get('total_repeat_fraction', 'N/A'):.1%}")
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
