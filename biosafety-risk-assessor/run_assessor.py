#!/usr/bin/env python3
"""CLI for the Biosafety Risk Assessor.

Usage:
    python run_assessor.py sequences.fa --output-dir results/
    python run_assessor.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_assessor.py",
        description=(
            "Biosafety Risk Assessor — screens synthetic biology sequences against "
            "select agent, virulence factor, antibiotic resistance, and toxin databases "
            "to recommend biosafety containment levels."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("sequences", metavar="FASTA", nargs="?", default=None,
                        help="DNA or protein FASTA of sequences to assess. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--no-select-agents", action="store_true",
                        help="Skip select agent database screening.")
    parser.add_argument("--no-virulence", action="store_true",
                        help="Skip VFDB virulence factor screening.")
    parser.add_argument("--no-resistance", action="store_true",
                        help="Skip CARD antibiotic resistance screening.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.sequences is None:
        _fatal("Provide a sequence FASTA or use --demo.")
    if args.sequences and not Path(args.sequences).exists():
        _fatal(f"FASTA not found: {args.sequences}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            sequence_fasta=Path(args.sequences) if args.sequences else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            screen_select_agents=not args.no_select_agents,
            screen_virulence=not args.no_virulence,
            screen_resistance=not args.no_resistance,
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
    assessments = result.get("assessments", [])
    failed = sum(1 for a in assessments if isinstance(a, dict) and not a.get("passes_screening"))
    print(f"\nBiosafety Risk Assessor v{result['pipeline_version']}")
    print(f"Sequences assessed: {len(assessments)}")
    print(f"Sequences with flags: {failed}")
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
