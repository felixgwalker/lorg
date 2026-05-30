#!/usr/bin/env python3
"""CLI for the CRISPR Array Designer.

Usage:
    python run_designer.py targets.fa --cas Cas12a --output-dir results/
    python run_designer.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_designer.py",
        description=(
            "CRISPR Array Designer — assembles a CRISPR array from multiple "
            "target sequences by identifying PAM sites, extracting spacers, "
            "checking uniqueness, and interspersing system-specific direct repeats."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("targets", metavar="TARGETS", nargs="?", default=None,
                        help="FASTA or TSV of target sequences. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic targets.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--cas", default="Cas12a", metavar="SYSTEM",
                        choices=["Cas12a", "Cas9", "Cas12b"],
                        help="Cas nuclease system.")
    parser.add_argument("--direct-repeat", default=None, metavar="SEQ",
                        help="Override direct repeat sequence.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.targets is None:
        _fatal("Provide a targets file or use --demo.")
    if not args.demo and args.targets and not Path(args.targets).exists():
        _fatal(f"Targets file not found: {args.targets}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            targets_path=Path(args.targets) if args.targets else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            cas_system=args.cas,
            direct_repeat=args.direct_repeat,
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
    spacers = result.get("spacers", [])
    print(f"\nCRISPR Array Designer v{result['pipeline_version']}")
    print(f"Spacers designed: {len(spacers)}")
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
