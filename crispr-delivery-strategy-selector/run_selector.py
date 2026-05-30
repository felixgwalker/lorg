#!/usr/bin/env python3
"""CLI for the CRISPR Delivery Strategy Selector.

Usage:
    python run_selector.py delivery_spec.json --output-dir results/
    python run_selector.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_selector.py",
        description=(
            "CRISPR Delivery Strategy Selector — ranks delivery modalities "
            "(plasmid, RNP, LNP, AAV, lentivirus) for a given cell type and "
            "payload using a rule-based compatibility scoring table."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("spec", metavar="SPEC_JSON", nargs="?", default=None,
                        help="Delivery specification JSON. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Use a synthetic HEK293 editing scenario.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip bar chart generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.spec is None:
        _fatal("Provide a spec JSON or use --demo.")
    if not args.demo and args.spec and not Path(args.spec).exists():
        _fatal(f"Spec JSON not found: {args.spec}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            spec_json=Path(args.spec) if args.spec else None,
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
    strategies = result.get("strategies", [])
    print(f"\nCRISPR Delivery Strategy Selector v{result['pipeline_version']}")
    print(f"Strategies ranked: {len(strategies)}")
    print(f"Recommended: {result.get('recommended', 'N/A')}")
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
