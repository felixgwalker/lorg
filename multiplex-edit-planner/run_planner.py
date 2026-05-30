#!/usr/bin/env python3
"""CLI for the Multiplex Edit Planner.

Usage:
    python run_planner.py manifest.json --output-dir results/
    python run_planner.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_planner.py",
        description=(
            "Multiplex Edit Planner — checks guide cross-reactivity, cut-window "
            "overlap, and translocation risk for a set of simultaneous edits, "
            "then batches guides into compatible delivery groups."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("manifest", metavar="MANIFEST", nargs="?", default=None,
                        help="JSON or TSV edit manifest. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on a synthetic 6-target manifest.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--max-guides-per-batch", type=int, default=4, metavar="INT",
                        help="Maximum guides per delivery batch.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.manifest is None:
        _fatal("Provide a manifest file or use --demo.")
    if not args.demo and args.manifest and not Path(args.manifest).exists():
        _fatal(f"Manifest not found: {args.manifest}")
    if args.max_guides_per_batch < 1:
        _fatal("--max-guides-per-batch must be >= 1.")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            manifest_path=Path(args.manifest) if args.manifest else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            max_guides_per_batch=args.max_guides_per_batch,
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
    batches = result.get("batches", [])
    print(f"\nMultiplex Edit Planner v{result['pipeline_version']}")
    print(f"Delivery batches: {len(batches)}")
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
