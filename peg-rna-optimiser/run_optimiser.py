#!/usr/bin/env python3
"""CLI for the pegRNA Optimiser.

Usage:
    python run_optimiser.py target.fa --edit edit.json --output-dir results/
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
            "pegRNA Optimiser — grid-searches PBS and RT template length "
            "combinations to find Pareto-optimal pegRNA designs trading "
            "predicted efficiency against synthesis complexity."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("fasta", metavar="FASTA", nargs="?", default=None,
                        help="Target locus FASTA. Omit with --demo.")
    parser.add_argument("--edit", metavar="JSON", default=None,
                        help="Edit specification JSON.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--pam", default="NGG", metavar="PAM",
                        help="PAM motif.")
    parser.add_argument("--pbs-range", default="8-15", metavar="MIN-MAX",
                        help="PBS length range as MIN-MAX.")
    parser.add_argument("--rt-range", default="10-16", metavar="MIN-MAX",
                        help="RT template length range as MIN-MAX.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip heatmap figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def _parse_range(s: str, label: str) -> tuple[int, int]:
    try:
        lo, hi = s.split("-")
        return int(lo), int(hi)
    except ValueError:
        print(f"ERROR: {label} must be MIN-MAX (e.g. 8-15), got: {s}", file=sys.stderr)
        sys.exit(1)


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.fasta is None:
        _fatal("Provide a FASTA file or use --demo.")
    if not args.demo and args.fasta and not Path(args.fasta).exists():
        _fatal(f"FASTA not found: {args.fasta}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    pbs_range = _parse_range(args.pbs_range, "--pbs-range")
    rt_range = _parse_range(args.rt_range, "--rt-range")
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            fasta_path=Path(args.fasta) if args.fasta else None,
            edit_json=Path(args.edit) if args.edit else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            pam=args.pam,
            pbs_range=pbs_range,
            rt_range=rt_range,
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
    print(f"\npegRNA Optimiser v{result['pipeline_version']}")
    print(f"Candidates evaluated: {len(result.get('candidates', []))}")
    print(f"Pareto-optimal designs: {len(result.get('pareto_front', []))}")
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
