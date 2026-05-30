#!/usr/bin/env python3
"""CLI for the Prime Edit Design Assistant.

Usage:
    python run_assistant.py target.fa --edit edit.json --output-dir results/
    python run_assistant.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_assistant.py",
        description=(
            "Prime Edit Design Assistant — designs pegRNAs for prime editing by "
            "scanning a target locus for PAM sites and enumerating RT template / "
            "PBS combinations, ranked by predicted activity features."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("fasta", metavar="FASTA", nargs="?", default=None,
                        help="Target locus FASTA (≥200 bp). Omit with --demo.")
    parser.add_argument("--edit", metavar="JSON", default=None,
                        help="Edit spec JSON: {position, ref, alt}.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--pam", default="NGG", metavar="PAM",
                        help="PAM motif (e.g. NGG, NNGRRT, TTTV).")
    parser.add_argument("--pbs-min", type=int, default=8, metavar="INT",
                        help="Minimum PBS length (nt).")
    parser.add_argument("--pbs-max", type=int, default=15, metavar="INT",
                        help="Maximum PBS length (nt).")
    parser.add_argument("--rt-min", type=int, default=10, metavar="INT",
                        help="Minimum RT template length (nt).")
    parser.add_argument("--rt-max", type=int, default=16, metavar="INT",
                        help="Maximum RT template length (nt).")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.fasta is None:
        _fatal("Provide a FASTA file or use --demo.")
    if not args.demo and args.fasta and not Path(args.fasta).exists():
        _fatal(f"FASTA not found: {args.fasta}")
    if args.edit and not Path(args.edit).exists():
        _fatal(f"Edit JSON not found: {args.edit}")
    if args.pbs_min >= args.pbs_max:
        _fatal("--pbs-min must be less than --pbs-max.")
    if args.rt_min >= args.rt_max:
        _fatal("--rt-min must be less than --rt-max.")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            fasta_path=Path(args.fasta) if args.fasta else None,
            edit_json=Path(args.edit) if args.edit else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            pam=args.pam,
            pbs_min=args.pbs_min,
            pbs_max=args.pbs_max,
            rt_min=args.rt_min,
            rt_max=args.rt_max,
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
    designs = result.get("designs", [])
    print(f"\nPrime Edit Design Assistant v{result['pipeline_version']}")
    print(f"pegRNA designs generated: {len(designs)}")
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
