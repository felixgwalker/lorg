#!/usr/bin/env python3
"""CLI for the Post-Mortem Damage Simulator.

Usage:
    python run_simulator.py reference.fa --n-reads 100000 --output-dir results/
    python run_simulator.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_simulator.py",
        description=(
            "Post-Mortem Damage Simulator — generates synthetic ancient DNA reads with "
            "realistic C→T / G→A deamination patterns, fragment length distributions, "
            "and nick frequencies using the Briggs or related damage model."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("reference", metavar="FASTA", nargs="?", default=None,
                        help="Reference genome FASTA. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--model", choices=["briggs", "uniform", "double_stranded",
                                             "single_stranded"],
                        default="briggs", help="Post-mortem damage model.")
    parser.add_argument("--n-reads", type=int, default=100000, metavar="INT",
                        help="Number of reads to simulate.")
    parser.add_argument("--mean-fragment-length", type=int, default=60, metavar="INT",
                        help="Mean fragment length in base pairs.")
    parser.add_argument("--deamination-rate-ss", type=float, default=0.68, metavar="FLOAT",
                        help="Single-stranded deamination rate (Briggs model).")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.reference is None:
        _fatal("Provide a reference FASTA or use --demo.")
    if args.reference and not Path(args.reference).exists():
        _fatal(f"Reference FASTA not found: {args.reference}")
    if args.mean_fragment_length <= 0:
        _fatal("--mean-fragment-length must be positive.")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            reference_fasta=Path(args.reference) if args.reference else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            model=args.model,
            n_reads=args.n_reads,
            mean_fragment_length=args.mean_fragment_length,
            deamination_rate_ss=args.deamination_rate_ss,
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
    reads = result.get("simulated_reads", [])
    print(f"\nPost-Mortem Damage Simulator v{result['pipeline_version']}")
    print(f"Reads simulated: {len(reads)}")
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
