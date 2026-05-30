#!/usr/bin/env python3
"""CLI for the Repair Pathway Bias Estimator.

Usage:
    python run_estimator.py cutsite.fa --cell-type iPSC --output-dir results/
    python run_estimator.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_estimator.py",
        description=(
            "Repair Pathway Bias Estimator — estimates NHEJ, MMEJ, and HDR "
            "pathway probabilities at a CRISPR cut site using microhomology "
            "enumeration and cell-type HDR bias tables."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("fasta", metavar="FASTA", nargs="?", default=None,
                        help="FASTA with ≥50 bp flanking each side of cut. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Use a synthetic 120 bp flanking sequence.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--cell-type", default="HEK293", metavar="TYPE",
                        help="Cell type for HDR bias lookup (e.g. HEK293, iPSC, neuron).")
    parser.add_argument("--cut-position", type=int, default=None, metavar="INT",
                        help="1-based cut position in FASTA (default: centre).")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip pie chart generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


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
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            fasta_path=Path(args.fasta) if args.fasta else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            cell_type=args.cell_type,
            cut_position=args.cut_position,
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
    probs = result.get("pathway_probabilities", {})
    print(f"\nRepair Pathway Bias Estimator v{result['pipeline_version']}")
    for k, v in probs.items():
        print(f"  {k}: {v:.3f}")
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
