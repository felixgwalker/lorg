#!/usr/bin/env python3
"""CLI for the Lineage Divergence Dater.

Usage:
    python run_dater.py alignment.fa --phylogeny tree.nwk --calibrations calibrations.json --output-dir results/
    python run_dater.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_dater.py",
        description=(
            "Lineage Divergence Dater — estimates split times between lineages using a "
            "molecular clock with calibrations, Bayesian dated-tips for ancient DNA, "
            "or simple pairwise genetic distance."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("alignment", metavar="FASTA", nargs="?", default=None,
                        help="FASTA multiple sequence alignment. Omit with --demo.")
    parser.add_argument("--phylogeny", metavar="NWK", default=None,
                        help="Newick phylogenetic tree.")
    parser.add_argument("--calibrations", metavar="JSON", default=None,
                        help="JSON of calibration points or radiocarbon-dated tip ages.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--method",
                        choices=["molecular_clock", "bayesian_dated_tips",
                                 "pairwise_distance"],
                        default="molecular_clock",
                        help="Dating method.")
    parser.add_argument("--mutation-rate", type=float, default=1.25e-8, metavar="FLOAT",
                        help="Per-base per-generation mutation rate.")
    parser.add_argument("--generation-time", type=float, default=30.0, metavar="YEARS",
                        help="Generation time in years.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.alignment is None:
        _fatal("Provide an alignment FASTA or use --demo.")
    if args.alignment and not Path(args.alignment).exists():
        _fatal(f"Alignment not found: {args.alignment}")
    if args.calibrations and not Path(args.calibrations).exists():
        _fatal(f"Calibrations not found: {args.calibrations}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            alignment=Path(args.alignment) if args.alignment else None,
            phylogeny=Path(args.phylogeny) if args.phylogeny else None,
            calibrations=Path(args.calibrations) if args.calibrations else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            method=args.method,
            mutation_rate=args.mutation_rate,
            generation_time=args.generation_time,
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
    dates = result.get("divergence_dates", [])
    print(f"\nLineage Divergence Dater v{result['pipeline_version']}")
    print(f"Divergence dates estimated: {len(dates)}")
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
