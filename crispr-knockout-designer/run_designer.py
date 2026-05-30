#!/usr/bin/env python3
"""CLI for the CRISPR Knockout Designer.

Usage:
    python run_designer.py gene.fa --output-dir results/
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
            "CRISPR Knockout Designer — designs sgRNAs targeting early coding "
            "exons for gene disruption, scored by Doench 2016-style on-target "
            "features and predicted frameshift efficiency."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("fasta", metavar="FASTA", nargs="?", default=None,
                        help="Gene or CDS FASTA. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic gene sequence.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--pam", default="NGG", metavar="PAM",
                        help="PAM motif.")
    parser.add_argument("--n-guides", type=int, default=5, metavar="INT",
                        help="Number of top guides to return.")
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
    if args.n_guides < 1:
        _fatal("--n-guides must be >= 1.")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            gene_fasta=Path(args.fasta) if args.fasta else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            pam=args.pam,
            n_guides=args.n_guides,
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
    guides = result.get("guides", [])
    print(f"\nCRISPR Knockout Designer v{result['pipeline_version']}")
    print(f"Guides designed: {len(guides)}")
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
