#!/usr/bin/env python3
"""CLI for the Contig Scaffolding Helper.

Usage:
    python run_helper.py contigs.fa --links reads.bam --output-dir results/
    python run_helper.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_helper.py",
        description=(
            "Contig Scaffolding Helper — orders and orients contigs into scaffolds "
            "using paired-end read links, Hi-C contacts, or reference-guided scaffolding, "
            "joining them with estimated gap sizes."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("contigs", metavar="FASTA", nargs="?", default=None,
                        help="Assembly contigs FASTA. Omit with --demo.")
    parser.add_argument("--links", metavar="BAM/COOL", default=None,
                        help="BAM or Hi-C contact file for scaffolding evidence.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--evidence-type",
                        choices=["paired_reads", "hi_c", "reference_guided"],
                        default="paired_reads",
                        help="Type of scaffolding evidence.")
    parser.add_argument("--min-link-support", type=int, default=3, metavar="INT",
                        help="Minimum links between contigs to join them.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.contigs is None:
        _fatal("Provide a contigs FASTA or use --demo.")
    if args.contigs and not Path(args.contigs).exists():
        _fatal(f"Contigs FASTA not found: {args.contigs}")
    if not args.demo and args.links and not Path(args.links).exists():
        _fatal(f"Links file not found: {args.links}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            contigs_fasta=Path(args.contigs) if args.contigs else None,
            links_file=Path(args.links) if args.links else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            evidence_type=args.evidence_type,
            min_link_support=args.min_link_support,
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
    scaffolds = result.get("scaffolds", [])
    print(f"\nContig Scaffolding Helper v{result['pipeline_version']}")
    print(f"Scaffolds built: {len(scaffolds)}")
    print(f"Contigs placed: {result.get('n_contigs_placed', 'N/A')}")
    print(f"Unplaced contigs: {len(result.get('unplaced_contigs', []))}")
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
