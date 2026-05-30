#!/usr/bin/env python3
"""CLI for the Safe Harbour Integration Finder.

Usage:
    python run_finder.py genome.fa --species human --output-dir results/
    python run_finder.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_finder.py",
        description=(
            "Safe Harbour Integration Finder — identifies genomic safe harbour "
            "sites by searching for AAVS1/H11/Rosa26 homologs and scoring "
            "intergenic regions by oncogene distance and regulatory context."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("genome", metavar="GENOME", nargs="?", default=None,
                        help="Genome FASTA (or chromosomal subset). Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Use a synthetic 10 Mb chromosome fragment.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--species", default="human", choices=["human", "mouse"],
                        help="Species for known-harbour lookup.")
    parser.add_argument("--regulatory-bed", default=None, metavar="BED",
                        help="BED of regulatory elements to avoid.")
    parser.add_argument("--oncogene-bed", default=None, metavar="BED",
                        help="BED of oncogene coordinates.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.genome is None:
        _fatal("Provide a genome FASTA or use --demo.")
    if not args.demo and args.genome and not Path(args.genome).exists():
        _fatal(f"Genome FASTA not found: {args.genome}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            genome_fasta=Path(args.genome) if args.genome else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            regulatory_bed=Path(args.regulatory_bed) if args.regulatory_bed else None,
            oncogene_bed=Path(args.oncogene_bed) if args.oncogene_bed else None,
            species=args.species,
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
    candidates = result.get("candidates", [])
    print(f"\nSafe Harbour Integration Finder v{result['pipeline_version']}")
    print(f"Candidates identified: {len(candidates)}")
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
