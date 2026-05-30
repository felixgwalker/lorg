#!/usr/bin/env python3
"""CLI for the Promoter Strength Estimator.

Usage:
    python run_estimator.py --annotation genes.gtf --fasta genome.fa --output-dir results/
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
            "Promoter Strength Estimator — scores gene promoters using TATA box strength, "
            "Inr element, CpG island overlap, GC content, TFBS density, and optional "
            "H3K4me3 ChIP-seq signal."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--annotation", metavar="GTF/BED", default=None,
                        help="GTF or BED of gene TSS positions.")
    parser.add_argument("--fasta", metavar="FASTA", default=None,
                        help="Reference genome FASTA.")
    parser.add_argument("--chip-h3k4me3", metavar="BIGWIG", default=None,
                        help="H3K4me3 ChIP-seq bigWig (optional).")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--promoter-window", type=int, default=2000, metavar="INT",
                        help="Upstream window (bp) to define the promoter region.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.annotation is None:
        _fatal("Provide --annotation or use --demo.")
    for p in [args.annotation, args.fasta, args.chip_h3k4me3]:
        if p and not Path(p).exists():
            _fatal(f"File not found: {p}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            gene_annotation=Path(args.annotation) if args.annotation else None,
            fasta=Path(args.fasta) if args.fasta else None,
            chip_h3k4me3=Path(args.chip_h3k4me3) if args.chip_h3k4me3 else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            promoter_window=args.promoter_window,
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
    strengths = result.get("promoter_strengths", [])
    print(f"\nPromoter Strength Estimator v{result['pipeline_version']}")
    print(f"Promoters scored: {len(strengths)}")
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
