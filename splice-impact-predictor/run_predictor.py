#!/usr/bin/env python3
"""CLI for the Splice Impact Predictor.

Usage:
    python run_predictor.py variants.vcf --fasta genome.fa --gtf genes.gtf --output-dir results/
    python run_predictor.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_predictor.py",
        description=(
            "Splice Impact Predictor — predicts the effect of variants on splicing by "
            "scoring donor and acceptor sites with position weight matrices and reporting "
            "delta scores between reference and alternate alleles."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="VCF of variants to evaluate. Omit with --demo.")
    parser.add_argument("--fasta", metavar="FASTA", default=None,
                        help="Reference genome FASTA.")
    parser.add_argument("--gtf", metavar="GTF", default=None,
                        help="Gene annotation GTF.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--window", type=int, default=200, metavar="INT",
                        help="Bases around each variant to search for splice sites.")
    parser.add_argument("--delta-threshold", type=float, default=2.0, metavar="FLOAT",
                        help="Absolute delta score to call a significant splicing effect.")
    parser.add_argument("--canonical-only", action="store_true",
                        help="Restrict analysis to canonical GT-AG splice sites.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.vcf is None:
        _fatal("Provide a VCF file or use --demo.")
    if not args.demo and args.vcf and not Path(args.vcf).exists():
        _fatal(f"VCF not found: {args.vcf}")
    if not args.demo and args.fasta and not Path(args.fasta).exists():
        _fatal(f"FASTA not found: {args.fasta}")
    if not args.demo and args.gtf and not Path(args.gtf).exists():
        _fatal(f"GTF not found: {args.gtf}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(levelname)s: %(message)s", stream=sys.stderr)
    validate_args(args)
    from src.pipeline import run_pipeline
    try:
        result = run_pipeline(
            vcf_path=Path(args.vcf) if args.vcf else None,
            fasta_path=Path(args.fasta) if args.fasta else None,
            gtf_path=Path(args.gtf) if args.gtf else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            window=args.window,
            delta_threshold=args.delta_threshold,
            canonical_only=args.canonical_only,
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
    scores = result.get("splice_scores", [])
    print(f"\nSplice Impact Predictor v{result['pipeline_version']}")
    print(f"Variants scored: {len(scores)}")
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
