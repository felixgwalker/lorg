#!/usr/bin/env python3
"""CLI for the Promoter Variant Scorer.

Usage:
    python run_scorer.py variants.vcf --fasta genome.fa --pwm-db jaspar.meme --output-dir results/
    python run_scorer.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_scorer.py",
        description=(
            "Promoter Variant Scorer — scores variants in promoter regions against a "
            "JASPAR-format PWM database to identify transcription factor binding site "
            "disruptions and de novo TFBS creations."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="VCF of promoter variants. Omit with --demo.")
    parser.add_argument("--fasta", metavar="FASTA", default=None,
                        help="Reference genome FASTA.")
    parser.add_argument("--pwm-db", metavar="MEME/TRANSFAC", default=None,
                        help="JASPAR-format PWM database.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--promoter-window", type=int, default=2000, metavar="INT",
                        help="Upstream window (bp) from TSS to define the promoter region.")
    parser.add_argument("--ic-threshold", type=float, default=8.0, metavar="FLOAT",
                        help="Minimum information content (bits) to include a PWM.")
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
    if not args.demo and args.pwm_db and not Path(args.pwm_db).exists():
        _fatal(f"PWM database not found: {args.pwm_db}")


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
            pwm_db=Path(args.pwm_db) if args.pwm_db else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            promoter_window=args.promoter_window,
            ic_threshold=args.ic_threshold,
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
    disruptions = result.get("disruptions", [])
    print(f"\nPromoter Variant Scorer v{result['pipeline_version']}")
    print(f"TFBS events identified: {len(disruptions)}")
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
