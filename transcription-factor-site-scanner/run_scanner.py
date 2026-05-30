#!/usr/bin/env python3
"""CLI for the Transcription Factor Site Scanner.

Usage:
    python run_scanner.py sequences.fa --pwm-db jaspar.meme --output-dir results/
    python run_scanner.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_scanner.py",
        description=(
            "Transcription Factor Site Scanner — scans FASTA sequences against a JASPAR "
            "or MEME-format PWM database, reporting TFBS hits with p-values against a "
            "background nucleotide model."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("sequences", metavar="FASTA", nargs="?", default=None,
                        help="FASTA of sequences to scan. Omit with --demo.")
    parser.add_argument("--pwm-db", metavar="MEME/JASPAR", default=None,
                        help="JASPAR or MEME-format PWM database.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--p-value-threshold", type=float, default=1e-4, metavar="FLOAT",
                        help="P-value cutoff for reporting TFBS hits.")
    parser.add_argument("--min-ic", type=float, default=8.0, metavar="FLOAT",
                        help="Minimum PWM information content (bits) to include.")
    parser.add_argument("--single-strand", action="store_true",
                        help="Scan forward strand only.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.sequences is None:
        _fatal("Provide a FASTA file or use --demo.")
    if args.sequences and not Path(args.sequences).exists():
        _fatal(f"FASTA not found: {args.sequences}")
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
            sequences=Path(args.sequences) if args.sequences else None,
            pwm_db=Path(args.pwm_db) if args.pwm_db else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            p_value_threshold=args.p_value_threshold,
            min_ic=args.min_ic,
            both_strands=not args.single_strand,
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
    hits = result.get("hits", [])
    print(f"\nTranscription Factor Site Scanner v{result['pipeline_version']}")
    print(f"Sequences scanned: {result.get('n_sequences_scanned', 'N/A')}")
    print(f"PWMs used: {result.get('n_pwms_used', 'N/A')}")
    print(f"TFBS hits: {len(hits)}")
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
