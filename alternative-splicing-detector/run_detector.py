#!/usr/bin/env python3
"""CLI for the Alternative Splicing Detector.

Usage:
    python run_detector.py --bam-a condA.bam --bam-b condB.bam --annotation genes.gtf --output-dir results/
    python run_detector.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_detector.py",
        description=(
            "Alternative Splicing Detector — detects differential exon skipping, "
            "intron retention, and alternative 5'/3' splice site usage between two "
            "RNA-seq conditions using PSI quantification."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bam-a", metavar="BAM", default=None,
                        help="RNA-seq BAM for condition A.")
    parser.add_argument("--bam-b", metavar="BAM", default=None,
                        help="RNA-seq BAM for condition B.")
    parser.add_argument("--annotation", metavar="GTF", default=None,
                        help="GTF gene annotation.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--event-types", nargs="*",
                        choices=["exon_skipping", "intron_retention", "alt_5_prime",
                                 "alt_3_prime", "mutually_exclusive",
                                 "alt_first_exon", "alt_last_exon"],
                        default=None, help="Event types to detect (default: all).")
    parser.add_argument("--min-delta-psi", type=float, default=0.1, metavar="FLOAT",
                        help="Minimum |ΔPSI| to report as differential.")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip figure generation.")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        default="INFO", help="Logging verbosity.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and args.bam_a is None:
        _fatal("Provide --bam-a and --bam-b, or use --demo.")
    for p in [args.bam_a, args.bam_b, args.annotation]:
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
            bam_a=Path(args.bam_a) if args.bam_a else None,
            bam_b=Path(args.bam_b) if args.bam_b else None,
            annotation=Path(args.annotation) if args.annotation else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            event_types=args.event_types,
            min_delta_psi=args.min_delta_psi,
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
    events = result.get("events", [])
    print(f"\nAlternative Splicing Detector v{result['pipeline_version']}")
    print(f"Splicing events detected: {len(events)}")
    print(f"Genes with events: {result.get('n_genes_with_events', 'N/A')}")
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
