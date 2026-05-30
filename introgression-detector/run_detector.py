#!/usr/bin/env python3
"""CLI for the Introgression Detector.

Usage:
    python run_detector.py variants.vcf --pop-map populations.tsv --output-dir results/
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
            "Introgression Detector — detects gene flow between populations using "
            "Patterson's D-statistic and f4-ratio in a four-population ABBA-BABA "
            "framework, with sliding-window localisation."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="Multi-population VCF. Omit with --demo.")
    parser.add_argument("--pop-map", metavar="TSV", default=None,
                        help="TSV mapping sample IDs to population labels.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--p1", default="pop1", metavar="POP",
                        help="Population P1 (no introgression expected).")
    parser.add_argument("--p2", default="pop2", metavar="POP",
                        help="Population P2 (no introgression expected).")
    parser.add_argument("--p3", default="pop3", metavar="POP",
                        help="Putative donor population P3.")
    parser.add_argument("--outgroup", default="outgroup", metavar="POP",
                        help="Outgroup population.")
    parser.add_argument("--test",
                        choices=["D_statistic", "f4_ratio", "RND_min", "Dfoil"],
                        default="D_statistic",
                        help="Statistical test.")
    parser.add_argument("--window-kb", type=int, default=50, metavar="INT",
                        help="Sliding window size in kilobases.")
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
    if not args.demo and args.pop_map and not Path(args.pop_map).exists():
        _fatal(f"Population map not found: {args.pop_map}")


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
            population_map=Path(args.pop_map) if args.pop_map else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            p1=args.p1,
            p2=args.p2,
            p3=args.p3,
            outgroup=args.outgroup,
            test=args.test,
            window_kb=args.window_kb,
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
    segments = result.get("segments", [])
    gw = result.get("genome_wide_result")
    print(f"\nIntrogression Detector v{result['pipeline_version']}")
    if gw:
        print(f"Genome-wide D: {gw.get('d_statistic', 'N/A')}")
    print(f"Introgressed segments detected: {len(segments)}")
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
