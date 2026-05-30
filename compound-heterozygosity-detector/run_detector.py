#!/usr/bin/env python3
"""CLI for the Compound Heterozygosity Detector.

Usage:
    python run_detector.py variants.vcf --ped family.ped --output-dir results/
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
            "Compound Heterozygosity Detector — identifies pairs of heterozygous "
            "variants in the same gene on opposite haplotypes, using phased genotypes, "
            "trio phase-by-transmission, or statistical phasing."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="Phased or unphased VCF. Omit with --demo.")
    parser.add_argument("--ped", metavar="PED", default=None,
                        help="PLINK PED file describing family structure.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--max-af", type=float, default=0.01, metavar="FLOAT",
                        help="Maximum gnomAD AF to consider a variant as a comp-het candidate.")
    parser.add_argument("--sample-id", metavar="ID", default=None,
                        help="Proband sample ID (required for multi-sample VCF without PED).")
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
    if args.ped and not Path(args.ped).exists():
        _fatal(f"PED file not found: {args.ped}")


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
            ped_path=Path(args.ped) if args.ped else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            max_af=args.max_af,
            sample_id=args.sample_id,
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
    pairs = result.get("comp_het_pairs", [])
    print(f"\nCompound Heterozygosity Detector v{result['pipeline_version']}")
    print(f"Comp-het pairs detected: {len(pairs)}")
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
