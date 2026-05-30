#!/usr/bin/env python3
"""CLI for the Rare Variant Prioritiser.

Usage:
    python run_prioritiser.py variants.vcf --output-dir results/
    python run_prioritiser.py --demo --hpo-terms HP:0001250 HP:0000256 --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_prioritiser.py",
        description=(
            "Rare Variant Prioritiser — combines allele frequency, CADD score, gene "
            "constraint, HPO phenotype matching, and ClinVar evidence into a weighted "
            "composite score for three-tier variant prioritisation."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="Annotated VCF with gnomAD AF and CADD scores. Omit with --demo.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--max-af", type=float, default=0.01, metavar="FLOAT",
                        help="Maximum gnomAD allele frequency to consider a variant rare.")
    parser.add_argument("--min-cadd", type=float, default=20.0, metavar="FLOAT",
                        help="Minimum CADD Phred score to include.")
    parser.add_argument("--hpo-terms", nargs="*", metavar="HPO_ID", default=None,
                        help="HPO term IDs for phenotype matching (e.g. HP:0001250).")
    parser.add_argument("--gene-panel", metavar="FILE", default=None,
                        help="Text file of gene symbols (one per line).")
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
    if args.gene_panel and not Path(args.gene_panel).exists():
        _fatal(f"Gene panel file not found: {args.gene_panel}")


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
            output_dir=Path(args.output_dir),
            demo=args.demo,
            max_af=args.max_af,
            min_cadd=args.min_cadd,
            hpo_terms=args.hpo_terms,
            gene_panel=Path(args.gene_panel) if args.gene_panel else None,
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
    variants = result.get("prioritised_variants", [])
    print(f"\nRare Variant Prioritiser v{result['pipeline_version']}")
    print(f"Variants prioritised: {len(variants)}")
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
