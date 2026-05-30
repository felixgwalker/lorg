#!/usr/bin/env python3
"""CLI for the Structural Variant Prioritiser.

Usage:
    python run_prioritiser.py svs.vcf --gene-annotation genes.bed --dosage-scores clingen_hi.tsv --output-dir results/
    python run_prioritiser.py --demo --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_prioritiser.py",
        description=(
            "Structural Variant Prioritiser — prioritises SVs by gene overlap, "
            "ClinGen dosage sensitivity (HI/TS), population frequency, and "
            "overlap with DECIPHER and ClinVar pathogenic SVs."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("vcf", metavar="VCF", nargs="?", default=None,
                        help="VCF or BED of structural variants. Omit with --demo.")
    parser.add_argument("--gene-annotation", metavar="BED/GTF", default=None,
                        help="Gene annotation BED or GTF.")
    parser.add_argument("--dosage-scores", metavar="TSV", default=None,
                        help="ClinGen dosage sensitivity scores TSV.")
    parser.add_argument("--demo", action="store_true",
                        help="Run on synthetic data without real input.")
    parser.add_argument("-o", "--output-dir", metavar="DIR", default="./results",
                        help="Output directory.")
    parser.add_argument("--sv-types", nargs="*", metavar="TYPE",
                        default=None,
                        help="SV types to include (e.g. DEL DUP INV).")
    parser.add_argument("--min-size", type=int, default=50, metavar="INT",
                        help="Minimum SV size in base pairs.")
    parser.add_argument("--max-af", type=float, default=0.01, metavar="FLOAT",
                        help="Maximum population AF to call an SV rare.")
    parser.add_argument("--haploinsufficiency-threshold", type=float, default=0.9,
                        metavar="FLOAT",
                        help="ClinGen HI score threshold for dosage-sensitive genes.")
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
    if not args.demo and args.gene_annotation and not Path(args.gene_annotation).exists():
        _fatal(f"Gene annotation not found: {args.gene_annotation}")
    if not args.demo and args.dosage_scores and not Path(args.dosage_scores).exists():
        _fatal(f"Dosage scores file not found: {args.dosage_scores}")


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
            gene_annotation=Path(args.gene_annotation) if args.gene_annotation else None,
            dosage_scores=Path(args.dosage_scores) if args.dosage_scores else None,
            output_dir=Path(args.output_dir),
            demo=args.demo,
            sv_types=args.sv_types,
            min_size=args.min_size,
            max_af=args.max_af,
            haploinsufficiency_threshold=args.haploinsufficiency_threshold,
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
    svs = result.get("prioritised_svs", [])
    print(f"\nStructural Variant Prioritiser v{result['pipeline_version']}")
    print(f"SVs prioritised: {len(svs)}")
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
