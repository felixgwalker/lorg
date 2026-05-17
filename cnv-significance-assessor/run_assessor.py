#!/usr/bin/env python3
"""
CLI orchestrator for the CNV Significance Assessor.

Run from the project root (cnv-significance-assessor/):

    python run_assessor.py variants.bed annotation.gff3 --output-dir results/

    python run_assessor.py variants.vcf annotation.gff3 \\
        --output-dir results/ \\
        --pop-db gnomad_sv.bed \\
        --dosage-scores pHaplo_pTriplo.csv \\
        --min-cnv-size 5000 \\
        --overlap-fraction 0.2 \\
        --pop-freq-cutoff 0.005 \\
        --sample-name "mammoth_sample_01"

Both BED and VCF CNV inputs are accepted; format is auto-detected.
The GFF3 annotation must use standard column layout (9 tab-separated fields).
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser for the CNV assessor CLI."""
    parser = argparse.ArgumentParser(
        prog="run_assessor.py",
        description=(
            "CNV Significance Assessor — annotates copy number variants with "
            "gene content, dosage sensitivity, and population frequency, then "
            "classifies each as LIKELY_BENIGN, VUS, or LIKELY_PATHOGENIC."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # ── Required positional arguments ────────────────────────────────────
    parser.add_argument(
        "cnv_input",
        metavar="CNV",
        help=(
            "CNV calls in BED or VCF format.  Format is auto-detected from "
            "the file extension or header.  "
            "BED: chrom, start, end, name[, score, strand, cnv_type].  "
            "VCF: SVTYPE in INFO field required."
        ),
    )
    parser.add_argument(
        "gff3",
        metavar="GFF3",
        help="Reference genome gene annotation in GFF3 format.",
    )

    # ── Output ───────────────────────────────────────────────────────────
    parser.add_argument(
        "-o", "--output-dir",
        metavar="DIR",
        default="./results",
        help="Directory for output files.  Created if absent.",
    )

    # ── Optional inputs ───────────────────────────────────────────────────
    parser.add_argument(
        "--pop-db",
        metavar="FILE",
        default=None,
        help=(
            "Population-level CNV frequency database in BED or VCF format "
            "(e.g. DGV, gnomAD-SV, or a custom panel).  When provided, "
            "CNVs with matching population records above --pop-freq-cutoff "
            "are classified as LIKELY_BENIGN."
        ),
    )
    parser.add_argument(
        "--dosage-scores",
        metavar="CSV",
        default=None,
        help=(
            "Haploinsufficiency and triplosensitivity scores in CSV format.  "
            "Expected columns (case-insensitive): gene/symbol, and any of "
            "pLI, pHaplo, pTriplo.  "
            "Compatible with gnomAD and ClinGen score files."
        ),
    )

    # ── Filter parameters ─────────────────────────────────────────────────
    parser.add_argument(
        "--min-cnv-size",
        metavar="INT",
        type=int,
        default=1000,
        help="Minimum CNV size in bp; smaller variants are excluded.",
    )
    parser.add_argument(
        "--overlap-fraction",
        metavar="FLOAT",
        type=float,
        default=0.1,
        help=(
            "Minimum reciprocal overlap fraction between a CNV and a gene "
            "for the gene to be counted as overlapping.  "
            "Applied to the shorter of the two intervals."
        ),
    )
    parser.add_argument(
        "--pop-freq-cutoff",
        metavar="FLOAT",
        type=float,
        default=0.01,
        help=(
            "Population frequency threshold.  CNVs with frequency above this "
            "value are classified as LIKELY_BENIGN regardless of other scores."
        ),
    )

    # ── Misc ──────────────────────────────────────────────────────────────
    parser.add_argument(
        "--sample-name",
        metavar="STR",
        default="",
        help="Label used in plot titles and summary headers.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip matplotlib figure generation (useful in headless environments).",
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress console summary output.  Output files are still written.",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging verbosity.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 1.0.0",
    )

    return parser


def validate_args(args: argparse.Namespace) -> None:
    """
    Validate argument values before pipeline execution.

    Raises SystemExit(1) with a descriptive message on any violation.
    """
    cnv_path = Path(args.cnv_input)
    if not cnv_path.exists():
        _fatal(f"CNV input file not found: {cnv_path}")

    gff_path = Path(args.gff3)
    if not gff_path.exists():
        _fatal(f"GFF3 annotation file not found: {gff_path}")

    if args.pop_db is not None and not Path(args.pop_db).exists():
        _fatal(f"Population database file not found: {args.pop_db}")

    if args.dosage_scores is not None and not Path(args.dosage_scores).exists():
        _fatal(f"Dosage score file not found: {args.dosage_scores}")

    if args.min_cnv_size < 0:
        _fatal(f"--min-cnv-size must be ≥ 0, got {args.min_cnv_size}.")

    if not (0.0 < args.overlap_fraction <= 1.0):
        _fatal(
            f"--overlap-fraction must be in (0.0, 1.0], got {args.overlap_fraction}."
        )

    if not (0.0 < args.pop_freq_cutoff < 1.0):
        _fatal(
            f"--pop-freq-cutoff must be in (0.0, 1.0), got {args.pop_freq_cutoff}."
        )


def print_summary(result: dict) -> None:
    """
    Print a formatted human-readable summary to stdout.

    Args:
        result: The structured dict returned by run_pipeline().
    """
    from src.output_writer import _build_summary_lines

    summary = result["summary"]
    scored  = result["scored"]
    lines = _build_summary_lines(summary, scored)
    print("\n".join(lines))

    out_files = result.get("output_files", {})
    if out_files:
        labels = {
            "annotated_csv":     "  Annotated CSV    ->",
            "gene_impact_csv":   "  Gene Impact CSV  ->",
            "significance_txt":  "  Summary TXT      ->",
            "significance_json": "  Summary JSON     ->",
            "ideogram_png":      "  Ideogram PNG     ->",
            "ideogram_svg":      "  Ideogram SVG     ->",
        }
        print("  Output files:")
        for key, label in labels.items():
            if key in out_files:
                print(f"  {label} {out_files[key]}")
        print()


def main(argv: list[str] | None = None) -> int:
    """
    CLI entry point.

    Args:
        argv: Argument list (defaults to sys.argv[1:] when None).

    Returns:
        Integer exit code (0 = success, 1 = error).
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    # Configure logging
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s: %(message)s",
        stream=sys.stderr,
    )

    validate_args(args)

    # Import here so logging is configured first
    from src.pipeline import run_pipeline

    try:
        result = run_pipeline(
            cnv_path=Path(args.cnv_input),
            gff_path=Path(args.gff3),
            output_dir=Path(args.output_dir),
            pop_db_path=Path(args.pop_db) if args.pop_db else None,
            dosage_path=Path(args.dosage_scores) if args.dosage_scores else None,
            min_cnv_size=args.min_cnv_size,
            overlap_fraction=args.overlap_fraction,
            pop_freq_cutoff=args.pop_freq_cutoff,
            sample_name=args.sample_name,
            no_plot=args.no_plot,
        )
    except FileNotFoundError as exc:
        _fatal(str(exc))
    except ValueError as exc:
        _fatal(str(exc))
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: Unexpected failure — {exc}", file=sys.stderr)
        logging.exception("Unhandled exception in pipeline")
        return 1

    if not args.quiet:
        print_summary(result)

    return 0


def _fatal(msg: str) -> None:
    """Print an error message to stderr and exit with code 1."""
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
