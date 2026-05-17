#!/usr/bin/env python3
"""
CLI orchestrator for the Ancient DNA Damage Classifier.

Run from the project root (ancient-dna-damage-classifier/):

    python run_classifier.py reads.bam --output-dir results/

    python run_classifier.py reads.bam \\
        --output-dir results/ \\
        --min-mapq 20 \\
        --min-length 35 \\
        --n-terminal 20 \\
        --prior-ancient 0.95 \\
        --sample-name "mammoth_lib01"

FASTQ input is not supported directly.  Pre-align reads with BWA-MEM and
provide the sorted, indexed BAM:

    bwa mem reference.fa reads.fastq | samtools sort -o reads.bam
    samtools index reads.bam
    python run_classifier.py reads.bam --output-dir results/
"""

import argparse
import logging
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser for the damage classifier CLI."""
    parser = argparse.ArgumentParser(
        prog="run_classifier.py",
        description=(
            "Ancient DNA Damage Classifier — characterises post-mortem DNA damage "
            "patterns and classifies reads as authentically ancient, contaminated, "
            "or ambiguous."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "bam",
        metavar="BAM",
        help=(
            "Path to aligned, coordinate-sorted, indexed BAM file. "
            "FASTQ input is not supported — pre-align with BWA-MEM first."
        ),
    )
    parser.add_argument(
        "-o", "--output-dir",
        metavar="DIR",
        default="./results",
        help="Directory for output files. Created if absent.",
    )
    parser.add_argument(
        "--min-mapq",
        metavar="INT",
        type=int,
        default=30,
        help="Minimum mapping quality (MAPQ); reads below this are excluded.",
    )
    parser.add_argument(
        "--min-length",
        metavar="INT",
        type=int,
        default=30,
        help="Minimum read length in bp; shorter reads are excluded.",
    )
    parser.add_argument(
        "--n-terminal",
        metavar="INT",
        type=int,
        default=25,
        help=(
            "Number of terminal positions to profile at each end. "
            "Must be <= min_length // 2."
        ),
    )
    parser.add_argument(
        "--prior-ancient",
        metavar="FLOAT",
        type=float,
        default=0.9,
        help="Library-level prior probability P(ancient). Must be in (0.0, 1.0).",
    )
    parser.add_argument(
        "--auth-threshold",
        metavar="FLOAT",
        type=float,
        default=0.85,
        help="Posterior P(ancient) >= this value -> classified as 'authentic'.",
    )
    parser.add_argument(
        "--cont-threshold",
        metavar="FLOAT",
        type=float,
        default=0.15,
        help="Posterior P(ancient) <= this value -> classified as 'contaminated'.",
    )
    parser.add_argument(
        "--sample-name",
        metavar="STR",
        default="",
        help="Label used in plot titles and reports. Defaults to BAM filename stem.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip matplotlib figure generation (useful in headless environments).",
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress console summary output. Output files are still written.",
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
    bam_path = Path(args.bam)

    # FASTQ detection — provide a helpful redirect
    fastq_suffixes = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}
    if bam_path.suffix.lower() in fastq_suffixes or ".fastq" in bam_path.name.lower():
        print(
            "ERROR: FASTQ input is not supported.\n"
            "Pre-align reads to a reference genome using BWA-MEM and provide the\n"
            "sorted, indexed BAM file:\n\n"
            "    bwa mem reference.fa reads.fastq | samtools sort -o reads.bam\n"
            "    samtools index reads.bam\n"
            "    python run_classifier.py reads.bam --output-dir results/",
            file=sys.stderr,
        )
        sys.exit(1)

    if not bam_path.exists():
        print(f"ERROR: BAM file not found: {bam_path}", file=sys.stderr)
        sys.exit(1)

    if bam_path.suffix.lower() != ".bam":
        print(
            f"ERROR: Input file does not have a .bam extension: {bam_path}\n"
            "Provide a sorted, indexed BAM file.",
            file=sys.stderr,
        )
        sys.exit(1)

    # BAM index check (.bai or .bam.bai)
    bai_1 = bam_path.with_suffix(".bai")
    bai_2 = Path(str(bam_path) + ".bai")
    if not bai_1.exists() and not bai_2.exists():
        print(
            f"ERROR: BAM index not found for {bam_path}.\n"
            "Run 'samtools index <bam>' before classifying.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.n_terminal > args.min_length // 2:
        print(
            f"ERROR: --n-terminal ({args.n_terminal}) must be <= "
            f"--min-length // 2 ({args.min_length // 2}). "
            "Increase --min-length or decrease --n-terminal to avoid overlapping "
            "the 5\u2032 and 3\u2032 terminal windows on short reads.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.cont_threshold >= args.auth_threshold:
        print(
            f"ERROR: --cont-threshold ({args.cont_threshold}) must be strictly "
            f"less than --auth-threshold ({args.auth_threshold}).",
            file=sys.stderr,
        )
        sys.exit(1)

    if not (0.0 < args.prior_ancient < 1.0):
        print(
            f"ERROR: --prior-ancient ({args.prior_ancient}) must be in the open "
            "interval (0.0, 1.0). Values of exactly 0 or 1 make the Bayesian "
            "update degenerate.",
            file=sys.stderr,
        )
        sys.exit(1)


def print_summary(result: dict) -> None:
    """
    Print a formatted human-readable summary to stdout.

    Args:
        result: The structured dict returned by run_pipeline().
    """
    from src.output_writer import _build_summary_lines

    profile    = result["profile"]
    model      = result["model"]
    summary    = result["summary"]
    frag_stats = result["frag_stats"]

    lines = _build_summary_lines(profile, model, summary, frag_stats)
    print("\n".join(lines))

    out_files = result.get("output_files", {})
    if out_files:
        print("  Output files:")
        labels = {
            "damage_csv":      "  Damage CSV      ->",
            "read_tsv":        "  Read TSV        ->",
            "summary_json":    "  Summary JSON    ->",
            "summary_txt":     "  Summary TXT     ->",
            "damage_plot_png": "  Plot PNG        ->",
            "damage_plot_svg": "  Plot SVG        ->",
        }
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
            bam_path=Path(args.bam),
            output_dir=Path(args.output_dir),
            min_mapq=args.min_mapq,
            min_length=args.min_length,
            n_terminal=args.n_terminal,
            prior_ancient=args.prior_ancient,
            auth_threshold=args.auth_threshold,
            cont_threshold=args.cont_threshold,
            sample_name=args.sample_name,
            no_plot=args.no_plot,
        )
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: Unexpected failure — {exc}", file=sys.stderr)
        logging.exception("Unhandled exception in pipeline")
        return 1

    if not args.quiet:
        print_summary(result)

    return 0


if __name__ == "__main__":
    sys.exit(main())
