"""Pipeline orchestrator for erv-risk-mapper."""

import sys
from pathlib import Path


def run_pipeline(
    genome: dict[str, str],
    output_dir: Path,
    no_plot: bool = False,
    plot_fmt: str = "png",
) -> dict:
    from .erv_detector import detect_erv_windows
    from .age_estimator import estimate_age
    from .risk_scorer import score_all
    from .report import write_bed, write_csv_report
    from .plot import plot_chromosome_density

    output_dir.mkdir(parents=True, exist_ok=True)

    hits = detect_erv_windows(genome)

    for i, hit in enumerate(hits):
        hit["age_mya"] = estimate_age(hit, seed_offset=i)

    hits = score_all(hits)

    bed_path = write_bed(hits, output_dir)
    csv_path = write_csv_report(hits, output_dir)

    plot_path = None
    if not no_plot:
        plot_path = plot_chromosome_density(hits, output_dir, fmt=plot_fmt)

    tier_counts = {"low": 0, "moderate": 0, "high": 0}
    for h in hits:
        tier_counts[h["risk_tier"]] = tier_counts.get(h["risk_tier"], 0) + 1

    return {
        "n_erv_hits": len(hits),
        "tier_low": tier_counts["low"],
        "tier_moderate": tier_counts["moderate"],
        "tier_high": tier_counts["high"],
        "bed": str(bed_path),
        "csv": str(csv_path),
        "plot": str(plot_path) if plot_path else None,
    }


def main() -> int:
    import argparse
    from .fasta_parser import parse_fasta, make_demo_genome

    parser = argparse.ArgumentParser(
        prog="erv-risk-mapper",
        description="Identify ERV elements and risk-score them for reactivation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--genome", type=Path, metavar="FASTA",
                        help="Input genome FASTA file")
    parser.add_argument("--output-dir", type=Path, default=Path("results"),
                        metavar="PATH", help="Output directory")
    parser.add_argument("--no-plot", action="store_true", help="Skip chromosome density plot")
    parser.add_argument("--plot-format", choices=["png", "svg"], default="png",
                        help="Plot file format")
    parser.add_argument("--demo", action="store_true",
                        help="Run demo: 100kb synthetic genome with 8 planted ERV elements")
    args = parser.parse_args()

    if not args.demo and args.genome is None:
        parser.error("--genome is required unless --demo is set")

    if args.demo:
        genome = make_demo_genome()
        total_bp = sum(len(v) for v in genome.values())
        print(f"Demo mode: {len(genome)} chromosomes, {total_bp:,} bp total")
    else:
        genome = parse_fasta(args.genome)
        total_bp = sum(len(v) for v in genome.values())
        print(f"Genome: {len(genome)} sequences, {total_bp:,} bp total")

    result = run_pipeline(
        genome=genome,
        output_dir=args.output_dir,
        no_plot=args.no_plot,
        plot_fmt=args.plot_format,
    )

    print(f"ERV elements detected: {result['n_erv_hits']}")
    print(f"  High risk:     {result['tier_high']}")
    print(f"  Moderate risk: {result['tier_moderate']}")
    print(f"  Low risk:      {result['tier_low']}")
    print(f"  BED file:      {result['bed']}")
    print(f"  CSV report:    {result['csv']}")
    if result["plot"]:
        print(f"  Density plot:  {result['plot']}")
    return 0
