"""Pipeline orchestrator for guide-rna-off-target-scorer."""

import sys
from pathlib import Path


def run_pipeline(
    guides: list[tuple[str, str]],
    genome: dict[str, str],
    output_dir: Path,
    no_plot: bool = False,
    plot_fmt: str = "png",
) -> dict:
    from .off_target_finder import find_off_targets
    from .cfd_scorer import score_hits
    from .report import write_off_target_table, write_specificity_summary
    from .plot import plot_manhattan

    output_dir.mkdir(parents=True, exist_ok=True)

    all_hits = []
    for guide_name, guide_seq in guides:
        hits = find_off_targets(guide_name, guide_seq, genome)
        hits = score_hits(hits)
        all_hits.extend(hits)

    ot_path = write_off_target_table(all_hits, output_dir)
    spec_path = write_specificity_summary(all_hits, guides, output_dir)

    plot_path = None
    if not no_plot:
        plot_path = plot_manhattan(all_hits, output_dir, fmt=plot_fmt)

    return {
        "n_guides": len(guides),
        "n_hits": len(all_hits),
        "off_target_table": str(ot_path),
        "specificity_summary": str(spec_path),
        "plot": str(plot_path) if plot_path else None,
    }


def main() -> int:
    import argparse
    from .guide_parser import parse_guides, validate_guide, make_demo_guides
    from .genome_indexer import load_fasta_genome, make_demo_genome

    parser = argparse.ArgumentParser(
        prog="guide-rna-off-target-scorer",
        description="Predict CRISPR off-target cleavage sites using CFD scoring.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--guides", type=Path, metavar="FILE",
                        help="Guide RNA sequences (FASTA or CSV)")
    parser.add_argument("--genome", type=Path, metavar="FASTA",
                        help="Genome FASTA file")
    parser.add_argument("--output-dir", type=Path, default=Path("results"),
                        metavar="PATH", help="Output directory")
    parser.add_argument("--no-plot", action="store_true", help="Skip plot generation")
    parser.add_argument("--plot-format", choices=["png", "svg"], default="png",
                        help="Plot output format")
    parser.add_argument("--demo", action="store_true",
                        help="Run with synthetic data (no input files required)")
    args = parser.parse_args()

    if not args.demo and (args.guides is None or args.genome is None):
        parser.error("--guides and --genome are required unless --demo is set")

    if args.demo:
        guides = make_demo_guides()
        genome = make_demo_genome()
        print("Running in demo mode: synthetic genome (5000bp) + 3 guide RNAs")
    else:
        raw_guides = parse_guides(args.guides)
        guides = [(name, validate_guide(seq)) for name, seq in raw_guides]
        genome = load_fasta_genome(args.genome)
        print(f"Loaded {len(guides)} guides, genome {sum(len(v) for v in genome.values())} bp")

    result = run_pipeline(
        guides=guides,
        genome=genome,
        output_dir=args.output_dir,
        no_plot=args.no_plot,
        plot_fmt=args.plot_format,
    )

    print(f"Guides processed: {result['n_guides']}")
    print(f"Total hits (<=3 mismatches): {result['n_hits']}")
    print(f"Off-target table:      {result['off_target_table']}")
    print(f"Specificity summary:   {result['specificity_summary']}")
    if result["plot"]:
        print(f"Manhattan plot:        {result['plot']}")
    return 0
