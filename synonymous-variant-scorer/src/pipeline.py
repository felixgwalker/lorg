"""Synonymous Variant Scorer — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__
from src.aggregator import aggregate_scores, ScoredVariant
from src.codon_scorer import score_codon_usage
from src.folding_scorer import score_folding
from src.mrna_scorer import score_mrna_stability
from src.plot import plot_mechanism_scores
from src.report import write_composite_csv, write_mechanism_csv
from src.splicing_scorer import score_splicing
from src.vcf_parser import generate_demo_variants, parse_vcf, SynonymousVariant

logger = logging.getLogger(__name__)


def run_pipeline(
    vcf_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    no_plot: bool = False,
) -> dict:
    """Execute the full Synonymous Variant Scorer pipeline."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=== Synonymous Variant Scorer v%s ===", __version__)

    if demo:
        logger.info("Demo mode: generating synthetic variants.")
        variants: list[SynonymousVariant] = generate_demo_variants(n=10)
    else:
        if vcf_path is None:
            raise ValueError("VCF path required when not in demo mode.")
        logger.info("Parsing VCF: %s", vcf_path)
        variants = parse_vcf(str(vcf_path))
    logger.info("Loaded %d synonymous variants.", len(variants))

    if not variants:
        logger.warning("No synonymous variants found. Outputs will be empty.")

    scored: list[ScoredVariant] = []
    for v in variants:
        spl = score_splicing(v)
        cub = score_codon_usage(v)
        mrna = score_mrna_stability(v)
        fold = score_folding(v)
        sv = aggregate_scores(
            variant_id=v.variant_id,
            chrom=v.chrom,
            pos=v.pos,
            ref_codon=v.ref_codon,
            alt_codon=v.alt_codon,
            gene=v.gene,
            transcript=v.transcript,
            splicing=spl,
            codon_usage=cub,
            mrna_stability=mrna,
            folding=fold,
        )
        scored.append(sv)
    logger.info("Scored %d variants.", len(scored))

    mech_csv = write_mechanism_csv(scored, output_dir)
    comp_csv = write_composite_csv(scored, output_dir)
    logger.info("Written: %s", mech_csv)
    logger.info("Written: %s", comp_csv)

    output_files: dict[str, str] = {
        "mechanism_csv": str(mech_csv),
        "composite_csv": str(comp_csv),
    }

    if no_plot:
        logger.info("Plot generation skipped (--no-plot).")
    else:
        try:
            png, svg = plot_mechanism_scores(scored, output_dir)
            output_files["plot_png"] = str(png)
            output_files["plot_svg"] = str(svg)
            logger.info("Plot written: %s", png)
        except Exception as exc:
            logger.warning("Plot generation failed: %s", exc)

    logger.info("Done. Output: %s", output_dir)
    return {
        "variants": variants,
        "scored": scored,
        "output_files": output_files,
        "pipeline_version": __version__,
    }


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_scorer.py"
    spec = importlib.util.spec_from_file_location("run_scorer", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_scorer.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
