"""Rare Variant Prioritiser — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    vcf_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    max_af: float = 0.01,
    min_cadd: float = 20.0,
    hpo_terms: list[str] | None = None,
    gene_panel: Path | None = None,
    no_plot: bool = False,
) -> dict:
    """Prioritise rare variants by combined evidence for pathogenicity.

    Combines allele frequency, CADD score, gene constraint, HPO phenotype
    matching, and ClinVar evidence into a weighted composite score, then
    assigns each variant to a three-tier priority ranking.

    Args:
        vcf_path: Annotated VCF (VEP/SnpEff) with gnomAD AF and CADD scores. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        max_af: Maximum gnomAD allele frequency to consider a variant rare.
        min_cadd: Minimum CADD Phred score to include.
        hpo_terms: List of HPO term IDs (e.g. ["HP:0001250"]) for phenotype matching.
        gene_panel: Text file of gene symbols, one per line.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'prioritised_variants', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Rare Variant Prioritiser v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_prioritiser.py"
    spec = importlib.util.spec_from_file_location("run_prioritiser", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_prioritiser.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
