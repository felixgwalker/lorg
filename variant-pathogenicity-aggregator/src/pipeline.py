"""Variant Pathogenicity Aggregator — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    vcf_path: Path | None,
    tool_scores: Path | None,
    output_dir: Path,
    demo: bool = False,
    clinvar_stars: int = 1,
    classification_threshold: int = 4,
    no_plot: bool = False,
) -> dict:
    """Aggregate ACMG/AMP evidence from multiple sources into a composite classification.

    Collects ClinVar pathogenicity evidence, in silico tool scores (CADD, REVEL,
    SpliceAI, etc.) and maps them to applicable ACMG/AMP criteria, then applies
    the weighted point-based classification scheme to assign a five-tier class.

    Args:
        vcf_path: Annotated VCF with ClinVar and tool score INFO fields. None if demo.
        tool_scores: Optional TSV of additional in silico tool scores per variant.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        clinvar_stars: Minimum ClinVar review star count to trust ClinVar evidence.
        classification_threshold: Minimum point total to reach likely-pathogenic.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'results', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Variant Pathogenicity Aggregator v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_aggregator.py"
    spec = importlib.util.spec_from_file_location("run_aggregator", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_aggregator.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
