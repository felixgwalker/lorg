"""Missense Impact Scorer — pipeline orchestrator."""

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
    conservation_tool: str = "phyloP",
    benign_threshold: float = 0.3,
    pathogenic_threshold: float = 0.7,
    no_plot: bool = False,
) -> dict:
    """Score functional impact of missense variants.

    Combines conservation scores (PhyloP/GERP), BLOSUM62 substitution cost,
    and physicochemical property changes into a composite impact score.
    Classifies each variant on a five-tier ACMG-aligned scale.

    Args:
        vcf_path: VCF of missense variants with HGVS_P annotations. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        conservation_tool: Conservation metric to use (phyloP, GERP, phastCons).
        benign_threshold: Composite score below which variants are called benign.
        pathogenic_threshold: Composite score above which variants are called pathogenic.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'scored_variants', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Missense Impact Scorer v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_scorer.py"
    spec = importlib.util.spec_from_file_location("run_scorer", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_scorer.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
