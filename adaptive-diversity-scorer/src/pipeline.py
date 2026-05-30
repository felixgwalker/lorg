"""Adaptive Diversity Scorer — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    vcf_path: Path | None,
    environment_data: Path | None,
    adaptive_loci_bed: Path | None,
    output_dir: Path,
    demo: bool = False,
    no_plot: bool = False,
) -> dict:
    """Score adaptive genetic diversity in a conservation population.

    Identifies putatively adaptive loci (Fst outliers, environmental associations),
    computes adaptive heterozygosity relative to neutral heterozygosity, and
    classifies the population's adaptive diversity level.

    Args:
        vcf_path: VCF of the population. None if demo.
        environment_data: TSV of environmental variables per sample. None optional.
        adaptive_loci_bed: Pre-identified adaptive loci BED. None optional.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'score', 'adaptive_loci', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Adaptive Diversity Scorer v%s ===", __version__)
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
