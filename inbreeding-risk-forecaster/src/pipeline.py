"""Inbreeding Risk Forecaster — pipeline orchestrator."""

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
    min_roh_kb: int = 500,
    generation_time: float = 5.0,
    n_generations_forecast: int = 50,
    no_plot: bool = False,
) -> dict:
    """Forecast inbreeding risk for a population from genomic data.

    Estimates current inbreeding (FROH from runs of homozygosity, Fis from
    heterozygosity), infers effective population size, and projects inbreeding
    accumulation over future generations.

    Args:
        vcf_path: VCF of the population to assess. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        min_roh_kb: Minimum ROH length (kb) to include in FROH calculation.
        generation_time: Generation time in years (for forecasting).
        n_generations_forecast: Number of generations to project inbreeding.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'stats', 'forecast', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Inbreeding Risk Forecaster v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_forecaster.py"
    spec = importlib.util.spec_from_file_location("run_forecaster", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_forecaster.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
