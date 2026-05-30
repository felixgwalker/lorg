"""Population Viability Genomics Estimator — pipeline orchestrator."""

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
    census_size: int | None = None,
    generation_time: float = 5.0,
    time_horizons: list[int] | None = None,
    no_plot: bool = False,
) -> dict:
    """Estimate population viability using genomic metrics.

    Integrates inbreeding level, effective population size, genetic load,
    and adaptive diversity to project extinction probability over multiple
    time horizons and estimate minimum viable population size.

    Args:
        vcf_path: VCF of the population. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        census_size: Current census population size (optional).
        generation_time: Generation time in years.
        time_horizons: Years to project viability (default [50, 100, 200]).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'metrics', 'projections', 'minimum_viable_population_genomic',
        'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Population Viability Genomics Estimator v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_estimator.py"
    spec = importlib.util.spec_from_file_location("run_estimator", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_estimator.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
