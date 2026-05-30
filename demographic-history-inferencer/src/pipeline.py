"""Demographic History Inferencer — pipeline orchestrator."""

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
    models: list[str] | None = None,
    generation_time: float = 30.0,
    mutation_rate: float = 1.25e-8,
    no_plot: bool = False,
) -> dict:
    """Infer demographic history from the site frequency spectrum.

    Fits parametric demographic models (constant, exponential growth, two-epoch,
    three-epoch) to the observed SFS using diffusion approximations or moment-based
    inference, selecting the best model by AIC.

    Args:
        vcf_path: VCF of SNPs for one population. None if demo.
        output_dir: Directory for output and optional plot.
        demo: Run on synthetic data without real inputs.
        models: List of demographic models to fit.
        generation_time: Assumed generation time in years.
        mutation_rate: Per-base per-generation mutation rate.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'best_model', 'all_fits', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Demographic History Inferencer v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_inferencer.py"
    spec = importlib.util.spec_from_file_location("run_inferencer", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_inferencer.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
