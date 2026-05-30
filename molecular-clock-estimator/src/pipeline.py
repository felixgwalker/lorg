"""Molecular Clock Estimator — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    alignment: Path | None,
    phylogeny: Path | None,
    calibrations: Path | None,
    output_dir: Path,
    demo: bool = False,
    clock_model: str = "relaxed_lognormal",
    substitution_model: str = "GTR+G",
    no_plot: bool = False,
) -> dict:
    """Estimate molecular clock rate and test strict vs. relaxed clock models.

    Fits strict and relaxed lognormal/exponential clock models to a calibrated
    alignment + phylogeny using Bayesian MCMC, comparing models by Bayes factors.

    Args:
        alignment: FASTA multiple sequence alignment. None if demo.
        phylogeny: Newick starting tree. None if demo.
        calibrations: JSON of calibration constraints (node, min/max age, prior). None if demo.
        output_dir: Directory for output and optional plot.
        demo: Run on synthetic data without real inputs.
        clock_model: Clock model (strict, relaxed_lognormal, relaxed_exponential).
        substitution_model: Nucleotide substitution model.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'estimates', 'best_model', 'calibration_points', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Molecular Clock Estimator v%s ===", __version__)
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
