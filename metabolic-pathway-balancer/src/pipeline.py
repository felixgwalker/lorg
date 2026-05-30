"""Metabolic Pathway Balancer — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    pathway_json: Path | None,
    output_dir: Path,
    demo: bool = False,
    objective: str = "maximise_yield",
    no_plot: bool = False,
) -> dict:
    """Balance a heterologous metabolic pathway for optimal flux and yield.

    Runs flux balance analysis (FBA) on the pathway, identifies flux-limiting
    steps (bottlenecks), cofactor imbalances, and toxic intermediate
    accumulation, recommending enzyme expression adjustments.

    Args:
        pathway_json: JSON of pathway reactions with stoichiometry and kinetics. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        objective: Optimisation objective (maximise_yield, maximise_productivity, minimise_byproducts).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'balanced_pathway', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Metabolic Pathway Balancer v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_balancer.py"
    spec = importlib.util.spec_from_file_location("run_balancer", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_balancer.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
