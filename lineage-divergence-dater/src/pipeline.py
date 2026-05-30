"""Lineage Divergence Dater — pipeline orchestrator."""

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
    method: str = "molecular_clock",
    mutation_rate: float = 1.25e-8,
    generation_time: float = 30.0,
    no_plot: bool = False,
) -> dict:
    """Date lineage divergence events from sequence alignments and calibrations.

    Estimates split times between lineages using a molecular clock (simple
    pairwise distance / rate), Bayesian dated-tips approach for ancient DNA
    samples with radiocarbon ages, or calibrated node dating.

    Args:
        alignment: FASTA multiple sequence alignment. None if demo.
        phylogeny: Newick phylogenetic tree. None if demo.
        calibrations: JSON of calibration points or dated tip ages. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        method: Dating method (molecular_clock, bayesian_dated_tips, pairwise_distance).
        mutation_rate: Per-base per-generation mutation rate.
        generation_time: Generation time in years.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'divergence_dates', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Lineage Divergence Dater v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_dater.py"
    spec = importlib.util.spec_from_file_location("run_dater", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_dater.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
