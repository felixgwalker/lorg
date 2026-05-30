"""Paralog Cluster Builder — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    proteome: Path | None,
    output_dir: Path,
    demo: bool = False,
    min_identity: float = 30.0,
    inflation: float = 2.0,
    no_plot: bool = False,
) -> dict:
    """Build clusters of paralogous genes from a single-species proteome.

    Performs all-vs-self BLAST, filters by identity and coverage, and clusters
    paralog pairs using the Markov Cluster Algorithm (MCL) to produce gene
    family clusters.

    Args:
        proteome: Single-species protein FASTA. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        min_identity: Minimum sequence identity (%) for paralog pairs.
        inflation: MCL inflation parameter controlling cluster granularity.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'clusters', 'pairs', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Paralog Cluster Builder v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_builder.py"
    spec = importlib.util.spec_from_file_location("run_builder", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_builder.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
