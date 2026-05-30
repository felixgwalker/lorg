"""Conserved Synteny Detector — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    ortholog_table: Path | None,
    gene_positions_a: Path | None,
    gene_positions_b: Path | None,
    output_dir: Path,
    demo: bool = False,
    min_anchors: int = 5,
    min_block_length_kb: int = 100,
    no_plot: bool = False,
) -> dict:
    """Detect conserved synteny blocks between two genomes.

    Uses ortholog gene pairs as anchors and a sliding-window collinearity
    algorithm to identify regions where gene order and orientation are
    conserved between two species.

    Args:
        ortholog_table: TSV of ortholog gene pairs between species A and B. None if demo.
        gene_positions_a: BED of gene positions in species A. None if demo.
        gene_positions_b: BED of gene positions in species B. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        min_anchors: Minimum anchor genes to form a synteny block.
        min_block_length_kb: Minimum block length in kilobases.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'synteny_blocks', 'anchors', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Conserved Synteny Detector v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_detector.py"
    spec = importlib.util.spec_from_file_location("run_detector", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_detector.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
