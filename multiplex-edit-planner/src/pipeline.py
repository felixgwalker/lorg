"""Multiplex Edit Planner — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    manifest_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    max_guides_per_batch: int = 4,
    no_plot: bool = False,
) -> dict:
    """Plan a multiplex CRISPR editing strategy for multiple genomic targets.

    Reads an edit manifest (per-target: guide sequence, locus coordinates, edit
    type) and checks for: (1) guide cross-reactivity by CFD score between all
    guide pairs, (2) overlapping cut windows (within 100 bp of each other on
    the same chromosome), (3) potential chromosomal translocations from pairs
    of DSBs on the same chromosome, and (4) delivery capacity constraints.
    Outputs a compatibility matrix and batched delivery plan.

    Args:
        manifest_path: JSON or TSV manifest of edits. None if demo.
        output_dir: Directory for compatibility matrix, batched plan, and plot.
        demo: Use a synthetic 6-target manifest.
        max_guides_per_batch: Maximum guides per delivery batch.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'compatibility_matrix' (list[list]), 'batches' (list),
        'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Multiplex Edit Planner v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_planner.py"
    spec = importlib.util.spec_from_file_location("run_planner", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_planner.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
