"""Microhomology Repair Predictor — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    fasta_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    cut_position: int | None = None,
    min_mh_length: int = 2,
    search_window: int = 30,
    no_plot: bool = False,
) -> dict:
    """Predict MMEJ repair outcomes from microhomology sequences flanking a DSB.

    Enumerates all microhomologies of length ≥ `min_mh_length` within
    `search_window` bp of each side of the cut site.  Scores each match by
    MH score = length² × (1 + GC_fraction), ranks predicted deletion products,
    and estimates relative frequency using a normalised softmax over scores.

    Args:
        fasta_path: FASTA with cut-site sequence (≥60 bp each side). None if demo.
        output_dir: Directory for ranked MH TSV and optional lollipop plot.
        demo: Use a synthetic flanking sequence.
        cut_position: 1-based cut position in FASTA (default: centre).
        min_mh_length: Minimum microhomology length to consider (default 2).
        search_window: Search window (bp) on each side of cut (default 30).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'mh_products' (list), 'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Microhomology Repair Predictor v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_predictor.py"
    spec = importlib.util.spec_from_file_location("run_predictor", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_predictor.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
