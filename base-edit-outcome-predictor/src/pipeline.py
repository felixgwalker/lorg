"""Base Edit Outcome Predictor — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    targets_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    editor: str = "CBE3",
    no_plot: bool = False,
) -> dict:
    """Predict base editing outcomes for CBE and ABE editors.

    Identifies editable cytosines (CBE) or adenines (ABE) within the editing
    window (positions 4–8 counting from the spacer 5′ end) and scores each
    by position and trinucleotide context using BE-Hive-derived efficiency
    weights.  Reports primary product probability, bystander edit probability,
    and indel frequency.

    Args:
        targets_path: TSV with columns: id, spacer (20 nt), target_sequence.
        output_dir: Directory for per-base outcome TSV and optional heatmap.
        demo: Use synthetic spacer/target pairs.
        editor: Base editor type ('CBE3', 'BE4max', 'ABE8e', 'ABEmax').
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'outcomes' (list), 'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Base Edit Outcome Predictor v%s ===", __version__)
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
