"""Guide RNA GC Optimiser — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    guides_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    gc_min: float = 0.30,
    gc_max: float = 0.70,
    no_plot: bool = False,
) -> dict:
    """Score and rank guide RNAs by GC content features.

    Scores each guide on: total GC fraction (optimal 30–70 %), seed-region
    GC (positions 1–12 from PAM, optimal 40–60 %), homopolymer run penalty
    (≥4 identical bases), and poly-T run penalty (terminates U6 transcription).
    Returns a composite GC optimality score and per-feature subscores.

    Args:
        guides_path: FASTA or TSV of guide sequences. None if demo.
        output_dir: Directory for scored TSV and optional bar chart.
        demo: Use synthetic guides.
        gc_min/gc_max: Acceptable total GC range for scoring.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'scored_guides' (list), 'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Guide RNA GC Optimiser v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_optimiser.py"
    spec = importlib.util.spec_from_file_location("run_optimiser", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_optimiser.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
