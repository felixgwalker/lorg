"""Guide RNA Specificity Ranker — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    guides_path: Path | None,
    genome_fasta: Path | None,
    offtargets_bed: Path | None,
    output_dir: Path,
    demo: bool = False,
    max_mismatches: int = 3,
    no_plot: bool = False,
) -> dict:
    """Rank guide RNAs by predicted specificity using CFD scoring.

    Enumerates all genome positions with up to `max_mismatches` mismatches
    relative to each guide (seed region + PAM-distal) and computes the
    Cutting Frequency Determination (CFD) score from published mismatch
    penalty matrices.  Guides are ranked by specificity score (1 − sum of
    off-target CFDs) and off-target count.

    Args:
        guides_path: FASTA or TSV of 20 nt guide sequences. None if demo.
        genome_fasta: Reference genome FASTA for off-target search.
        offtargets_bed: Pre-computed off-target BED (alternative to genome scan).
        output_dir: Directory for ranked TSV and optional specificity plot.
        demo: Use synthetic guides and mock off-target data.
        max_mismatches: Maximum mismatches to enumerate (default 3).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'ranked_guides' (list), 'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Guide RNA Specificity Ranker v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_ranker.py"
    spec = importlib.util.spec_from_file_location("run_ranker", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_ranker.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
