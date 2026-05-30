"""Chromatin Accessibility Scorer — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    bam_path: Path | None,
    peaks_bed: Path | None,
    output_dir: Path,
    demo: bool = False,
    normalisation: str = "RPM",
    q_value_threshold: float = 0.05,
    no_plot: bool = False,
) -> dict:
    """Score chromatin accessibility at genomic regions from ATAC-seq BAM data.

    Counts ATAC-seq fragment insertions within peak regions, normalises by
    library size, and classifies each region as open, intermediate, or closed
    based on enrichment over background.

    Args:
        bam_path: ATAC-seq BAM file. None if demo.
        peaks_bed: BED of peak regions to score. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        normalisation: Normalisation method (RPM, RPKM, TMM).
        q_value_threshold: FDR threshold for open chromatin classification.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'scores', 'n_open_regions', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Chromatin Accessibility Scorer v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_scorer.py"
    spec = importlib.util.spec_from_file_location("run_scorer", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_scorer.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
