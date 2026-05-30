"""Promoter Strength Estimator — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    gene_annotation: Path | None,
    fasta: Path | None,
    chip_h3k4me3: Path | None,
    output_dir: Path,
    demo: bool = False,
    promoter_window: int = 2000,
    no_plot: bool = False,
) -> dict:
    """Estimate promoter strength from sequence features and optional ChIP-seq data.

    Scores each gene's promoter region using TATA box strength, initiator element
    (Inr) score, CpG island overlap, GC content, TFBS density, and optional
    H3K4me3/Pol II ChIP-seq signal.

    Args:
        gene_annotation: GTF or BED of gene TSS positions. None if demo.
        fasta: Reference genome FASTA. None if demo.
        chip_h3k4me3: H3K4me3 ChIP-seq bigWig (optional). None optional.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        promoter_window: Upstream window (bp) to define the promoter region.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'promoter_strengths', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Promoter Strength Estimator v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_estimator.py"
    spec = importlib.util.spec_from_file_location("run_estimator", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_estimator.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
