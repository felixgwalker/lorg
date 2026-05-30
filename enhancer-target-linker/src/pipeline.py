"""Enhancer Target Linker — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    enhancers: Path | None,
    gene_tss: Path | None,
    activity_matrix: Path | None,
    hic_matrix: Path | None,
    output_dir: Path,
    demo: bool = False,
    method: str = "activity_by_contact",
    max_distance_bp: int = 500000,
    no_plot: bool = False,
) -> dict:
    """Link enhancers to their likely target genes using activity-by-contact or correlation.

    Scores enhancer-gene pairs within a distance window using the ABC model
    (activity × contact), expression correlation, or simple proximity, and
    reports the top-scoring target gene per enhancer.

    Args:
        enhancers: BED of enhancer elements with activity scores. None if demo.
        gene_tss: BED of gene TSS positions. None if demo.
        activity_matrix: Matrix of enhancer activity across samples (for correlation). None optional.
        hic_matrix: Hi-C contact matrix (for ABC model). None optional.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        method: Linking method (activity_by_contact, correlation, distance, hi_c).
        max_distance_bp: Maximum enhancer-TSS distance to consider.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'links', 'n_enhancers', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Enhancer Target Linker v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_linker.py"
    spec = importlib.util.spec_from_file_location("run_linker", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_linker.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
