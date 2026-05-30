"""Enhancer Conservation Analyser — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    enhancers: Path | None,
    conservation_bigwig: Path | None,
    alignment_maf: Path | None,
    output_dir: Path,
    demo: bool = False,
    min_phastcons: float = 0.4,
    no_plot: bool = False,
) -> dict:
    """Analyse evolutionary conservation of enhancer elements.

    Retrieves phastCons/PhyloP scores for each enhancer, assesses sequence
    conservation via multiple alignment coverage, and classifies each element
    as highly conserved, moderately conserved, or lineage-specific.

    Args:
        enhancers: BED of enhancer elements. None if demo.
        conservation_bigwig: PhastCons or PhyloP bigWig file. None optional.
        alignment_maf: MAF multiple alignment file. None optional.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        min_phastcons: Minimum mean phastCons to call moderate conservation.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'conserved_enhancers', 'n_elements_assessed', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Enhancer Conservation Analyser v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_analyser.py"
    spec = importlib.util.spec_from_file_location("run_analyser", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_analyser.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
