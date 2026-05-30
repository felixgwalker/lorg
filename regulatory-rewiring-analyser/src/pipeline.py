"""Regulatory Rewiring Analyser — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    elements_a: Path | None,
    elements_b: Path | None,
    ortholog_table: Path | None,
    output_dir: Path,
    demo: bool = False,
    min_conservation: float = 0.7,
    no_plot: bool = False,
) -> dict:
    """Analyse regulatory rewiring between two species by comparing enhancers and promoters.

    Identifies orthologous genes whose regulatory elements (enhancers, promoters)
    have been gained, lost, conserved, or relocated between two species, using
    sequence conservation and synteny context.

    Args:
        elements_a: BED of regulatory elements in species A. None if demo.
        elements_b: BED of regulatory elements in species B. None if demo.
        ortholog_table: TSV of ortholog gene pairs. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        min_conservation: Minimum sequence identity to call an element conserved.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'rewirings', 'n_elements_assessed', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Regulatory Rewiring Analyser v%s ===", __version__)
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
