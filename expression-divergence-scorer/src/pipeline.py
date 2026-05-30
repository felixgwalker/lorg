"""Expression Divergence Scorer — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    expression_a: Path | None,
    expression_b: Path | None,
    ortholog_table: Path | None,
    output_dir: Path,
    demo: bool = False,
    metric: str = "log_fold_change",
    no_plot: bool = False,
) -> dict:
    """Score expression divergence between orthologous genes in two species.

    Compares expression profiles of orthologous gene pairs across matched
    tissues or conditions using log fold change, Jensen-Shannon Index, tau
    (tissue specificity index), or Euclidean distance.

    Args:
        expression_a: TPM/count matrix for species A (genes × tissues). None if demo.
        expression_b: TPM/count matrix for species B (genes × tissues). None if demo.
        ortholog_table: TSV of ortholog gene pairs. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        metric: Divergence metric (log_fold_change, tau, JSI, euclidean).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'divergence_scores', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Expression Divergence Scorer v%s ===", __version__)
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
