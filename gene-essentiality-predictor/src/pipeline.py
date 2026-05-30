"""Gene Essentiality Predictor — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    gene_list: Path | None,
    depmap_scores: Path | None,
    constraint_data: Path | None,
    output_dir: Path,
    demo: bool = False,
    cell_line_context: str | None = None,
    no_plot: bool = False,
) -> dict:
    """Predict gene essentiality by aggregating CRISPR, RNAi, and constraint evidence.

    Combines DepMap CRISPR fitness scores, RNAi dependency scores, gnomAD constraint
    metrics (LOEUF, pLI), and optional cell-line context to produce a composite
    essentiality score and classification.

    Args:
        gene_list: Text file of gene IDs/symbols to assess. None if demo.
        depmap_scores: DepMap CRISPR gene effect score TSV. None optional.
        constraint_data: gnomAD constraint TSV. None optional.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        cell_line_context: DepMap cell line ID for context-specific essentiality.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'predictions', 'n_essential', 'n_context_dependent', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Gene Essentiality Predictor v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_predictor.py"
    spec = importlib.util.spec_from_file_location("run_predictor", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_predictor.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
