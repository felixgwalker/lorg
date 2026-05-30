"""Coexpression Module Finder — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    expression_matrix: Path | None,
    trait_data: Path | None,
    output_dir: Path,
    demo: bool = False,
    method: str = "WGCNA",
    min_module_size: int = 30,
    no_plot: bool = False,
) -> dict:
    """Find coexpression modules from an expression matrix.

    Constructs a weighted gene coexpression network (WGCNA-style) or uses
    clique-based / k-means clustering to identify groups of co-expressed genes,
    optionally correlating module eigengenes with sample traits.

    Args:
        expression_matrix: Normalised expression matrix (genes × samples). None if demo.
        trait_data: TSV of sample traits for module-trait correlation. None optional.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        method: Module detection method (WGCNA, clique_based, kmeans).
        min_module_size: Minimum number of genes per module.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'modules', 'module_trait_associations', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Coexpression Module Finder v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_finder.py"
    spec = importlib.util.spec_from_file_location("run_finder", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_finder.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
