"""Gene Regulatory Network Builder — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    expression_matrix: Path | None,
    tf_list: Path | None,
    output_dir: Path,
    demo: bool = False,
    method: str = "GENIE3",
    min_edge_weight: float = 0.01,
    no_plot: bool = False,
) -> dict:
    """Build a gene regulatory network from expression data.

    Infers regulatory edges between transcription factors and target genes
    using GENIE3 (random forest importance), ARACNE (mutual information), or
    Pearson correlation, and reports a ranked edge list with hub genes.

    Args:
        expression_matrix: Gene expression matrix (genes × samples). None if demo.
        tf_list: Text file of transcription factor gene IDs. None optional.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        method: Network inference method (GENIE3, correlation, ARACNE).
        min_edge_weight: Minimum edge weight to include in output.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'edges', 'nodes', 'n_tf_regulators', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Gene Regulatory Network Builder v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_builder.py"
    spec = importlib.util.spec_from_file_location("run_builder", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_builder.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
