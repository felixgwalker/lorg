"""Ancestral Gene Content Reconstructor — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    presence_absence_table: Path | None,
    phylogeny: Path | None,
    output_dir: Path,
    demo: bool = False,
    method: str = "parsimony",
    no_plot: bool = False,
) -> dict:
    """Reconstruct ancestral gene content at internal nodes of a species phylogeny.

    Applies Dollo parsimony or Bayesian gain/loss model to infer the set of genes
    present in each ancestral genome, reporting posterior probabilities for each
    gene's presence at each internal node.

    Args:
        presence_absence_table: TSV of gene presence/absence per species. None if demo.
        phylogeny: Newick species tree. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        method: Reconstruction method (parsimony, bayesian).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'ancestral_genomes', 'n_genes_assessed', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Ancestral Gene Content Reconstructor v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_reconstructor.py"
    spec = importlib.util.spec_from_file_location("run_reconstructor", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_reconstructor.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
