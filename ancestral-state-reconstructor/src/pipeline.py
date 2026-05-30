"""Ancestral State Reconstructor — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    alignment: Path | None,
    phylogeny: Path | None,
    output_dir: Path,
    demo: bool = False,
    method: str = "maximum_likelihood",
    substitution_model: str = "GTR+G",
    no_plot: bool = False,
) -> dict:
    """Reconstruct ancestral nucleotide or amino acid sequences at phylogenetic nodes.

    Applies marginal maximum likelihood, parsimony, or Bayesian ancestral state
    reconstruction to a multiple sequence alignment and phylogeny, reporting
    the reconstructed sequence at each internal node with posterior probabilities.

    Args:
        alignment: FASTA multiple sequence alignment. None if demo.
        phylogeny: Newick phylogenetic tree. None if demo.
        output_dir: Directory for output FASTA and TSV.
        demo: Run on synthetic data without real inputs.
        method: Reconstruction method (parsimony, maximum_likelihood, bayesian).
        substitution_model: Substitution model for ML (e.g. GTR+G, JC).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'ancestral_sequences', 'site_states', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Ancestral State Reconstructor v%s ===", __version__)
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
