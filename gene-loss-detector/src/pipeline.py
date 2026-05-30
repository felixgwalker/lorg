"""Gene Loss Detector — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    ortholog_table: Path | None,
    species_proteomes: dict[str, Path] | None,
    phylogeny: Path | None,
    output_dir: Path,
    demo: bool = False,
    min_species_with_gene: int = 3,
    no_plot: bool = False,
) -> dict:
    """Detect gene losses on specific branches of a phylogenetic tree.

    Identifies genes present in the majority of species but absent or
    pseudogenised in one or more lineages, inferring the branch of loss
    using parsimony on the supplied species tree.

    Args:
        ortholog_table: TSV of ortholog presence/absence per species. None if demo.
        species_proteomes: Dict of species name → protein FASTA path. None if demo.
        phylogeny: Newick species tree. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        min_species_with_gene: Minimum species that must have the gene to assess loss.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'gene_losses', 'n_genes_assessed', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Gene Loss Detector v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_detector.py"
    spec = importlib.util.spec_from_file_location("run_detector", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_detector.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
