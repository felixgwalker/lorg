"""Ortholog Mapper — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    query_proteome: Path | None,
    target_proteomes: list[Path] | None,
    output_dir: Path,
    demo: bool = False,
    method: str = "reciprocal_best_hits",
    min_identity: float = 30.0,
    synteny_support: bool = False,
    no_plot: bool = False,
) -> dict:
    """Map orthologous genes between a query and one or more target species.

    Uses reciprocal best BLAST hits (RBH), OMA, or Inparanoid-style algorithm
    to identify orthologs, classifies relationships (1:1, 1:N, N:1, N:N), and
    optionally validates with synteny context.

    Args:
        query_proteome: FASTA of query species protein sequences. None if demo.
        target_proteomes: List of target species protein FASTAs. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        method: Ortholog detection method (reciprocal_best_hits, OMA, inparanoid).
        min_identity: Minimum sequence identity (%) to accept a hit.
        synteny_support: Require synteny evidence to confirm orthologs.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'ortholog_groups', 'n_query_genes', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Ortholog Mapper v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_mapper.py"
    spec = importlib.util.spec_from_file_location("run_mapper", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_mapper.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
