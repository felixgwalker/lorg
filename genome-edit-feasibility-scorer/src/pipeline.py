"""Genome Edit Feasibility Scorer — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    spec_json: Path | None,
    output_dir: Path,
    demo: bool = False,
    no_plot: bool = False,
) -> dict:
    """Score the overall feasibility of a genome editing project.

    Reads a project specification JSON and computes a composite feasibility
    score from weighted components: PAM density at the target locus (from
    FASTA if provided), GC content of the protospacer region, chromatin
    accessibility proxy (repeat density / provided ATAC signal), gene
    essentiality score (user-supplied or DepMap lookup), and delivery
    suitability for the specified cell type.

    Args:
        spec_json: JSON with keys: locus_fasta, edit_type, cell_type,
            delivery_method, optional chromatin_score and essentiality_score.
        output_dir: Directory for feasibility report JSON and optional radar plot.
        demo: Use a synthetic project specification.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'feasibility_score' (float), 'components' (dict),
        'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Genome Edit Feasibility Scorer v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_scorer.py"
    spec = importlib.util.spec_from_file_location("run_scorer", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_scorer.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
