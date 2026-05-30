"""Repair Pathway Bias Estimator — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    fasta_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    cell_type: str = "HEK293",
    cut_position: int | None = None,
    no_plot: bool = False,
) -> dict:
    """Estimate DSB repair pathway probabilities at a CRISPR cut site.

    Enumerates microhomologies (≥2 bp) within ±30 bp of the cut site to
    score MMEJ likelihood.  HDR probability is modulated by a cell-type
    lookup table (dividing vs. non-dividing, S/G2-phase fraction estimates).
    NHEJ probability is set as the remainder.  The output is a pathway
    probability vector plus a ranked table of predicted MMEJ deletion products.

    Args:
        fasta_path: FASTA with ≥50 bp flanking each side of the cut site.
        output_dir: Directory for pathway JSON, MMEJ TSV, and optional pie chart.
        demo: Use a synthetic 120 bp sequence.
        cell_type: Cell type identifier for HDR bias lookup.
        cut_position: 1-based cut position within FASTA (default: centre).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'pathway_probabilities' (dict), 'mmej_products' (list),
        'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Repair Pathway Bias Estimator v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_estimator.py"
    spec = importlib.util.spec_from_file_location("run_estimator", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_estimator.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
