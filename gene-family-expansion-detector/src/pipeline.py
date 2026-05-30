"""Gene Family Expansion Detector — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    family_table: Path | None,
    phylogeny: Path | None,
    output_dir: Path,
    demo: bool = False,
    min_fold_expansion: float = 2.0,
    p_value_threshold: float = 0.05,
    no_plot: bool = False,
) -> dict:
    """Detect gene family expansions on specific lineages.

    Models gene family size evolution under a birth-death process (CAFE-style)
    to identify families that have expanded significantly on particular branches,
    tested against the genome-wide rate of gene gain/loss.

    Args:
        family_table: TSV of gene family sizes per species. None if demo.
        phylogeny: Newick species tree with branch lengths. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        min_fold_expansion: Minimum fold change to report as expanded.
        p_value_threshold: FDR-corrected p-value threshold for significance.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'expansions', 'n_families_assessed', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Gene Family Expansion Detector v%s ===", __version__)
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
