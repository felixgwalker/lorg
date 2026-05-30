"""Conservation Priority Ranker — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    vcfs: list[Path] | None,
    population_metadata: Path | None,
    output_dir: Path,
    demo: bool = False,
    no_plot: bool = False,
) -> dict:
    """Rank multiple populations by conservation genomic priority.

    Combines inbreeding level, effective population size, adaptive diversity,
    unique allele fraction, and threat status into a composite priority score,
    assigning each population to a conservation tier.

    Args:
        vcfs: List of population VCFs. None if demo.
        population_metadata: TSV of per-population metadata (size, threat status, etc.). None optional.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'ranked_populations', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Conservation Priority Ranker v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_ranker.py"
    spec = importlib.util.spec_from_file_location("run_ranker", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_ranker.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
