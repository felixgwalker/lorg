"""Genetic Rescue Candidate Selector — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    recipient_vcf: Path | None,
    donor_vcfs: list[Path] | None,
    population_map: Path | None,
    output_dir: Path,
    demo: bool = False,
    recipient_pop: str = "recipient",
    no_plot: bool = False,
) -> dict:
    """Select optimal donor populations for genetic rescue of an inbred population.

    Evaluates each candidate donor population by expected heterozygosity gain in
    the recipient, kinship distance (too close → inbreeding; too far → outbreeding
    depression), ecotype compatibility, and geographic proximity.

    Args:
        recipient_vcf: VCF of the recipient population. None if demo.
        donor_vcfs: List of VCFs for candidate donor populations. None if demo.
        population_map: TSV mapping samples to populations. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        recipient_pop: Population label for the recipient.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'candidates', 'recipient_population_id', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Genetic Rescue Candidate Selector v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_selector.py"
    spec = importlib.util.spec_from_file_location("run_selector", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_selector.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
