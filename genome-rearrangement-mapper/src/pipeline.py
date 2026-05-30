"""Genome Rearrangement Mapper — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    synteny_blocks: Path | None,
    output_dir: Path,
    demo: bool = False,
    rearrangement_types: list[str] | None = None,
    no_plot: bool = False,
) -> dict:
    """Map chromosomal rearrangements between two genomes from synteny block data.

    Identifies inversions, translocations, fusions, and fissions by analysing
    the orientation and chromosomal assignment of synteny blocks, reporting
    breakpoint coordinates and the type of each rearrangement.

    Args:
        synteny_blocks: TSV of synteny blocks between species A and B. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        rearrangement_types: Types to detect (inversion, translocation, fusion, fission).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'rearrangements', 'n_breakpoints', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Genome Rearrangement Mapper v%s ===", __version__)
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
