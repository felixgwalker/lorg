"""DNA Fragmentation Profiler — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    bam_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    context_bases: int = 25,
    no_plot: bool = False,
) -> dict:
    """Profile DNA fragmentation and post-mortem damage in an ancient DNA BAM.

    Computes the fragment length distribution and 5'/3' terminal deamination
    profiles (C→T and G→A frequencies) characteristic of ancient DNA, and
    classifies the fragmentation pattern.

    Args:
        bam_path: BAM file of aligned reads. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        context_bases: Number of terminal bases to include in deamination profile.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'fragment_distribution', 'deamination_profile', 'fragmentation_pattern',
        'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== DNA Fragmentation Profiler v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_profiler.py"
    spec = importlib.util.spec_from_file_location("run_profiler", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_profiler.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
