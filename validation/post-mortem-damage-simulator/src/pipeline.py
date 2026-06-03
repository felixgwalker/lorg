"""Post-Mortem Damage Simulator — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    reference_fasta: Path | None,
    output_dir: Path,
    demo: bool = False,
    model: str = "briggs",
    n_reads: int = 100000,
    mean_fragment_length: int = 60,
    deamination_rate_ss: float = 0.68,
    no_plot: bool = False,
) -> dict:
    """Simulate post-mortem damage in ancient DNA reads.

    Generates synthetic FASTQ reads with realistic C→T (5' end) and G→A (3' end)
    deamination patterns, fragment length distributions, and nick frequencies
    according to the Briggs et al. or other damage models.

    Args:
        reference_fasta: Reference genome FASTA to sample reads from. None if demo.
        output_dir: Directory for output FASTQ and optional plot.
        demo: Run on synthetic data without real inputs.
        model: Damage model (briggs, uniform, double_stranded, single_stranded).
        n_reads: Number of reads to simulate.
        mean_fragment_length: Mean fragment length in base pairs.
        deamination_rate_ss: Single-stranded deamination rate for Briggs model.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'simulated_reads', 'damage_parameters', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Post-Mortem Damage Simulator v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_simulator.py"
    spec = importlib.util.spec_from_file_location("run_simulator", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_simulator.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
