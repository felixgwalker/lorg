"""Edit Outcome Simulator — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    target_fasta: Path | None,
    output_dir: Path,
    demo: bool = False,
    guide: str | None = None,
    cas_variant: str = "SpCas9",
    n_simulations: int = 10000,
    no_plot: bool = False,
) -> dict:
    """Simulate CRISPR indel outcome distribution at a target site.

    Uses an inDelphi-style approach: enumerates 1-bp templated insertions
    (copying the nucleotide at position +1 relative to the cut), and deletions
    of 1–30 bp weighted by microhomology scores (length² × GC factor).
    Predicts frameshift frequency for coding targets when reading frame is
    provided.  Returns a probability distribution over indel outcomes.

    Args:
        target_fasta: FASTA of ≥60 bp flanking the cut site. None if demo.
        output_dir: Directory for indel distribution TSV and optional stacked bar.
        demo: Use a synthetic target sequence.
        guide: 20 nt guide sequence (used to locate cut site in FASTA).
        cas_variant: Cas nuclease variant for cut position offset.
        n_simulations: Number of Monte Carlo draws for indel sampling.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'indel_distribution' (list), 'frameshift_rate' (float),
        'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Edit Outcome Simulator v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_simulator.py"
    spec = importlib.util.spec_from_file_location("run_simulator", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_simulator.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
