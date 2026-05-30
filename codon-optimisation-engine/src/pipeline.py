"""Codon Optimisation Engine — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    protein_fasta: Path | None,
    output_dir: Path,
    demo: bool = False,
    host_organism: str = "Homo sapiens",
    strategy: str = "CAI_maximised",
    avoid_restriction_sites: list[str] | None = None,
    target_gc_min: float = 0.40,
    target_gc_max: float = 0.65,
    no_plot: bool = False,
) -> dict:
    """Optimise protein-coding sequences for expression in a target host organism.

    Replaces each codon with the most suitable synonym for the target host based
    on codon usage tables, maximising CAI while maintaining GC content within
    bounds and avoiding specified restriction enzyme recognition sites.

    Args:
        protein_fasta: Protein FASTA of sequences to optimise. None if demo.
        output_dir: Directory for output FASTA and optional plot.
        demo: Run on synthetic data without real inputs.
        host_organism: Target host organism for codon usage table.
        strategy: Optimisation strategy (most_frequent, harmonised, CAI_maximised).
        avoid_restriction_sites: List of restriction enzyme names to avoid.
        target_gc_min: Minimum GC content of the optimised sequence.
        target_gc_max: Maximum GC content of the optimised sequence.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'optimised_sequences', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Codon Optimisation Engine v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_engine.py"
    spec = importlib.util.spec_from_file_location("run_engine", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_engine.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
