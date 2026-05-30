"""pegRNA Optimiser — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    fasta_path: Path | None,
    edit_json: Path | None,
    output_dir: Path,
    demo: bool = False,
    pam: str = "NGG",
    pbs_range: tuple[int, int] = (8, 15),
    rt_range: tuple[int, int] = (10, 16),
    no_plot: bool = False,
) -> dict:
    """Optimise pegRNA design by grid search over PBS and RT length combinations.

    Enumerates all (PBS length, RT length) pairs in the specified ranges,
    computes a feature score for each (PBS GC, RT GC, RT MFE approximation,
    spacer score), and returns Pareto-optimal candidates trading efficiency
    against synthesis complexity.

    Args:
        fasta_path: Target locus FASTA. None if demo.
        edit_json: Edit specification JSON.
        output_dir: Directory for ranked TSV and optional heatmap plot.
        demo: Use synthetic data.
        pam: PAM motif.
        pbs_range: (min, max) PBS length range.
        rt_range: (min, max) RT template length range.
        no_plot: Skip plot generation.

    Returns:
        dict with keys 'candidates' (list), 'pareto_front' (list),
        'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== pegRNA Optimiser v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_optimiser.py"
    spec = importlib.util.spec_from_file_location("run_optimiser", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_optimiser.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
