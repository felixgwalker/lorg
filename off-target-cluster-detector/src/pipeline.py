"""Off-Target Cluster Detector — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    offtargets_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    window_size: int = 100_000,
    min_cluster_size: int = 3,
    no_plot: bool = False,
) -> dict:
    """Detect genomic clusters of CRISPR off-target sites.

    Reads a BED or TSV of off-target sites (chromosome, position, CFD score)
    and applies a sliding window (default 100 kb) to count site density.
    Regions above the density threshold are further resolved by DBSCAN-style
    single-linkage clustering with an inter-site distance threshold.  Clusters
    are annotated by their overlap with regulatory, genic, and repeat regions
    if an annotation BED is provided.

    Args:
        offtargets_path: BED/TSV of off-target sites with score column. None if demo.
        output_dir: Directory for cluster BED, summary TSV, and optional Manhattan plot.
        demo: Use synthetic off-target site list.
        window_size: Sliding window size in bp for density scan.
        min_cluster_size: Minimum sites per cluster to report.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'clusters' (list), 'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Off-Target Cluster Detector v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_detector.py"
    spec = importlib.util.spec_from_file_location("run_detector", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_detector.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
