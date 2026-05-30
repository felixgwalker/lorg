"""Compound Heterozygosity Detector — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    vcf_path: Path | None,
    ped_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    max_af: float = 0.01,
    sample_id: str | None = None,
    no_plot: bool = False,
) -> dict:
    """Detect compound heterozygous variant pairs in the same gene.

    Identifies pairs of heterozygous variants in the same gene that are on
    opposite haplotypes, using phased genotypes, trio phase-by-transmission,
    or statistical phasing as available.

    Args:
        vcf_path: Phased or unphased VCF (multi-sample or singleton). None if demo.
        ped_path: PLINK PED file describing family structure. None for singletons.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        max_af: Maximum gnomAD AF to consider a variant as a comp-het candidate.
        sample_id: Proband sample ID (required for multi-sample VCF without PED).
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'comp_het_pairs', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Compound Heterozygosity Detector v%s ===", __version__)
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
