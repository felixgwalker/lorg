"""Structural Variant Prioritiser — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    vcf_path: Path | None,
    gene_annotation: Path | None,
    dosage_scores: Path | None,
    output_dir: Path,
    demo: bool = False,
    sv_types: list[str] | None = None,
    min_size: int = 50,
    max_af: float = 0.01,
    haploinsufficiency_threshold: float = 0.9,
    no_plot: bool = False,
) -> dict:
    """Prioritise structural variants by gene overlap, dosage sensitivity, and database evidence.

    Intersects SVs with gene annotations and dosage sensitivity scores (ClinGen HI/TS),
    then scores each SV on rarity, gene content, exon disruption, and overlap with
    DECIPHER or ClinVar pathogenic SVs.

    Args:
        vcf_path: VCF or BED of structural variants. None if demo.
        gene_annotation: Gene annotation BED or GTF. None if demo.
        dosage_scores: ClinGen dosage sensitivity scores TSV. None if demo.
        output_dir: Directory for output TSV and optional plot.
        demo: Run on synthetic data without real inputs.
        sv_types: SV types to include (DEL, DUP, INV, BND, INS, CNV).
        min_size: Minimum SV size in base pairs.
        max_af: Maximum population AF for a rare SV call.
        haploinsufficiency_threshold: ClinGen HI score threshold.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'prioritised_svs', 'output_files', 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== Structural Variant Prioritiser v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    run_path = Path(__file__).parent.parent / "run_prioritiser.py"
    spec = importlib.util.spec_from_file_location("run_prioritiser", run_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_prioritiser.py at {run_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
