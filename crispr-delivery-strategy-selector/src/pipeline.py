"""CRISPR Delivery Strategy Selector — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__

logger = logging.getLogger(__name__)


def run_pipeline(
    spec_json: Path | None,
    output_dir: Path,
    demo: bool = False,
    no_plot: bool = False,
) -> dict:
    """Rank CRISPR delivery strategies for a given cell type and editing payload.

    Reads a delivery specification JSON (cell type, organism, editing payload:
    RNP/plasmid/mRNA, target tissue, payload size in bp) and scores five
    delivery modalities — plasmid transfection, ribonucleoprotein (RNP)
    electroporation, lipid nanoparticle (LNP), adeno-associated virus (AAV),
    and lentiviral transduction — using a rule-based compatibility table.
    Scores are normalised to 0–1 and include notes on efficiency, risk,
    and size constraints.

    Args:
        spec_json: JSON with keys: cell_type, organism, payload_type,
            payload_size_bp, target_tissue, allow_integration (bool).
        output_dir: Directory for strategy ranking TSV, notes JSON, and optional plot.
        demo: Use a synthetic HEK293 editing scenario.
        no_plot: Skip figure generation.

    Returns:
        dict with keys 'strategies' (list), 'output_files' (dict), 'pipeline_version'.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== CRISPR Delivery Strategy Selector v%s ===", __version__)
    raise NotImplementedError("run_pipeline is not yet implemented.")


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_selector.py"
    spec = importlib.util.spec_from_file_location("run_selector", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_selector.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
