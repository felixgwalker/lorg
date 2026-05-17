"""
CNV Significance Assessor — top-level pipeline orchestrator.

Calls all analysis modules in dependency order and returns a structured
result dict.  Also exposes main() so that main.py can delegate via:

    from pipeline import main; sys.exit(main())
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__
from src.annotator import annotate_cnvs, load_population_cnvs
from src.classifier import classify_cnvs
from src.cnv_parser import parse_cnvs
from src.config import (
    DEFAULT_MIN_CNV_SIZE,
    DEFAULT_OVERLAP_FRACTION,
    DEFAULT_POP_FREQ_CUTOFF,
)
from src.dosage_loader import load_dosage_scores
from src.gff_parser import parse_gff3
from src.output_writer import (
    write_annotated_csv,
    write_gene_impact_report,
    write_significance_json,
    write_significance_txt,
)
from src.plot_generator import generate_ideogram_plot

logger = logging.getLogger(__name__)


def run_pipeline(
    cnv_path: Path,
    gff_path: Path,
    output_dir: Path,
    pop_db_path: Path | None = None,
    dosage_path: Path | None = None,
    min_cnv_size: int = DEFAULT_MIN_CNV_SIZE,
    overlap_fraction: float = DEFAULT_OVERLAP_FRACTION,
    pop_freq_cutoff: float = DEFAULT_POP_FREQ_CUTOFF,
    sample_name: str = "",
    no_plot: bool = False,
) -> dict:
    """
    Execute the complete CNV Significance Assessor pipeline.

    Steps:
        1. Parse CNVs      — load BED or VCF input, apply size filter.
        2. Parse GFF3      — load gene bodies and regulatory elements.
        3. Load dosage     — load pLI/pHaplo/pTriplo scores (optional).
        4. Load pop DB     — load population CNV frequency database (optional).
        5. Annotate        — intersect CNVs with genes; look up pop frequency.
        6. Classify        — rule-based significance scoring.
        7. Write outputs   — annotated CSV, summary TXT/JSON, gene impact CSV.
        8. Generate plot   — chromosome ideogram PNG + SVG (unless --no-plot).

    Args:
        cnv_path:         Path to CNV calls (BED or VCF).
        gff_path:         Path to gene annotation (GFF3).
        output_dir:       Directory for all output files (created if absent).
        pop_db_path:      Optional population CNV frequency database (BED/VCF).
        dosage_path:      Optional haploinsufficiency/triplosensitivity CSV.
        min_cnv_size:     Minimum CNV size in bp (smaller CNVs are excluded).
        overlap_fraction: Minimum gene–CNV reciprocal overlap fraction.
        pop_freq_cutoff:  CNVs above this population frequency → LIKELY_BENIGN.
        sample_name:      Label for plot titles and summary headers.
        no_plot:          If True, skip matplotlib rendering.

    Returns:
        Structured result dict with keys:
            "cnvs"             → list[CNVRecord]
            "genes"            → list[GeneRecord]
            "annotated"        → list[AnnotatedCNV]
            "scored"           → list[ScoredCNV]
            "summary"          → ClassificationSummary
            "output_files"     → dict[str, str] mapping type to path
            "pipeline_version" → str
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not sample_name:
        sample_name = Path(cnv_path).stem

    logger.info("=== CNV Significance Assessor v%s ===", __version__)
    logger.info("CNV input    : %s", cnv_path)
    logger.info("GFF3 input   : %s", gff_path)
    logger.info("Output dir   : %s", output_dir)
    logger.info(
        "Parameters   : min_size=%d bp, overlap_frac=%.2f, pop_cutoff=%.4f",
        min_cnv_size, overlap_fraction, pop_freq_cutoff,
    )

    # ── Step 1: Parse CNVs ───────────────────────────────────────────────
    logger.info("Step 1/7 — Parsing CNV input...")
    cnvs = parse_cnvs(Path(cnv_path), min_size=min_cnv_size)
    logger.info("  %d CNVs loaded.", len(cnvs))

    if not cnvs:
        logger.warning("No CNVs passed the size filter.  Outputs will be empty.")

    # ── Step 2: Parse GFF3 ───────────────────────────────────────────────
    logger.info("Step 2/7 — Parsing GFF3 annotation...")
    genes = parse_gff3(Path(gff_path))
    logger.info("  %d annotation records loaded.", len(genes))

    # ── Step 3: Load dosage scores (optional) ────────────────────────────
    dosage_scores: dict = {}
    if dosage_path is not None:
        logger.info("Step 3/7 — Loading dosage sensitivity scores...")
        dosage_scores = load_dosage_scores(Path(dosage_path))
        logger.info("  %d genes with dosage scores.", len(dosage_scores))
    else:
        logger.info("Step 3/7 — Dosage scores not provided; skipping.")

    # ── Step 4: Load population database (optional) ───────────────────────
    pop_cnvs = None
    if pop_db_path is not None:
        logger.info("Step 4/7 — Loading population CNV database...")
        pop_cnvs = load_population_cnvs(Path(pop_db_path), min_size=0)
    else:
        logger.info("Step 4/7 — Population database not provided; skipping.")

    # ── Step 5: Annotate ─────────────────────────────────────────────────
    logger.info("Step 5/7 — Annotating CNVs...")
    annotated = annotate_cnvs(
        cnvs=cnvs,
        genes=genes,
        dosage_scores=dosage_scores,
        pop_cnvs=pop_cnvs,
        overlap_fraction=overlap_fraction,
    )

    # ── Step 6: Classify ─────────────────────────────────────────────────
    logger.info("Step 6/7 — Classifying CNVs...")
    scored, summary = classify_cnvs(annotated, pop_freq_cutoff=pop_freq_cutoff)
    logger.info(
        "  LIKELY_BENIGN: %d  VUS: %d  LIKELY_PATHOGENIC: %d",
        summary.n_likely_benign, summary.n_vus, summary.n_likely_pathogenic,
    )

    # ── Step 7: Write outputs ────────────────────────────────────────────
    logger.info("Step 7/7 — Writing output files...")

    args_dict = {
        "min_cnv_size":     min_cnv_size,
        "overlap_fraction": overlap_fraction,
        "pop_freq_cutoff":  pop_freq_cutoff,
        "pop_db":           str(pop_db_path) if pop_db_path else None,
        "dosage_csv":       str(dosage_path) if dosage_path else None,
    }

    csv_path  = write_annotated_csv(scored, output_dir)
    gene_path = write_gene_impact_report(scored, output_dir)
    txt_path  = write_significance_txt(
        summary, scored, output_dir, str(cnv_path), args_dict
    )
    json_path = write_significance_json(
        summary, scored, output_dir, str(cnv_path), args_dict
    )

    output_files: dict[str, str] = {
        "annotated_csv":          str(csv_path),
        "gene_impact_csv":        str(gene_path),
        "significance_txt":       str(txt_path),
        "significance_json":      str(json_path),
    }

    # ── Step 8: Plot (optional) ──────────────────────────────────────────
    if no_plot:
        logger.info("Plot generation skipped (--no-plot).")
    else:
        try:
            png_path, svg_path = generate_ideogram_plot(
                scored=scored,
                output_dir=output_dir,
                sample_name=sample_name,
            )
            output_files["ideogram_png"] = str(png_path)
            output_files["ideogram_svg"] = str(svg_path)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Plot generation failed: %s. Continuing without plot.", exc)

    logger.info("Done.  Output directory: %s", output_dir)

    return {
        "cnvs":             cnvs,
        "genes":            genes,
        "annotated":        annotated,
        "scored":           scored,
        "summary":          summary,
        "output_files":     output_files,
        "pipeline_version": __version__,
    }


def main(argv: list[str] | None = None) -> int:
    """
    Entry point delegated to from main.py.

    Imports run_assessor.main() so that all CLI logic lives in one place.

    Args:
        argv: Argument list (defaults to sys.argv[1:] when None).

    Returns:
        Integer exit code.
    """
    import importlib.util
    import os

    ra_path = Path(__file__).parent.parent / "run_assessor.py"
    spec = importlib.util.spec_from_file_location("run_assessor", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_assessor.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
