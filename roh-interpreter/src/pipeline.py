"""ROH Interpreter — pipeline orchestrator."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__
from src.froh_calculator import compute_froh
from src.ne_estimator import estimate_ne, estimate_ne_from_roh
from src.plot import plot_roh, plot_roh_karyotype, plot_froh_distribution
from src.report import write_froh_csv, write_roh_bed, write_roh_report
from src.roh_detector import detect_roh, ROHSegment
from src.vcf_parser import generate_demo_genotypes, GenotypeData, parse_vcf

logger = logging.getLogger(__name__)


def run_pipeline(
    vcf_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    window_size: int = 50,
    homo_threshold: float = 0.95,
    genome_length: int = 2_700_000_000,
    generation_time: float = 25.0,
    no_plot: bool = False,
) -> dict:
    """Execute the full ROH Interpreter pipeline."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=== ROH Interpreter v%s ===", __version__)

    if demo:
        logger.info("Demo mode: generating synthetic genotype data.")
        gdata_list: list[GenotypeData] = generate_demo_genotypes()
    else:
        if vcf_path is None:
            raise ValueError("VCF path required when not in demo mode.")
        logger.info("Parsing VCF: %s", vcf_path)
        gdata_list = parse_vcf(str(vcf_path))
    logger.info(
        "Loaded genotype data for %d individual-chromosome combinations.",
        len(gdata_list),
    )

    all_segments: list[ROHSegment] = []
    for gd in gdata_list:
        segs = detect_roh(gd, window_size=window_size, homo_threshold=homo_threshold)
        all_segments.extend(segs)
    logger.info("Detected %d ROH segments total.", len(all_segments))

    individual_ids = sorted({gd.individual_id for gd in gdata_list})
    froh_results = [
        compute_froh(ind_id, all_segments, genome_length=genome_length)
        for ind_id in individual_ids
    ]
    logger.info("Computed FROH for %d individuals.", len(froh_results))

    # Per-individual NeEstimate objects (new API)
    ne_estimates_per_ind = []
    for ind_id in individual_ids:
        ind_segs = [s for s in all_segments if s.individual_id == ind_id]
        ne_est = estimate_ne_from_roh(
            ind_segs,
            generation_time_years=generation_time,
            genome_length=genome_length,
            individual_id=ind_id,
        )
        ne_estimates_per_ind.append(ne_est)

    # Legacy flat list of _NePoint objects for the trajectory plot
    ne_estimates_flat = estimate_ne(
        froh_results,
        generation_time_years=generation_time,
        genome_length=genome_length,
    )

    # --- Write output files ---
    bed_path = write_roh_bed(all_segments, output_dir)
    froh_path = write_froh_csv(froh_results, output_dir)
    txt_path, report_csv_path = write_roh_report(
        all_segments, froh_results, ne_estimates_per_ind, output_dir
    )
    logger.info("Written: %s", bed_path)
    logger.info("Written: %s", froh_path)
    logger.info("Written: %s", txt_path)
    logger.info("Written: %s", report_csv_path)

    output_files: dict[str, str] = {
        "roh_bed": str(bed_path),
        "froh_csv": str(froh_path),
        "report_txt": str(txt_path),
        "report_csv": str(report_csv_path),
    }

    if no_plot:
        logger.info("Plot generation skipped (--no-plot).")
    else:
        try:
            png, svg = plot_roh(all_segments, froh_results, ne_estimates_flat, output_dir)
            output_files["plot_png"] = str(png)
            output_files["plot_svg"] = str(svg)
            logger.info("Plot written: %s", png)
        except Exception as exc:
            logger.warning("Plot generation failed: %s", exc)

        try:
            # Karyotype plot — use max chrom position as chromosome length
            chrom_sizes: dict[str, int] = {}
            for seg in all_segments:
                chrom_sizes[seg.chrom] = max(
                    chrom_sizes.get(seg.chrom, 0), seg.end_pos
                )
            kp_png, kp_svg = plot_roh_karyotype(all_segments, chrom_sizes, output_dir)
            output_files["karyotype_png"] = str(kp_png)
            output_files["karyotype_svg"] = str(kp_svg)
            logger.info("Karyotype plot written: %s", kp_png)
        except Exception as exc:
            logger.warning("Karyotype plot failed: %s", exc)

        if len(froh_results) > 1:
            try:
                fd_png, fd_svg = plot_froh_distribution(froh_results, output_dir)
                output_files["froh_dist_png"] = str(fd_png)
                output_files["froh_dist_svg"] = str(fd_svg)
                logger.info("FROH distribution plot written: %s", fd_png)
            except Exception as exc:
                logger.warning("FROH distribution plot failed: %s", exc)

    logger.info("Done. Output: %s", output_dir)
    return {
        "genotype_data": gdata_list,
        "segments": all_segments,
        "froh_results": froh_results,
        "ne_estimates": ne_estimates_per_ind,
        "ne_estimates_flat": ne_estimates_flat,
        "output_files": output_files,
        "pipeline_version": __version__,
    }


def main(argv: list[str] | None = None) -> int:
    """Entry point delegated from main.py."""
    import importlib.util
    ra_path = Path(__file__).parent.parent / "run_interpreter.py"
    spec = importlib.util.spec_from_file_location("run_interpreter", ra_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_interpreter.py at {ra_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
