"""
Ancient DNA Damage Classifier — top-level pipeline orchestrator.

Calls all analysis modules in dependency order and returns a structured
result dict.  Also exposes main() so that main.py can delegate via:

    from pipeline import main; sys.exit(main())
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from src import __version__
from src.classifier import classify_reads
from src.config import (
    AUTH_THRESHOLD,
    BACKGROUND_ERROR_RATE,
    CONT_THRESHOLD,
    DEFAULT_MIN_LENGTH,
    DEFAULT_MIN_MAPQ,
    DEFAULT_N_TERMINAL,
    DEFAULT_PRIOR_ANCIENT,
)
from src.damage_profiler import generate_demo_profile, profile_damage
from src.decay_model import fit_model
from src.fragment_length import compute_fragment_lengths
from src.output_writer import (
    write_damage_csv,
    write_read_tsv,
    write_summary_json,
    write_summary_txt,
)
from src.plot_generator import generate_damage_plot

logger = logging.getLogger(__name__)


def run_pipeline(
    bam_path: Path,
    output_dir: Path,
    min_mapq: int = DEFAULT_MIN_MAPQ,
    min_length: int = DEFAULT_MIN_LENGTH,
    n_terminal: int = DEFAULT_N_TERMINAL,
    prior_ancient: float = DEFAULT_PRIOR_ANCIENT,
    auth_threshold: float = AUTH_THRESHOLD,
    cont_threshold: float = CONT_THRESHOLD,
    sample_name: str = "",
    no_plot: bool = False,
) -> dict:
    """
    Execute the complete Ancient DNA Damage Classifier pipeline.

    Steps:
        1. Profile damage   — traverse BAM, compute per-position C→T/G→A arrays
                              and collect per-read terminal mismatch features.
        2. Fit decay model  — fit geometric decay to both terminus profiles,
                              derive library_deamination_rate.
        3. Classify reads   — Bayesian per-read authenticity classification.
        4. Fragment lengths — extract template/query length distribution.
        5. Write outputs    — CSV, TSV, JSON, TXT files.
        6. Generate plot    — two-panel damage profile PNG and SVG (unless no_plot).

    Args:
        bam_path:        Path to coordinate-sorted, indexed BAM file.
        output_dir:      Directory for all output files (created if absent).
        min_mapq:        Minimum mapping quality filter.
        min_length:      Minimum read length filter (bp).
        n_terminal:      Number of terminal positions to profile.
        prior_ancient:   Bayesian prior P(ancient) for the library.
        auth_threshold:  Posterior threshold for "authentic" label.
        cont_threshold:  Posterior threshold for "contaminated" label.
        sample_name:     Optional label for plot titles and reports.
        no_plot:         If True, skip matplotlib rendering.

    Returns:
        Structured result dict with keys:
            "profile"           -> DamageProfile
            "model"             -> ModelResult
            "classifications"   -> list[ReadClassification]
            "summary"           -> ClassificationSummary
            "frag_stats"        -> FragmentLengthStats
            "output_files"      -> dict[str, str] mapping type to path
            "pipeline_version"  -> str
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    bam_path = Path(bam_path)
    if not sample_name:
        sample_name = bam_path.stem

    logger.info("=== Ancient DNA Damage Classifier v%s ===", __version__)
    logger.info("Input BAM    : %s", bam_path)
    logger.info("Output dir   : %s", output_dir)
    logger.info(
        "Filters      : min_mapq=%d, min_length=%d, n_terminal=%d",
        min_mapq, min_length, n_terminal,
    )

    # ── Step 1: Damage profiling ─────────────────────────────────────────────
    logger.info("Step 1/5 — Profiling damage from BAM...")
    damage_profile = profile_damage(
        bam_path=bam_path,
        min_mapq=min_mapq,
        min_length=min_length,
        n_terminal=n_terminal,
    )
    logger.info(
        "  %d reads passed filters (of %d total).",
        damage_profile.n_reads_passed,
        damage_profile.n_reads_total,
    )

    # ── Step 2: Decay model fitting ─────────────────────────────────────────
    logger.info("Step 2/5 — Fitting geometric decay model...")
    model = fit_model(damage_profile)
    logger.info(
        "  5′ C→T: amplitude=%.4f, rate=%.4f, R²=%.4f (%s signal)",
        model.five_prime.amplitude, model.five_prime.rate,
        model.five_prime.r_squared, model.five_prime.signal_quality,
    )
    logger.info(
        "  3′ G→A: amplitude=%.4f, rate=%.4f, R²=%.4f (%s signal)",
        model.three_prime.amplitude, model.three_prime.rate,
        model.three_prime.r_squared, model.three_prime.signal_quality,
    )
    logger.info("  Library deamination rate: %.4f", model.library_deamination_rate)

    # ── Step 3: Bayesian classification ─────────────────────────────────────
    logger.info("Step 3/5 — Classifying reads...")
    classifications, class_summary = classify_reads(
        read_features=damage_profile.read_features,
        library_deamination_rate=model.library_deamination_rate,
        prior_ancient=prior_ancient,
        auth_threshold=auth_threshold,
        cont_threshold=cont_threshold,
    )
    logger.info(
        "  Authentic: %d (%.1f%%)  Contaminated: %d (%.1f%%)  Ambiguous: %d (%.1f%%)",
        class_summary.n_authentic,    class_summary.fraction_authentic * 100,
        class_summary.n_contaminated, class_summary.fraction_contaminated * 100,
        class_summary.n_ambiguous,    class_summary.fraction_ambiguous * 100,
    )

    # ── Step 4: Fragment lengths ─────────────────────────────────────────────
    logger.info("Step 4/5 — Computing fragment length statistics...")
    frag_stats = compute_fragment_lengths(damage_profile.read_features)
    logger.info(
        "  Mean: %.1f bp  Median: %.1f bp  Std: %.1f bp",
        frag_stats.mean, frag_stats.median, frag_stats.std,
    )

    # ── Step 5: Write output files ───────────────────────────────────────────
    logger.info("Step 5/5 — Writing output files...")

    args_dict = {
        "min_mapq": min_mapq,
        "min_length": min_length,
        "n_terminal": n_terminal,
        "prior_ancient": prior_ancient,
        "auth_threshold": auth_threshold,
        "cont_threshold": cont_threshold,
    }

    csv_path  = write_damage_csv(damage_profile, output_dir)
    tsv_path  = write_read_tsv(classifications, output_dir)
    json_path = write_summary_json(
        damage_profile, model, class_summary, frag_stats,
        output_dir, str(bam_path), args_dict,
    )
    txt_path  = write_summary_txt(damage_profile, model, class_summary, frag_stats, output_dir)

    output_files: dict[str, str] = {
        "damage_csv":    str(csv_path),
        "read_tsv":      str(tsv_path),
        "summary_json":  str(json_path),
        "summary_txt":   str(txt_path),
    }

    # ── Step 6: Plot (optional) ──────────────────────────────────────────────
    if no_plot:
        logger.info("Plot generation skipped (--no-plot).")
    else:
        try:
            png_path, svg_path = generate_damage_plot(
                profile=damage_profile,
                model=model,
                output_dir=output_dir,
                sample_name=sample_name,
            )
            output_files["damage_plot_png"] = str(png_path)
            output_files["damage_plot_svg"] = str(svg_path)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Plot generation failed: %s. Continuing without plot.", exc)

    logger.info("Done. Output directory: %s", output_dir)

    return {
        "profile":          damage_profile,
        "model":            model,
        "classifications":  classifications,
        "summary":          class_summary,
        "frag_stats":       frag_stats,
        "output_files":     output_files,
        "pipeline_version": __version__,
    }


def run_demo_pipeline(
    output_dir: Path,
    n_terminal: int = DEFAULT_N_TERMINAL,
    prior_ancient: float = DEFAULT_PRIOR_ANCIENT,
    auth_threshold: float = AUTH_THRESHOLD,
    cont_threshold: float = CONT_THRESHOLD,
    sample_name: str = "demo_synthetic",
    no_plot: bool = False,
) -> dict:
    """
    Execute the full pipeline using a synthetic ancient DNA damage profile.

    Generates a realistic simulated damage profile (no BAM file required) and
    runs all downstream steps — decay model fitting, Bayesian classification,
    fragment length analysis, output writing, and optionally plot generation.

    Args:
        output_dir:      Directory for all output files (created if absent).
        n_terminal:      Number of terminal positions to simulate and profile.
        prior_ancient:   Bayesian prior P(ancient) for the library.
        auth_threshold:  Posterior threshold for "authentic" label.
        cont_threshold:  Posterior threshold for "contaminated" label.
        sample_name:     Label for plot titles and reports.
        no_plot:         If True, skip matplotlib rendering.

    Returns:
        Structured result dict (same schema as run_pipeline()).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=== Ancient DNA Damage Classifier v%s  [DEMO MODE] ===", __version__)
    logger.info("Output dir   : %s", output_dir)
    logger.info("n_terminal=%d  prior_ancient=%.2f", n_terminal, prior_ancient)

    # ── Step 1: Generate synthetic damage profile ────────────────────────────
    logger.info("Step 1/5 — Generating synthetic damage profile...")
    damage_profile = generate_demo_profile(n_terminal=n_terminal)
    logger.info(
        "  %d synthetic reads generated.",
        damage_profile.n_reads_passed,
    )

    # ── Step 2: Decay model fitting ──────────────────────────────────────────
    logger.info("Step 2/5 — Fitting geometric decay model...")
    model = fit_model(damage_profile)
    logger.info(
        "  5′ C→T: amplitude=%.4f, rate=%.4f, R²=%.4f (%s signal)",
        model.five_prime.amplitude, model.five_prime.rate,
        model.five_prime.r_squared, model.five_prime.signal_quality,
    )
    logger.info(
        "  3′ G→A: amplitude=%.4f, rate=%.4f, R²=%.4f (%s signal)",
        model.three_prime.amplitude, model.three_prime.rate,
        model.three_prime.r_squared, model.three_prime.signal_quality,
    )
    logger.info("  Library deamination rate: %.4f", model.library_deamination_rate)

    # ── Step 3: Bayesian classification ─────────────────────────────────────
    logger.info("Step 3/5 — Classifying reads...")
    classifications, class_summary = classify_reads(
        read_features=damage_profile.read_features,
        library_deamination_rate=model.library_deamination_rate,
        prior_ancient=prior_ancient,
        auth_threshold=auth_threshold,
        cont_threshold=cont_threshold,
    )
    logger.info(
        "  Authentic: %d (%.1f%%)  Contaminated: %d (%.1f%%)  Ambiguous: %d (%.1f%%)",
        class_summary.n_authentic,    class_summary.fraction_authentic * 100,
        class_summary.n_contaminated, class_summary.fraction_contaminated * 100,
        class_summary.n_ambiguous,    class_summary.fraction_ambiguous * 100,
    )

    # ── Step 4: Fragment lengths ─────────────────────────────────────────────
    logger.info("Step 4/5 — Computing fragment length statistics...")
    frag_stats = compute_fragment_lengths(damage_profile.read_features)
    logger.info(
        "  Mean: %.1f bp  Median: %.1f bp  Std: %.1f bp",
        frag_stats.mean, frag_stats.median, frag_stats.std,
    )

    # ── Step 5: Write output files ───────────────────────────────────────────
    logger.info("Step 5/5 — Writing output files...")

    args_dict = {
        "demo": True,
        "n_terminal": n_terminal,
        "prior_ancient": prior_ancient,
        "auth_threshold": auth_threshold,
        "cont_threshold": cont_threshold,
    }

    csv_path  = write_damage_csv(damage_profile, output_dir)
    tsv_path  = write_read_tsv(classifications, output_dir)
    json_path = write_summary_json(
        damage_profile, model, class_summary, frag_stats,
        output_dir, "DEMO_SYNTHETIC", args_dict,
    )
    txt_path  = write_summary_txt(damage_profile, model, class_summary, frag_stats, output_dir)

    output_files: dict[str, str] = {
        "damage_csv":    str(csv_path),
        "read_tsv":      str(tsv_path),
        "summary_json":  str(json_path),
        "summary_txt":   str(txt_path),
    }

    # ── Step 6: Plot (optional) ──────────────────────────────────────────────
    if no_plot:
        logger.info("Plot generation skipped (--no-plot).")
    else:
        try:
            png_path, svg_path = generate_damage_plot(
                profile=damage_profile,
                model=model,
                output_dir=output_dir,
                sample_name=sample_name,
            )
            output_files["damage_plot_png"] = str(png_path)
            output_files["damage_plot_svg"] = str(svg_path)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Plot generation failed: %s. Continuing without plot.", exc)

    logger.info("Demo run complete. Output directory: %s", output_dir)

    return {
        "profile":          damage_profile,
        "model":            model,
        "classifications":  classifications,
        "summary":          class_summary,
        "frag_stats":       frag_stats,
        "output_files":     output_files,
        "pipeline_version": __version__,
    }


def main(argv: list[str] | None = None) -> int:
    """
    Entry point delegated to from main.py.

    Imports run_classifier.main() so that all CLI logic lives in one place.

    Args:
        argv: Argument list (defaults to sys.argv[1:] when None).

    Returns:
        Integer exit code.
    """
    # Import here to avoid circular import at module load time
    import importlib.util
    import os

    # run_classifier.py lives one level above src/
    rc_path = Path(__file__).parent.parent / "run_classifier.py"
    spec = importlib.util.spec_from_file_location("run_classifier", rc_path)
    if spec is None or spec.loader is None:
        print(f"ERROR: Cannot find run_classifier.py at {rc_path}", file=sys.stderr)
        return 1
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module.main(argv)
