"""
Output file writer for the Ancient DNA Damage Classifier.

Writes five output files:
  - damage_frequencies.csv   — per-position C→T and G→A rates and counts
  - read_classifications.tsv — per-read Bayesian classification results
  - summary_report.json      — machine-readable structured report
  - summary_report.txt       — human-readable console-style summary
"""

from __future__ import annotations

import csv
import json
import logging
from datetime import datetime
from pathlib import Path

import numpy as np

from src import __version__
from src.classifier import ClassificationSummary, ReadClassification
from src.config import (
    OUTFILE_DAMAGE_CSV,
    OUTFILE_READS_TSV,
    OUTFILE_SUMMARY_JSON,
    OUTFILE_SUMMARY_TXT,
)
from src.damage_profiler import DamageProfile
from src.decay_model import ModelResult
from src.fragment_length import FragmentLengthStats

logger = logging.getLogger(__name__)


def write_damage_csv(profile: DamageProfile, output_dir: Path) -> Path:
    """
    Write per-position damage frequency table to CSV.

    Columns: position, ct_rate, ct_count, c_count, ga_rate, ga_count, g_count.
    Positions are 1-based; index 1 = most terminal.  The C→T columns reflect 5'
    terminal positions; the G→A columns reflect the corresponding 3' terminal
    positions (index 1 = most 3'-terminal).

    Args:
        profile:    DamageProfile with per-position arrays.
        output_dir: Directory to write into.

    Returns:
        Path to the written CSV file.
    """
    out_path = output_dir / OUTFILE_DAMAGE_CSV
    ct_rate = profile.ct_rate
    ga_rate = profile.ga_rate

    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "position",
            "ct_rate", "ct_count", "c_count",
            "ga_rate", "ga_count", "g_count",
        ])
        for i in range(profile.n_terminal):
            writer.writerow([
                i + 1,
                f"{ct_rate[i]:.6f}", int(profile.ct_count[i]), int(profile.c_count[i]),
                f"{ga_rate[i]:.6f}", int(profile.ga_count[i]), int(profile.g_count[i]),
            ])

    logger.info("Damage frequency table written: %s", out_path)
    return out_path


def write_read_tsv(
    classifications: list[ReadClassification],
    output_dir: Path,
) -> Path:
    """
    Write per-read classification table to TSV.

    Columns: read_id, classification, posterior_ancient, posterior_contaminated,
             ct_terminal, ga_terminal, read_length, uninformative.

    Args:
        classifications: List from classifier.classify_reads().
        output_dir:      Directory to write into.

    Returns:
        Path to the written TSV file.
    """
    out_path = output_dir / OUTFILE_READS_TSV

    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "read_id", "classification",
            "posterior_ancient", "posterior_contaminated",
            "ct_terminal", "ga_terminal",
            "read_length", "uninformative",
        ])
        for rc in classifications:
            writer.writerow([
                rc.read_id,
                rc.classification,
                f"{rc.posterior_ancient:.6f}",
                f"{rc.posterior_contaminated:.6f}",
                rc.ct_terminal,
                rc.ga_terminal,
                rc.read_length,
                "true" if rc.uninformative else "false",
            ])

    logger.info("Read classifications written: %s", out_path)
    return out_path


def write_summary_json(
    profile: DamageProfile,
    model: ModelResult,
    summary: ClassificationSummary,
    frag_stats: FragmentLengthStats,
    output_dir: Path,
    bam_path: str,
    args_dict: dict,
) -> Path:
    """
    Write machine-readable summary statistics to JSON.

    Args:
        profile:     DamageProfile from damage_profiler.
        model:       ModelResult from decay_model.
        summary:     ClassificationSummary from classifier.
        frag_stats:  FragmentLengthStats from fragment_length.
        output_dir:  Directory to write into.
        bam_path:    String path of the input BAM (for metadata).
        args_dict:   Dict of CLI arguments for reproducibility.

    Returns:
        Path to the written JSON file.
    """
    out_path = output_dir / OUTFILE_SUMMARY_JSON

    report = {
        "meta": {
            "tool": "Ancient DNA Damage Classifier",
            "version": __version__,
            "analysis_date": datetime.now().isoformat()[:10],
            "bam_input": bam_path,
            "parameters": args_dict,
        },
        "read_stats": {
            "n_reads_total": profile.n_reads_total,
            "n_reads_passed": profile.n_reads_passed,
            "n_reads_no_md": profile.n_reads_no_md,
            "pass_rate": (
                round(profile.n_reads_passed / profile.n_reads_total, 4)
                if profile.n_reads_total > 0 else 0.0
            ),
        },
        "damage_model": {
            "five_prime": {
                "amplitude": round(model.five_prime.amplitude, 6),
                "rate": round(model.five_prime.rate, 6),
                "background": round(model.five_prime.background, 6),
                "r_squared": round(model.five_prime.r_squared, 4),
                "converged": model.five_prime.converged,
                "signal_quality": model.five_prime.signal_quality,
            },
            "three_prime": {
                "amplitude": round(model.three_prime.amplitude, 6),
                "rate": round(model.three_prime.rate, 6),
                "background": round(model.three_prime.background, 6),
                "r_squared": round(model.three_prime.r_squared, 4),
                "converged": model.three_prime.converged,
                "signal_quality": model.three_prime.signal_quality,
            },
            "library_deamination_rate": round(model.library_deamination_rate, 6),
            "overall_signal_quality": model.overall_signal_quality,
        },
        "classification": {
            "n_authentic": summary.n_authentic,
            "n_contaminated": summary.n_contaminated,
            "n_ambiguous": summary.n_ambiguous,
            "n_total": summary.n_total,
            "fraction_authentic": round(summary.fraction_authentic, 4),
            "fraction_contaminated": round(summary.fraction_contaminated, 4),
            "fraction_ambiguous": round(summary.fraction_ambiguous, 4),
            "mean_posterior_ancient": round(summary.mean_posterior_ancient, 4),
            "overall_authenticity_estimate": round(summary.fraction_authentic, 4),
        },
        "fragment_length": {
            "mean": round(frag_stats.mean, 2),
            "median": round(frag_stats.median, 2),
            "std": round(frag_stats.std, 2),
            "min": frag_stats.min_len,
            "max": frag_stats.max_len,
            "n_reads": frag_stats.n_reads,
            "is_paired_end": frag_stats.is_paired,
        },
        "damage_rates_terminal_position_1": {
            "ct_rate_5prime": round(float(profile.ct_rate[0]), 6) if profile.n_terminal > 0 else None,
            "ga_rate_3prime": round(float(profile.ga_rate[0]), 6) if profile.n_terminal > 0 else None,
        },
    }

    out_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    logger.info("Summary JSON written: %s", out_path)
    return out_path


def write_summary_txt(
    profile: DamageProfile,
    model: ModelResult,
    summary: ClassificationSummary,
    frag_stats: FragmentLengthStats,
    output_dir: Path,
) -> Path:
    """
    Write human-readable summary report to plain text.

    Args:
        profile:     DamageProfile from damage_profiler.
        model:       ModelResult from decay_model.
        summary:     ClassificationSummary from classifier.
        frag_stats:  FragmentLengthStats from fragment_length.
        output_dir:  Directory to write into.

    Returns:
        Path to the written TXT file.
    """
    out_path = output_dir / OUTFILE_SUMMARY_TXT
    lines = _build_summary_lines(profile, model, summary, frag_stats)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    logger.info("Summary TXT written: %s", out_path)
    return out_path


def _build_summary_lines(
    profile: DamageProfile,
    model: ModelResult,
    summary: ClassificationSummary,
    frag_stats: FragmentLengthStats,
) -> list[str]:
    """Build the list of lines for the human-readable summary."""
    sep = "=" * 58
    lines = [
        "",
        sep,
        "     ANCIENT DNA DAMAGE CLASSIFIER  v" + __version__,
        sep,
        "",
        f"  Reads total          : {profile.n_reads_total:>10,}",
        f"  Reads passed filters : {profile.n_reads_passed:>10,}",
        f"  Reads (no MD tags)   : {profile.n_reads_no_md:>10,}",
        "",
        "  ── Damage Signal ──────────────────────────────────",
        "",
        f"  5′ C→T  amplitude    : {model.five_prime.amplitude:.4f}",
        f"  5′ C→T  decay rate   : {model.five_prime.rate:.4f}",
        f"  5′ C→T  R²           : {model.five_prime.r_squared:.4f}  ({model.five_prime.signal_quality})",
        "",
        f"  3′ G→A  amplitude    : {model.three_prime.amplitude:.4f}",
        f"  3′ G→A  decay rate   : {model.three_prime.rate:.4f}",
        f"  3′ G→A  R²           : {model.three_prime.r_squared:.4f}  ({model.three_prime.signal_quality})",
        "",
        f"  Library deamination  : {model.library_deamination_rate:.4f}",
        f"  Overall signal       : {model.overall_signal_quality.upper()}",
    ]

    if profile.n_terminal > 0:
        lines += [
            "",
            f"  5′ C→T at position 1 : {profile.ct_rate[0]:.4f}",
            f"  3′ G→A at position 1 : {profile.ga_rate[0]:.4f}",
        ]

    lines += [
        "",
        "  ── Read Classification ────────────────────────────",
        "",
        f"  {'Label':<16}  {'Count':>8}  {'Fraction':>9}",
        f"  {'-'*16}  {'-'*8}  {'-'*9}",
        f"  {'Authentic':<16}  {summary.n_authentic:>8,}  {summary.fraction_authentic:>8.1%}",
        f"  {'Contaminated':<16}  {summary.n_contaminated:>8,}  {summary.fraction_contaminated:>8.1%}",
        f"  {'Ambiguous':<16}  {summary.n_ambiguous:>8,}  {summary.fraction_ambiguous:>8.1%}",
        f"  {'Total':<16}  {summary.n_total:>8,}",
        "",
        f"  Mean posterior P(ancient)  : {summary.mean_posterior_ancient:.4f}",
        f"  Authenticity estimate      : {summary.fraction_authentic:.1%}",
        "",
        "  ── Fragment Length ────────────────────────────────",
        "",
        f"  Mean length   : {frag_stats.mean:.1f} bp",
        f"  Median length : {frag_stats.median:.1f} bp",
        f"  Std deviation : {frag_stats.std:.1f} bp",
        f"  Range         : {frag_stats.min_len} – {frag_stats.max_len} bp",
        f"  Data type     : {'Paired-end' if frag_stats.is_paired else 'Single-end'}",
        "",
        sep,
        "",
    ]

    return lines
