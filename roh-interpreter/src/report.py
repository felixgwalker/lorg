"""Write ROH BED catalog, FROH CSV, and text / CSV analysis report."""

from __future__ import annotations

import datetime
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from src.froh_calculator import FROHResult
from src.roh_detector import ROHSegment

if TYPE_CHECKING:
    from src.ne_estimator import NeEstimate


def write_roh_bed(segments: list[ROHSegment], output_dir: Path) -> Path:
    """Write ROH segments as a tab-separated BED-like file."""
    rows = [
        {
            "chrom": s.chrom,
            "start": s.start_pos,
            "end": s.end_pos,
            "name": f"{s.individual_id}_{s.length_class}",
            "score": int(s.mean_homozygosity * 1000),
            "strand": ".",
            "individual_id": s.individual_id,
            "n_snps": s.n_snps,
            "length_bp": s.length_bp,
            "length_class": s.length_class,
            "mean_homozygosity": s.mean_homozygosity,
        }
        for s in segments
    ]
    df = pd.DataFrame(rows)
    out = Path(output_dir) / "roh_catalog.bed"
    df.to_csv(out, sep="\t", index=False)
    return out


def write_froh_csv(froh_results: list[FROHResult], output_dir: Path) -> Path:
    """Write FROH per-individual summary CSV."""
    rows = [
        {
            "individual_id": fr.individual_id,
            "n_roh": fr.n_roh,
            "total_roh_bp": fr.total_roh_bp,
            "froh": fr.froh,
            "froh_short": fr.froh_short,
            "froh_medium": fr.froh_medium,
            "froh_long": fr.froh_long,
            "n_short": fr.n_short,
            "n_medium": fr.n_medium,
            "n_long": fr.n_long,
            "bp_short": fr.bp_short,
            "bp_medium": fr.bp_medium,
            "bp_long": fr.bp_long,
        }
        for fr in froh_results
    ]
    df = pd.DataFrame(rows)
    out = Path(output_dir) / "froh_summary.csv"
    df.to_csv(out, index=False)
    return out


def write_roh_report(
    roh_segments: list[ROHSegment],
    froh: FROHResult | list[FROHResult],
    ne_estimate: "NeEstimate | list[NeEstimate] | None",
    output_dir: Path,
) -> tuple[Path, Path]:
    """Write a human-readable text report and a machine-readable CSV summary.

    Parameters
    ----------
    roh_segments:
        All ROH segments (one or more individuals).
    froh:
        A single :class:`~src.froh_calculator.FROHResult` or a list of them.
    ne_estimate:
        A single :class:`~src.ne_estimator.NeEstimate`, a list, or ``None``.
    output_dir:
        Directory for output files.  Created if absent.

    Returns
    -------
    (txt_path, csv_path) — Path objects for the text and CSV report files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Normalise to lists
    froh_list: list[FROHResult] = froh if isinstance(froh, list) else [froh]
    if ne_estimate is None:
        ne_list: list[NeEstimate] = []
    elif isinstance(ne_estimate, list):
        ne_list = ne_estimate
    else:
        ne_list = [ne_estimate]

    # --- Text report ---
    txt_path = output_dir / "roh_report.txt"
    lines: list[str] = [
        "=" * 70,
        "ROH Interpreter — Analysis Report",
        f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "=" * 70,
        "",
        f"Total ROH segments detected : {len(roh_segments)}",
    ]

    # Length-class breakdown
    for cls in ("short", "medium", "long"):
        n = sum(1 for s in roh_segments if s.length_class == cls)
        bp = sum(s.length_bp for s in roh_segments if s.length_class == cls)
        lines.append(f"  {cls:8s}: {n:4d} segments  ({bp / 1e6:,.2f} Mb total)")
    lines.append("")

    # Per-individual FROH summary
    lines.append("-" * 70)
    lines.append("FROH per individual:")
    lines.append(
        f"  {'ID':<20} {'FROH':>8} {'short':>9} {'medium':>9} {'long':>9}"
        f"  {'n_ROH':>6}"
    )
    lines.append("  " + "-" * 64)
    for fr in froh_list:
        lines.append(
            f"  {fr.individual_id:<20} {fr.froh:>8.4f}"
            f" {fr.froh_short:>9.4f} {fr.froh_medium:>9.4f} {fr.froh_long:>9.4f}"
            f"  {fr.n_roh:>6}"
        )
    lines.append("")

    # Ne estimates
    if ne_list:
        lines.append("-" * 70)
        lines.append("Effective population size (Ne) estimates:")
        lines.append(
            f"  {'ID':<20} {'Ne_recent':>10} {'Ne_moderate':>12}"
            f" {'Ne_ancient':>10} {'CI_low':>8} {'CI_high':>8}"
        )
        lines.append("  " + "-" * 72)
        for ne in ne_list:
            lines.append(
                f"  {ne.individual_id:<20} {ne.ne_recent:>10.1f}"
                f" {ne.ne_moderate:>12.1f} {ne.ne_ancient:>10.1f}"
                f" {ne.ci_low:>8.1f} {ne.ci_high:>8.1f}"
            )
        lines.append(
            f"  (generation time = {ne_list[0].generation_time_years} years)"
        )
        lines.append("")

    lines.append("=" * 70)

    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # --- CSV summary ---
    csv_rows = []
    froh_by_id = {fr.individual_id: fr for fr in froh_list}
    ne_by_id = {ne.individual_id: ne for ne in ne_list}
    all_ids = list({s.individual_id for s in roh_segments} | froh_by_id.keys())
    for ind_id in sorted(all_ids):
        fr = froh_by_id.get(ind_id)
        ne = ne_by_id.get(ind_id)
        row: dict = {"individual_id": ind_id}
        if fr is not None:
            row.update({
                "n_roh": fr.n_roh,
                "total_roh_bp": fr.total_roh_bp,
                "froh": fr.froh,
                "froh_short": fr.froh_short,
                "froh_medium": fr.froh_medium,
                "froh_long": fr.froh_long,
                "n_short": fr.n_short,
                "n_medium": fr.n_medium,
                "n_long": fr.n_long,
            })
        if ne is not None:
            row.update({
                "ne_recent": ne.ne_recent,
                "ne_moderate": ne.ne_moderate,
                "ne_ancient": ne.ne_ancient,
                "ne_ci_low": ne.ci_low,
                "ne_ci_high": ne.ci_high,
                "generation_time_years": ne.generation_time_years,
            })
        csv_rows.append(row)

    csv_path = output_dir / "roh_report.csv"
    pd.DataFrame(csv_rows).to_csv(csv_path, index=False)

    return txt_path, csv_path
