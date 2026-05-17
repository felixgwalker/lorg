"""
Output writers for the CNV Significance Assessor.

Produces four output files:

  cnv_annotated.csv          — full annotated CNV table
  significance_summary.txt   — human-readable classification summary
  significance_summary.json  — machine-readable summary
  gene_impact_report.csv     — all genes spanned by high-significance CNVs
"""

from __future__ import annotations

import csv
import json
import logging
from datetime import datetime
from pathlib import Path

from src import __version__
from src.classifier import ClassificationSummary, ScoredCNV, TIER_PATHOGENIC, TIER_VUS
from src.config import (
    OUTFILE_ANNOTATED_CSV,
    OUTFILE_GENE_IMPACT_CSV,
    OUTFILE_SIGNIFICANCE_JSON,
    OUTFILE_SIGNIFICANCE_TXT,
)

logger = logging.getLogger(__name__)


# ── Annotated CNV table ───────────────────────────────────────────────────

def write_annotated_csv(scored: list[ScoredCNV], output_dir: Path) -> Path:
    """
    Write the full annotated CNV table to CSV.

    Columns:
        cnv_id, chrom, start, end, cnv_type, size_bp,
        n_genes, gene_names, n_regulatory,
        dosage_metric, max_dosage_sensitivity,
        pop_frequency, pop_match_count,
        size_score, gene_score, dosage_score, pop_modifier, total_score,
        significance_tier, classification_reason

    Args:
        scored:     List of ScoredCNV from classifier.classify_cnvs().
        output_dir: Directory to write into.

    Returns:
        Path to the written CSV file.
    """
    out_path = output_dir / OUTFILE_ANNOTATED_CSV

    fieldnames = [
        "cnv_id", "chrom", "start", "end", "cnv_type", "size_bp",
        "n_genes", "gene_names", "n_regulatory",
        "dosage_metric", "max_dosage_sensitivity",
        "pop_frequency", "pop_match_count",
        "size_score", "gene_score", "dosage_score", "pop_modifier", "total_score",
        "significance_tier", "classification_reason",
    ]

    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for s in scored:
            rec = s.annotated.record
            ann = s.annotated
            pop_freq_str = f"{ann.pop_frequency:.6g}" if ann.pop_frequency is not None else ""
            writer.writerow({
                "cnv_id":                rec.cnv_id,
                "chrom":                 rec.chrom,
                "start":                 rec.start,
                "end":                   rec.end,
                "cnv_type":              rec.cnv_type,
                "size_bp":               rec.size,
                "n_genes":               ann.n_genes,
                "gene_names":            ";".join(ann.gene_names),
                "n_regulatory":          len(ann.overlapping_regulatory),
                "dosage_metric":         ann.dosage_metric,
                "max_dosage_sensitivity": f"{ann.max_dosage_sensitivity:.4f}",
                "pop_frequency":         pop_freq_str,
                "pop_match_count":       ann.pop_match_count,
                "size_score":            s.size_score,
                "gene_score":            s.gene_score,
                "dosage_score":          s.dosage_score,
                "pop_modifier":          s.pop_modifier,
                "total_score":           s.total_score,
                "significance_tier":     s.significance_tier,
                "classification_reason": s.classification_reason,
            })

    logger.info("Annotated CNV table written: %s", out_path)
    return out_path


# ── Gene impact report ────────────────────────────────────────────────────

def write_gene_impact_report(scored: list[ScoredCNV], output_dir: Path) -> Path:
    """
    Write a CSV listing all genes spanned by VUS or LIKELY_PATHOGENIC CNVs.

    Columns:
        gene_name, gene_id, chrom, gene_start, gene_end, strand,
        cnv_id, cnv_type, cnv_chrom, cnv_start, cnv_end,
        significance_tier, max_dosage_sensitivity, dosage_metric

    Genes may appear multiple times if spanned by more than one high-significance CNV.

    Args:
        scored:     List of ScoredCNV.
        output_dir: Output directory.

    Returns:
        Path to the written CSV.
    """
    out_path = output_dir / OUTFILE_GENE_IMPACT_CSV

    fieldnames = [
        "gene_name", "gene_id", "chrom", "gene_start", "gene_end", "strand",
        "cnv_id", "cnv_type", "cnv_start", "cnv_end",
        "significance_tier", "max_dosage_sensitivity", "dosage_metric",
    ]

    high_sig_tiers = {TIER_VUS, TIER_PATHOGENIC}

    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for s in scored:
            if s.significance_tier not in high_sig_tiers:
                continue
            rec = s.annotated.record
            ann = s.annotated
            for gene in ann.overlapping_genes:
                writer.writerow({
                    "gene_name":              gene.gene_name,
                    "gene_id":                gene.gene_id,
                    "chrom":                  gene.chrom,
                    "gene_start":             gene.start,
                    "gene_end":               gene.end,
                    "strand":                 gene.strand,
                    "cnv_id":                 rec.cnv_id,
                    "cnv_type":               rec.cnv_type,
                    "cnv_start":              rec.start,
                    "cnv_end":                rec.end,
                    "significance_tier":      s.significance_tier,
                    "max_dosage_sensitivity": f"{ann.max_dosage_sensitivity:.4f}",
                    "dosage_metric":          ann.dosage_metric,
                })

    n_rows = sum(
        len(s.annotated.overlapping_genes)
        for s in scored if s.significance_tier in high_sig_tiers
    )
    logger.info("Gene impact report written: %s  (%d gene×CNV rows)", out_path, n_rows)
    return out_path


# ── Significance summary (TXT) ────────────────────────────────────────────

def write_significance_txt(
    summary: ClassificationSummary,
    scored: list[ScoredCNV],
    output_dir: Path,
    input_cnv_path: str = "",
    args_dict: dict | None = None,
) -> Path:
    """
    Write a human-readable significance summary to plain text.

    Args:
        summary:        ClassificationSummary from classifier.
        scored:         List of ScoredCNV.
        output_dir:     Output directory.
        input_cnv_path: String path of the CNV input (for metadata header).
        args_dict:      CLI argument dict for reproducibility notes.

    Returns:
        Path to the written TXT file.
    """
    out_path = output_dir / OUTFILE_SIGNIFICANCE_TXT
    lines = _build_summary_lines(summary, scored, input_cnv_path, args_dict)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    logger.info("Significance summary (TXT) written: %s", out_path)
    return out_path


def _build_summary_lines(
    summary: ClassificationSummary,
    scored: list[ScoredCNV],
    input_cnv_path: str = "",
    args_dict: dict | None = None,
) -> list[str]:
    """Build the list of lines for the human-readable summary."""
    sep = "=" * 60
    lines = [
        "",
        sep,
        "     CNV SIGNIFICANCE ASSESSOR  v" + __version__,
        sep,
        "",
    ]

    if input_cnv_path:
        lines.append(f"  Input             : {input_cnv_path}")
    lines.append(f"  Analysis date     : {datetime.now().isoformat()[:10]}")
    lines.append("")

    lines += [
        "  ── Classification Summary ─────────────────────────────",
        "",
        f"  {'Tier':<25}  {'Count':>7}  {'Fraction':>9}",
        f"  {'-'*25}  {'-'*7}  {'-'*9}",
        f"  {'LIKELY_BENIGN':<25}  {summary.n_likely_benign:>7,}  {summary.fraction_benign:>8.1%}",
        f"  {'VUS':<25}  {summary.n_vus:>7,}  {summary.fraction_vus:>8.1%}",
        f"  {'LIKELY_PATHOGENIC':<25}  {summary.n_likely_pathogenic:>7,}  {summary.fraction_pathogenic:>8.1%}",
        f"  {'-'*25}  {'-'*7}",
        f"  {'Total':<25}  {summary.n_total:>7,}",
        "",
    ]

    # CNV type breakdown
    type_counts: dict[str, dict[str, int]] = {}
    for s in scored:
        ctype = s.annotated.record.cnv_type
        tier  = s.significance_tier
        type_counts.setdefault(ctype, {}).setdefault(tier, 0)
        type_counts[ctype][tier] += 1

    if type_counts:
        lines += [
            "  ── By CNV Type ────────────────────────────────────────",
            "",
            f"  {'Type':<10}  {'Benign':>7}  {'VUS':>7}  {'Pathogenic':>10}",
            f"  {'-'*10}  {'-'*7}  {'-'*7}  {'-'*10}",
        ]
        for ctype in sorted(type_counts):
            tc = type_counts[ctype]
            lines.append(
                f"  {ctype:<10}  "
                f"{tc.get('LIKELY_BENIGN', 0):>7}  "
                f"{tc.get('VUS', 0):>7}  "
                f"{tc.get('LIKELY_PATHOGENIC', 0):>10}"
            )
        lines.append("")

    # Top pathogenic/VUS CNVs
    high_sig = [
        s for s in scored
        if s.significance_tier in {TIER_PATHOGENIC, TIER_VUS}
    ]
    high_sig.sort(key=lambda s: (-s.total_score, -s.annotated.record.size))

    if high_sig:
        n_show = min(10, len(high_sig))
        lines += [
            f"  ── Top {n_show} High-Significance CNVs ───────────────────────",
            "",
            f"  {'ID':<20}  {'Tier':<18}  {'Score':>5}  {'Genes':>5}  {'MaxDS':>5}",
            f"  {'-'*20}  {'-'*18}  {'-'*5}  {'-'*5}  {'-'*5}",
        ]
        for s in high_sig[:n_show]:
            rec = s.annotated.record
            lines.append(
                f"  {rec.cnv_id[:20]:<20}  "
                f"{s.significance_tier:<18}  "
                f"{s.total_score:>5}  "
                f"{s.annotated.n_genes:>5}  "
                f"{s.annotated.max_dosage_sensitivity:>5.2f}"
            )
        lines.append("")

    lines += [sep, ""]
    return lines


# ── Significance summary (JSON) ───────────────────────────────────────────

def write_significance_json(
    summary: ClassificationSummary,
    scored: list[ScoredCNV],
    output_dir: Path,
    input_cnv_path: str = "",
    args_dict: dict | None = None,
) -> Path:
    """
    Write a machine-readable significance summary to JSON.

    Args:
        summary:        ClassificationSummary from classifier.
        scored:         List of ScoredCNV.
        output_dir:     Output directory.
        input_cnv_path: String path of the CNV input (for metadata).
        args_dict:      CLI argument dict for reproducibility.

    Returns:
        Path to the written JSON file.
    """
    out_path = output_dir / OUTFILE_SIGNIFICANCE_JSON

    # Tier breakdowns by CNV type
    type_breakdown: dict[str, dict[str, int]] = {}
    for s in scored:
        ctype = s.annotated.record.cnv_type
        tier  = s.significance_tier
        type_breakdown.setdefault(ctype, {"LIKELY_BENIGN": 0, "VUS": 0, "LIKELY_PATHOGENIC": 0})
        type_breakdown[ctype][tier] = type_breakdown[ctype].get(tier, 0) + 1

    report = {
        "meta": {
            "tool":          "CNV Significance Assessor",
            "version":       __version__,
            "analysis_date": datetime.now().isoformat()[:10],
            "input_cnv":     input_cnv_path,
            "parameters":    args_dict or {},
        },
        "classification_summary": {
            "n_total":              summary.n_total,
            "n_likely_benign":      summary.n_likely_benign,
            "n_vus":                summary.n_vus,
            "n_likely_pathogenic":  summary.n_likely_pathogenic,
            "fraction_benign":      summary.fraction_benign,
            "fraction_vus":         summary.fraction_vus,
            "fraction_pathogenic":  summary.fraction_pathogenic,
        },
        "by_cnv_type":  type_breakdown,
        "score_statistics": {
            "mean_total_score": round(
                sum(s.total_score for s in scored) / len(scored), 3
            ) if scored else 0.0,
            "max_total_score": max((s.total_score for s in scored), default=0),
            "n_with_gene_overlap": sum(1 for s in scored if s.annotated.n_genes > 0),
            "n_with_pop_frequency": sum(
                1 for s in scored if s.annotated.pop_frequency is not None
            ),
            "n_with_high_ds": sum(
                1 for s in scored
                if s.annotated.max_dosage_sensitivity >= 0.9
            ),
        },
    }

    out_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    logger.info("Significance summary (JSON) written: %s", out_path)
    return out_path
