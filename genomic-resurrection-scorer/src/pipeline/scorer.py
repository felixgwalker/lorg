"""
Genomic Resurrection Scorer — pipeline orchestrator.

Reads a metrics JSON file, runs all five scoring layers in sequence, combines
sub-scores into a weighted feasibility index, and returns a structured report.

Usage (run from project root):
    python run_scorer.py data/thylacine_case_study/metrics.json
    python run_scorer.py data/thylacine_case_study/metrics.json --output src/frontend/public/reports/report.json
"""

import json
import argparse
import sys
from datetime import datetime
from pathlib import Path

from .config import LAYER_WEIGHTS, GRADE_THRESHOLDS
from .layers.ancient_dna_quality import score_adq
from .layers.genomic_completeness import score_gc
from .layers.divergence import score_divergence
from .layers.edit_burden import score_edit_burden
from .layers.ethical_ecological import score_ethical_ecological


def load_metrics(path: str) -> dict:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def assign_grade(score: float) -> tuple[str, str]:
    for threshold, grade, label in GRADE_THRESHOLDS:
        if score >= threshold:
            return grade, label
    return "F", "Currently Infeasible"


def calculate_feasibility_index(layer_results: dict[str, dict]) -> float:
    return round(
        sum(
            layer_results[layer]["score"] * LAYER_WEIGHTS[layer]
            for layer in LAYER_WEIGHTS
        ),
        1,
    )


def run_scorer(metrics_path: str, output_path: str | None = None) -> dict:
    metrics = load_metrics(metrics_path)

    layer_results: dict[str, dict] = {
        "ancient_dna_quality":  score_adq(metrics["ancient_dna_quality"]),
        "genomic_completeness": score_gc(metrics["genomic_completeness"]),
        "divergence":           score_divergence(metrics["divergence"]),
        "edit_burden":          score_edit_burden(metrics["edit_burden"]),
        "ethical_ecological":   score_ethical_ecological(metrics["ethical_ecological"]),
    }

    feasibility_index = calculate_feasibility_index(layer_results)
    grade, grade_label = assign_grade(feasibility_index)

    report: dict = {
        "meta": {
            **metrics.get("meta", {}),
            "analysis_date":  datetime.now().isoformat()[:10],
            "scorer_version": "1.0.0",
        },
        "overall": {
            "feasibility_index": feasibility_index,
            "grade":             grade,
            "grade_label":       grade_label,
            "layer_weights":     LAYER_WEIGHTS,
        },
        "layers": {
            layer: {**result, "weight": LAYER_WEIGHTS[layer]}
            for layer, result in layer_results.items()
        },
    }

    if output_path:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"Report written -> {out}")

    return report


def print_summary(report: dict) -> None:
    meta    = report["meta"]
    overall = report["overall"]
    layers  = report["layers"]

    LAYER_LABELS = {
        "ancient_dna_quality":  "Ancient DNA Quality",
        "genomic_completeness": "Genomic Completeness",
        "divergence":           "Divergence",
        "edit_burden":          "Edit Burden",
        "ethical_ecological":   "Ethical / Ecological",
    }

    print()
    print("=" * 54)
    print("        GENOMIC RESURRECTION SCORER  v1.0.0")
    print("=" * 54)
    print(f"  Extinct species : {meta.get('common_name_extinct', meta.get('species_extinct', '—'))}")
    print(f"  Proxy species   : {meta.get('common_name_proxy',   meta.get('species_proxy',   '—'))}")
    print(f"  Analysis date   : {meta.get('analysis_date', '—')}")
    print()
    print(f"  FEASIBILITY INDEX : {overall['feasibility_index']:.1f} / 100")
    print(f"  Grade             : {overall['grade']} - {overall['grade_label']}")
    print()
    print(f"  {'Layer':<30}  {'Score':>6}  {'Grade':>5}  {'Weight':>6}")
    print(f"  {'-'*30}  {'-'*6}  {'-'*5}  {'-'*6}")
    for key, data in layers.items():
        label  = LAYER_LABELS.get(key, key)
        weight = int(data["weight"] * 100)
        print(f"  {label:<30}  {data['score']:>5.1f}   {data['grade']:>4}   {weight:>4}%")
    print()
    for key, data in layers.items():
        if data.get("flags"):
            print(f"  [!] {LAYER_LABELS.get(key, key)} flags:")
            for flag in data["flags"]:
                print(f"       - {flag}")
    print()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Genomic Resurrection Scorer — feasibility assessment pipeline"
    )
    parser.add_argument("metrics", help="Path to input metrics JSON file")
    parser.add_argument("--output", "-o", metavar="FILE",
                        help="Write full report to this JSON file")
    parser.add_argument("--quiet", "-q", action="store_true",
                        help="Suppress console summary output")
    args = parser.parse_args(argv)

    report = run_scorer(args.metrics, args.output)
    if not args.quiet:
        print_summary(report)

    return 0


if __name__ == "__main__":
    sys.exit(main())
