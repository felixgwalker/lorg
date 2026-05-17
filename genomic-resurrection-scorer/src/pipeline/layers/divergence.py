"""
Layer 3: Divergence Quantification

Measures genome-wide and functional-region sequence divergence between the
extinct species and its living proxy. Higher divergence means more editing is
required and raises the biological uncertainty around gene-function transfer.

The score is *inversely* proportional to divergence: a score of 100 indicates
identical sequences; a score near 0 indicates extreme divergence.

Calibration reference points (genome-wide SNP rate → expected score):
  Mammoth / Asian elephant  ~0.6%  →  ~94
  Human / Chimpanzee        ~1.2%  →  ~88
  Human / Orangutan         ~3.1%  →  ~71
  Thylacine / T. devil     ~11.8%  →  ~27

Score range: 0–100 (higher = lower divergence / more feasible).
"""

import math
from ..config import DIV_WEIGHTS

# Scaling constant for the exponential decay function.
# Chosen so that human/chimp (~1.2%) → ~88 and mammoth/elephant (~0.6%) → ~94.
_DECAY_K = 11.0


def _divergence_score(snp_rate: float) -> float:
    """Convert SNP rate (0–1) to 0–100 feasibility score via exponential decay."""
    return round(max(5.0, 100.0 * math.exp(-_DECAY_K * snp_rate)), 1)


def score_divergence(metrics: dict) -> dict:
    """
    Compute the Divergence layer score.

    Expected metrics keys:
        coding_snp_rate      (float) SNP rate in annotated coding regions
        regulatory_snp_rate  (float) SNP rate in annotated regulatory regions
        genome_wide_snp_rate (float) SNP rate across the whole genome
        total_snps           (int)   Total SNPs called
        total_indels         (int)   Total indels called
        total_svs            (int)   Total structural variants called
        mya_divergence       (float) Estimated divergence time (million years ago)
    """
    components = {
        "coding_divergence":      _divergence_score(metrics["coding_snp_rate"]),
        "regulatory_divergence":  _divergence_score(metrics["regulatory_snp_rate"]),
        "genome_wide_divergence": _divergence_score(metrics["genome_wide_snp_rate"]),
    }

    score = round(sum(components[k] * DIV_WEIGHTS[k] for k in components), 1)

    flags: list[str] = []
    if metrics["genome_wide_snp_rate"] > 0.05:
        flags.append("HIGH_GENOME_WIDE_DIVERGENCE")
    if metrics["coding_snp_rate"] > 0.05:
        flags.append("HIGH_CODING_DIVERGENCE")
    if metrics["regulatory_snp_rate"] > 0.08:
        flags.append("HIGH_REGULATORY_DIVERGENCE")
    if metrics["total_svs"] > 5000:
        flags.append("HIGH_STRUCTURAL_VARIANT_COUNT")

    return {
        "score": score,
        "grade": _grade(score),
        "components": components,
        "interpretation": _interpret(score, metrics),
        "flags": flags,
        "context": {
            "total_snps":    metrics["total_snps"],
            "total_indels":  metrics["total_indels"],
            "total_svs":     metrics["total_svs"],
            "mya_divergence": metrics["mya_divergence"],
        },
    }


def _interpret(score: float, m: dict) -> str:
    gw = m["genome_wide_snp_rate"] * 100
    coding = m["coding_snp_rate"] * 100
    reg = m["regulatory_snp_rate"] * 100
    mya = m["mya_divergence"]
    snps = m["total_snps"]

    if score >= 80:
        return (
            f"Low divergence from proxy species ({gw:.2f}% genome-wide). "
            f"The extinct and proxy genomes are highly similar, making the proxy "
            f"a strong functional stand-in for most gene networks."
        )
    elif score >= 55:
        return (
            f"Moderate divergence ({gw:.2f}% genome-wide, diverged ~{mya:.1f} Mya). "
            f"Coding regions ({coding:.2f}%) and regulatory elements ({reg:.2f}%) "
            f"show meaningful but tractable differences. Targeted editing of key "
            f"functional loci is the primary strategy."
        )
    elif score >= 35:
        return (
            f"High divergence ({gw:.2f}% genome-wide, diverged ~{mya:.1f} Mya). "
            f"{snps:,} total SNPs identified. Coding ({coding:.2f}%) and regulatory "
            f"({reg:.2f}%) divergence pose substantial challenges to faithful "
            f"gene-network reconstruction."
        )
    else:
        return (
            f"Extremely high divergence ({gw:.2f}% genome-wide, ~{mya:.1f} Mya). "
            f"This far exceeds well-studied de-extinction pairs such as "
            f"mammoth/elephant (~0.6%). The proxy genome is a distant structural "
            f"template only; functional equivalence cannot be assumed for most loci. "
            f"This is the primary biological constraint on feasibility."
        )


def _grade(score: float) -> str:
    if score >= 80: return "A"
    if score >= 65: return "B"
    if score >= 50: return "C"
    if score >= 35: return "D"
    return "F"
