"""
Natural language explanation generator for the three-layer HGT risk result.

Produces a single cautious, scientifically worded paragraph summarising what
drove the score, what reduces it, and what data was missing.  The language
avoids deterministic claims: 'may', 'suggests', 'consistent with', 'cannot
be excluded', etc.

The generator is purely rule-based: it reads feature scores and selects
sentence fragments from predefined templates.  This keeps the output
interpretable and auditable.
"""

from __future__ import annotations

from src.models import LayerResult, ScoreBand

# Thresholds for classifying a feature's contribution
_HIGH_THRESHOLD   = 0.60    # score >= this → positive contributor
_LOW_THRESHOLD    = 0.20    # available score <= this → risk reducer
_IMPORTANT_UNAVAILABLE = {  # features whose absence warrants a comment
    "amr_content", "virulence_flags", "plasmid_context",
    "transposase_proximity", "conjugative_element", "is_element_match",
}

_FEATURE_LABELS: dict[str, str] = {
    "is_element_match":      "mobile-element adjacency",
    "integron_association":  "integron/attC-site association",
    "conjugative_element":   "conjugative element homology",
    "plasmid_context":       "plasmid-like context",
    "transposase_proximity": "transposase/integrase proximity",
    "repeat_density":        "elevated repeat density",
    "gc_deviation":          "GC content deviation from host",
    "codon_usage_distance":  "codon usage divergence from host",
    "taxonomic_distance":    "large donor–recipient taxonomic distance",
    "promoter_plausibility": "prokaryotic promoter-like sequences",
    "sequence_complexity":   "multi-gene sequence complexity",
    "prophage_context":      "prophage-associated context",
    "amr_content":           "AMR-associated content",
    "virulence_flags":       "virulence/toxin-associated annotations",
    "gene_completeness":     "high proportion of complete ORFs",
    "payload_count":         "multiple candidate functional genes",
}

_BAND_OPENER: dict[str, str] = {
    "low":      "The sequence shows minimal sequence-level indicators of HGT risk.",
    "moderate": "The sequence shows some sequence-level indicators consistent with elevated HGT potential.",
    "high":     "The sequence shows multiple significant sequence-level indicators of HGT risk.",
    "very_high":"The sequence shows strong sequence-level indicators of HGT risk across multiple assessment layers.",
}


def generate_explanation(
    transfer_layer: LayerResult,
    establishment_layer: LayerResult,
    consequence_layer: LayerResult,
    hgt_risk_index: float,
    score_band: ScoreBand,
    overall_completeness: float,
) -> tuple[str, list[str], list[str], list[str]]:
    """
    Generate a natural language explanation for the three-layer result.

    Returns
    -------
    explanation              : full paragraph explanation string
    top_contributors         : feature labels with score >= HIGH_THRESHOLD
    risk_reducers            : available feature labels with score <= LOW_THRESHOLD
    missing_important_features : important features that were unavailable
    """
    all_layers = [transfer_layer, establishment_layer, consequence_layer]
    all_features = [f for lr in all_layers for f in lr.feature_results]

    # Identify key groups
    top_contributors: list[str] = [
        _FEATURE_LABELS.get(f.feature_name, f.feature_name)
        for f in all_features
        if f.available and f.score is not None and f.score >= _HIGH_THRESHOLD
    ]
    risk_reducers: list[str] = [
        _FEATURE_LABELS.get(f.feature_name, f.feature_name)
        for f in all_features
        if f.available and f.score is not None and f.score <= _LOW_THRESHOLD
    ]
    missing_important: list[str] = [
        _FEATURE_LABELS.get(f.feature_name, f.feature_name)
        for f in all_features
        if not f.available and f.feature_name in _IMPORTANT_UNAVAILABLE
    ]

    # Build sentence fragments
    parts: list[str] = [_BAND_OPENER.get(score_band.value, "")]

    # Layer-level commentary
    layer_comments: list[str] = []
    if transfer_layer.active_weight_sum > 0:
        t = transfer_layer.layer_score
        if t >= 0.60:
            layer_comments.append(
                f"the transfer opportunity assessment is elevated (score {t:.2f}), "
                "suggesting sequence-level features consistent with mobilisation"
            )
        elif t >= 0.30:
            layer_comments.append(
                f"the transfer opportunity layer shows moderate indicators (score {t:.2f})"
            )

    if establishment_layer.active_weight_sum > 0:
        e = establishment_layer.layer_score
        if e >= 0.60:
            layer_comments.append(
                f"establishment compatibility with the assessed host appears plausible "
                f"(score {e:.2f})"
            )
        elif e <= 0.25:
            layer_comments.append(
                f"establishment in the assessed host may be reduced by compositional "
                f"incompatibility (score {e:.2f})"
            )

    if consequence_layer.active_weight_sum > 0:
        c = consequence_layer.layer_score
        if c >= 0.60:
            layer_comments.append(
                f"the functional consequence layer indicates a potentially significant "
                f"payload burden (score {c:.2f})"
            )

    if layer_comments:
        parts.append("Specifically, " + "; ".join(layer_comments) + ".")

    # Top contributors
    if top_contributors:
        if len(top_contributors) == 1:
            parts.append(
                f"The primary driver is {top_contributors[0]}."
            )
        else:
            listed = ", ".join(top_contributors[:-1]) + f", and {top_contributors[-1]}"
            parts.append(
                f"Elevated risk is driven primarily by {listed}."
            )

    # Risk reducers
    if risk_reducers:
        if len(risk_reducers) == 1:
            parts.append(
                f"This is partly offset by low scores for {risk_reducers[0]}."
            )
        else:
            listed = ", ".join(risk_reducers[:-1]) + f", and {risk_reducers[-1]}"
            parts.append(
                f"Risk is partially reduced by low scores for {listed}."
            )

    # Missing data commentary
    if missing_important:
        listed = ", ".join(missing_important)
        parts.append(
            f"Note: {listed} could not be assessed (database(s) not available); "
            "the index may underestimate true risk."
        )

    # Low completeness warning
    if overall_completeness < 0.50:
        parts.append(
            f"Only {overall_completeness:.0%} of features were assessed; "
            "interpret this result with caution."
        )

    # Closing disclaimer
    parts.append(
        "This assessment reflects sequence-level indicators only and does not "
        "constitute a determination that HGT will or will not occur."
    )

    explanation = "  ".join(p for p in parts if p)
    return explanation, top_contributors, risk_reducers, missing_important
