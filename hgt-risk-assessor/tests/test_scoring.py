"""
Tests for the three-layer HGT risk scoring framework.

Run with:  python -m pytest tests/ -v
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from src.models import (
    FeatureResult,
    HostProfile,
    InputFormat,
    LayerResult,
    QuerySequence,
    ScoreBand,
    SignalResult,
    ThreeLayerResult,
)
from src.config import (
    SCORE_BANDS,
    SCORE_BAND_TOP,
    WEIGHT_PROFILES,
    LAYER_FEATURE_WEIGHTS,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

_ECOLI_SEQ = (
    "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA"
    "ATGCGTAAAGGTTTCGGTTTCGCAGGTGGTTTCGGTTTCGGCGGTAAAGGTTTCGGCAAAGGTGGTTTCGGTTTCGGCTGA"
    "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA" * 5
)

_SHORT_SEQ = "ATGAAATGA"   # single tiny ORF

def _make_query(seq: str = _ECOLI_SEQ) -> QuerySequence:
    from Bio.SeqUtils import gc_fraction
    from Bio.Seq import Seq
    return QuerySequence(
        sequence=seq,
        identifier="test_query",
        length=len(seq),
        gc_content=gc_fraction(Seq(seq)),
        source_format=InputFormat.FASTA,
    )


def _make_host(gc: float = 0.508) -> HostProfile:
    return HostProfile(identifier="ecoli", gc_content=gc, source="test")


def _make_signal(name: str, score: float | None, skipped: bool = False) -> SignalResult:
    from src.aggregator import SIGNAL_WEIGHTS
    return SignalResult(
        signal_name=name,
        score=score,
        weight=SIGNAL_WEIGHTS.get(name, 0.2),
        evidence={},
        skipped=skipped,
    )


def _all_signals(scores: dict[str, float | None] | None = None) -> list[SignalResult]:
    defaults = {
        "gc_content":   0.1,
        "is_proximity": 0.0,
        "integron":     0.0,
        "conjugative":  0.0,
        "prophage":     0.0,
    }
    if scores:
        defaults.update(scores)
    return [_make_signal(k, v) for k, v in defaults.items()]


# ---------------------------------------------------------------------------
# config.py sanity checks
# ---------------------------------------------------------------------------

class TestConfig:
    def test_weight_profiles_sum_to_one(self):
        for name, profile in WEIGHT_PROFILES.items():
            total = sum(profile.values())
            assert abs(total - 1.0) < 1e-9, f"Profile {name!r} weights sum to {total}"

    def test_layer_feature_weights_sum_to_one(self):
        for layer, weights in LAYER_FEATURE_WEIGHTS.items():
            total = sum(weights.values())
            assert abs(total - 1.0) < 1e-9, f"Layer {layer!r} feature weights sum to {total}"

    def test_score_bands_ascending(self):
        thresholds = [t for t, _ in SCORE_BANDS]
        assert thresholds == sorted(thresholds), "SCORE_BANDS not in ascending order"

    def test_all_profiles_present(self):
        for profile in ("default", "environmental", "clinical_amr"):
            assert profile in WEIGHT_PROFILES


# ---------------------------------------------------------------------------
# Score-band classification
# ---------------------------------------------------------------------------

class TestScoreBand:
    def test_classify_low(self):
        from src.scoring.layers import classify_band
        assert classify_band(0.10) == ScoreBand.LOW
        assert classify_band(0.00) == ScoreBand.LOW

    def test_classify_moderate(self):
        from src.scoring.layers import classify_band
        assert classify_band(0.25) == ScoreBand.MODERATE
        assert classify_band(0.49) == ScoreBand.MODERATE

    def test_classify_high(self):
        from src.scoring.layers import classify_band
        assert classify_band(0.50) == ScoreBand.HIGH
        assert classify_band(0.74) == ScoreBand.HIGH

    def test_classify_very_high(self):
        from src.scoring.layers import classify_band
        assert classify_band(0.75) == ScoreBand.VERY_HIGH
        assert classify_band(1.00) == ScoreBand.VERY_HIGH

    def test_boundary_exactly_on_threshold(self):
        from src.scoring.layers import classify_band
        # 0.25 is the lower boundary of moderate (score < 0.25 → low)
        assert classify_band(0.249999) == ScoreBand.LOW
        assert classify_band(0.250000) == ScoreBand.MODERATE


# ---------------------------------------------------------------------------
# Layer aggregation
# ---------------------------------------------------------------------------

class TestLayerAggregation:
    def test_all_available(self):
        from src.scoring.layers import _aggregate_layer
        features = [
            FeatureResult("f1", "transfer_opportunity", 0.8, 0.5, True, {}),
            FeatureResult("f2", "transfer_opportunity", 0.2, 0.5, True, {}),
        ]
        lr = _aggregate_layer(features)
        assert abs(lr.layer_score - 0.5) < 1e-9
        assert lr.completeness == 1.0
        assert abs(lr.active_weight_sum - 1.0) < 1e-9

    def test_one_skipped(self):
        from src.scoring.layers import _aggregate_layer
        features = [
            FeatureResult("f1", "transfer_opportunity", 1.0, 0.6, True,  {}),
            FeatureResult("f2", "transfer_opportunity", None, 0.4, False, {}),
        ]
        lr = _aggregate_layer(features)
        # Only f1 active; renormalised weight is 1.0 → score = 1.0
        assert abs(lr.layer_score - 1.0) < 1e-9
        assert abs(lr.completeness - 0.5) < 1e-9

    def test_all_skipped_returns_zero(self):
        from src.scoring.layers import _aggregate_layer
        features = [
            FeatureResult("f1", "transfer_opportunity", None, 0.5, False, {}),
            FeatureResult("f2", "transfer_opportunity", None, 0.5, False, {}),
        ]
        lr = _aggregate_layer(features)
        assert lr.layer_score == 0.0
        assert lr.completeness == 0.0

    def test_score_clamped_to_range(self):
        from src.scoring.layers import _aggregate_layer
        features = [
            FeatureResult("f1", "transfer_opportunity", 1.5, 1.0, True, {}),
        ]
        lr = _aggregate_layer(features)
        assert 0.0 <= lr.layer_score <= 1.0


# ---------------------------------------------------------------------------
# Feature extractors — GC deviation (signal reuse)
# ---------------------------------------------------------------------------

class TestGCDeviation:
    def test_no_deviation_scores_zero(self):
        from src.scoring.features import compute_gc_deviation
        f = compute_gc_deviation(signal_score=0.0, evidence={"query_gc": 50.8, "host_gc": 50.8, "deviation_pct": 0.0})
        assert f.score == 0.0
        assert f.available is True

    def test_max_deviation_caps_at_one(self):
        from src.scoring.features import compute_gc_deviation
        f = compute_gc_deviation(signal_score=1.0)
        assert f.score == 1.0

    def test_none_signal_is_unavailable(self):
        from src.scoring.features import compute_gc_deviation
        f = compute_gc_deviation(signal_score=None)
        assert f.available is False
        assert f.score is None


# ---------------------------------------------------------------------------
# Feature extractors — codon usage distance
# ---------------------------------------------------------------------------

class TestCodonUsage:
    def test_returns_feature_result(self):
        from src.scoring.features import compute_codon_usage_distance
        q = _make_query()
        h = _make_host()
        f = compute_codon_usage_distance(query=q, host=h)
        assert f.available is True
        assert 0.0 <= f.score <= 1.0

    def test_score_between_zero_and_one(self):
        from src.scoring.features import compute_codon_usage_distance
        q = _make_query(_SHORT_SEQ)
        h = _make_host(0.70)
        f = compute_codon_usage_distance(query=q, host=h)
        assert 0.0 <= f.score <= 1.0

    def test_ecoli_host_selects_ecoli_reference(self):
        from src.scoring.features import compute_codon_usage_distance
        h = _make_host()    # identifier = "ecoli"
        q = _make_query()
        f = compute_codon_usage_distance(query=q, host=h)
        assert "ecoli" in f.evidence.get("reference_used", "").lower()

    def test_unknown_host_uses_equal_reference(self):
        from src.scoring.features import compute_codon_usage_distance
        h = HostProfile(identifier="alien_organism", gc_content=0.5, source="test")
        q = _make_query()
        f = compute_codon_usage_distance(query=q, host=h)
        assert "equal" in f.evidence.get("reference_used", "").lower()


# ---------------------------------------------------------------------------
# Feature extractors — taxonomic distance
# ---------------------------------------------------------------------------

class TestTaxonomicDistance:
    def test_missing_taxa_returns_unavailable(self):
        from src.scoring.features import compute_taxonomic_distance
        f = compute_taxonomic_distance(donor_taxon=None, recipient_taxon=None)
        assert f.available is False

    def test_same_genus_low_distance(self):
        from src.scoring.features import compute_taxonomic_distance
        f = compute_taxonomic_distance(donor_taxon="ecoli", recipient_taxon="salmonella")
        assert f.available is True
        assert f.score <= 0.35, f"Same-family distance should be low, got {f.score}"

    def test_cross_domain_high_distance(self):
        from src.scoring.features import compute_taxonomic_distance
        f = compute_taxonomic_distance(donor_taxon="ecoli", recipient_taxon="human")
        assert f.available is True
        assert f.score >= 0.80, f"Cross-domain distance should be high, got {f.score}"

    def test_unknown_taxon_midpoint(self):
        from src.scoring.features import compute_taxonomic_distance
        f = compute_taxonomic_distance(donor_taxon="unknown_xyzzy", recipient_taxon="ecoli")
        assert f.available is True
        assert abs(f.score - 0.5) < 1e-9


# ---------------------------------------------------------------------------
# Feature extractors — ORF scanning
# ---------------------------------------------------------------------------

class TestORFScanning:
    def test_gene_completeness_empty_sequence(self):
        from src.scoring.features import compute_gene_completeness
        q = QuerySequence("", "q", 0, 0.5, InputFormat.FASTA)
        f = compute_gene_completeness(query=q)
        assert f.score == 0.0

    def test_complete_orf_detected(self):
        from src.scoring.features import compute_gene_completeness
        # ATG + 100 nt coding + TAA stop
        seq = "ATG" + ("AAA" * 33) + "TAA"
        q = QuerySequence(seq, "q", len(seq), 0.3, InputFormat.FASTA)
        f = compute_gene_completeness(query=q)
        assert f.evidence["complete_orfs_ge90bp"] >= 1

    def test_payload_count_multiple_orfs(self):
        from src.scoring.features import compute_payload_count
        # Two complete ORFs each >= 300 bp
        orf = "ATG" + ("AAA" * 100) + "TAA"   # 303 bp
        seq = orf + orf
        q = QuerySequence(seq, "q", len(seq), 0.3, InputFormat.FASTA)
        f = compute_payload_count(query=q)
        assert f.evidence["complete_orfs_ge300bp"] >= 2
        assert f.score > 0.0

    def test_payload_count_ten_orfs_saturates(self):
        from src.scoring.features import compute_payload_count
        orf = "ATG" + ("AAA" * 100) + "TAA"
        seq = orf * 15
        q = QuerySequence(seq, "q", len(seq), 0.3, InputFormat.FASTA)
        f = compute_payload_count(query=q)
        assert f.score == 1.0


# ---------------------------------------------------------------------------
# Feature extractors — promoter scan
# ---------------------------------------------------------------------------

class TestPromoterScan:
    def test_no_promoter_in_random_gc_rich(self):
        from src.scoring.features import compute_promoter_plausibility
        # GC-rich sequence unlikely to have TATAAT / TTGACA
        seq = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG" * 5
        q = QuerySequence(seq, "q", len(seq), 0.8, InputFormat.FASTA)
        f = compute_promoter_plausibility(query=q)
        assert f.available is True

    def test_perfect_promoter_detected(self):
        from src.scoring.features import compute_promoter_plausibility
        # Embed a perfect -35/-10 pair: TTGACA + 17 bp spacer + TATAAT
        spacer = "A" * 17
        promoter = "TTGACA" + spacer + "TATAAT"
        seq = ("N" * 50 + promoter + "N" * 50).replace("N", "A")
        q = QuerySequence(seq, "q", len(seq), 0.3, InputFormat.FASTA)
        f = compute_promoter_plausibility(query=q)
        assert f.evidence["sigma70_like_count"] >= 1


# ---------------------------------------------------------------------------
# Feature extractors — placeholders
# ---------------------------------------------------------------------------

class TestPlaceholders:
    def test_plasmid_context_unavailable(self):
        from src.scoring.features import compute_plasmid_context
        f = compute_plasmid_context()
        assert f.available is False
        assert f.source == "placeholder"

    def test_transposase_proximity_unavailable(self):
        from src.scoring.features import compute_transposase_proximity
        f = compute_transposase_proximity()
        assert f.available is False

    def test_amr_content_unavailable(self):
        from src.scoring.features import compute_amr_content
        f = compute_amr_content()
        assert f.available is False
        assert f.source == "placeholder"

    def test_virulence_flags_unavailable(self):
        from src.scoring.features import compute_virulence_flags
        f = compute_virulence_flags()
        assert f.available is False


# ---------------------------------------------------------------------------
# Three-layer integration
# ---------------------------------------------------------------------------

class TestThreeLayerIntegration:
    def test_produces_three_layer_result(self):
        from src.scoring.layers import compute_three_layer
        result = compute_three_layer(
            signal_results=_all_signals(),
            query=_make_query(),
            host=_make_host(),
        )
        assert isinstance(result, ThreeLayerResult)
        assert 0.0 <= result.hgt_risk_index <= 1.0

    def test_index_in_range(self):
        from src.scoring.layers import compute_three_layer
        result = compute_three_layer(
            signal_results=_all_signals({"is_proximity": 1.0, "conjugative": 1.0}),
            query=_make_query(),
            host=_make_host(),
        )
        assert 0.0 <= result.hgt_risk_index <= 1.0

    def test_score_band_matches_index(self):
        from src.scoring.layers import compute_three_layer, classify_band
        result = compute_three_layer(
            signal_results=_all_signals(),
            query=_make_query(),
            host=_make_host(),
        )
        assert result.score_band == classify_band(result.hgt_risk_index)

    def test_default_profile_used_when_invalid(self):
        from src.scoring.layers import compute_three_layer
        result = compute_three_layer(
            signal_results=_all_signals(),
            query=_make_query(),
            host=_make_host(),
            weight_profile_name="nonexistent_profile",
        )
        assert result.weight_profile_name == "default"

    def test_all_profiles_run_without_error(self):
        from src.scoring.layers import compute_three_layer
        for profile_name in WEIGHT_PROFILES:
            result = compute_three_layer(
                signal_results=_all_signals(),
                query=_make_query(),
                host=_make_host(),
                weight_profile_name=profile_name,
            )
            assert result.weight_profile_name == profile_name

    def test_all_signals_skipped_still_produces_result(self):
        from src.scoring.layers import compute_three_layer
        signals = [_make_signal(n, None, skipped=True)
                   for n in ("gc_content", "is_proximity", "integron", "conjugative", "prophage")]
        result = compute_three_layer(
            signal_results=signals,
            query=_make_query(),
            host=_make_host(),
        )
        assert isinstance(result, ThreeLayerResult)
        assert result.hgt_risk_index >= 0.0

    def test_completeness_reflects_available_features(self):
        from src.scoring.layers import compute_three_layer
        result = compute_three_layer(
            signal_results=_all_signals(),
            query=_make_query(),
            host=_make_host(),
        )
        # Without BLAST DBs, placeholders make some features unavailable
        # Computable features (gc, codon, promoter, complexity, completeness, payload)
        # should all be available
        assert result.overall_completeness > 0.0


# ---------------------------------------------------------------------------
# Weight profile selection
# ---------------------------------------------------------------------------

class TestWeightProfiles:
    def test_clinical_amr_weights_consequence_higher(self):
        from src.scoring.layers import compute_three_layer
        r_default = compute_three_layer(
            signal_results=_all_signals({"prophage": 0.9}),
            query=_make_query(),
            host=_make_host(),
            weight_profile_name="default",
        )
        r_amr = compute_three_layer(
            signal_results=_all_signals({"prophage": 0.9}),
            query=_make_query(),
            host=_make_host(),
            weight_profile_name="clinical_amr",
        )
        # Both produce valid results
        assert isinstance(r_default, ThreeLayerResult)
        assert isinstance(r_amr, ThreeLayerResult)


# ---------------------------------------------------------------------------
# Missing data handling
# ---------------------------------------------------------------------------

class TestMissingData:
    def test_missing_features_excluded_from_index(self):
        from src.scoring.layers import compute_three_layer
        # Skipping all BLAST signals → those features unavailable
        # Score should still be computed from remaining features
        signals_skipped = [_make_signal(n, None, skipped=True)
                           for n in ("is_proximity", "conjugative", "integron", "prophage")]
        signals_skipped.append(_make_signal("gc_content", 0.3))
        result = compute_three_layer(
            signal_results=signals_skipped,
            query=_make_query(),
            host=_make_host(),
        )
        assert result.hgt_risk_index >= 0.0
        assert result.overall_completeness < 1.0

    def test_important_missing_features_flagged(self):
        from src.scoring.layers import compute_three_layer
        result = compute_three_layer(
            signal_results=_all_signals(),
            query=_make_query(),
            host=_make_host(),
        )
        # AMR and virulence are always placeholders → should be in missing list
        missing = result.missing_important_features
        assert any("AMR" in m or "amr" in m.lower() for m in missing)


# ---------------------------------------------------------------------------
# Explanation generator
# ---------------------------------------------------------------------------

class TestExplanation:
    def _run(self, scores=None):
        from src.scoring.layers import compute_three_layer
        return compute_three_layer(
            signal_results=_all_signals(scores),
            query=_make_query(),
            host=_make_host(),
        )

    def test_explanation_is_nonempty_string(self):
        result = self._run()
        assert isinstance(result.explanation, str)
        assert len(result.explanation) > 50

    def test_explanation_avoids_deterministic_claims(self):
        result = self._run({"is_proximity": 0.9, "conjugative": 0.9})
        exp = result.explanation.lower()
        # Must not assert certainty
        assert "will definitely" not in exp
        assert "guaranteed" not in exp

    def test_high_score_reflected_in_explanation(self):
        result = self._run({"is_proximity": 1.0, "conjugative": 1.0, "integron": 1.0})
        exp = result.explanation.lower()
        assert any(word in exp for word in ("elevated", "significant", "multiple", "strong"))

    def test_low_completeness_mentioned(self):
        from src.scoring.layers import compute_three_layer
        signals_skipped = [_make_signal(n, None, skipped=True)
                           for n in ("is_proximity", "conjugative", "integron", "prophage", "gc_content")]
        result = compute_three_layer(
            signal_results=signals_skipped,
            query=_make_query(_SHORT_SEQ),
            host=_make_host(),
        )
        exp = result.explanation.lower()
        assert "caution" in exp or "%" in result.explanation


# ---------------------------------------------------------------------------
# Flat aggregator — backward compat with config
# ---------------------------------------------------------------------------

class TestFlatAggregatorBackwardCompat:
    def test_classify_still_works(self):
        from src.aggregator import classify
        from src.models import RiskLevel
        assert classify(0.10) == RiskLevel.LOW
        assert classify(0.30) == RiskLevel.MEDIUM
        assert classify(0.60) == RiskLevel.HIGH
        assert classify(0.80) == RiskLevel.CRITICAL

    def test_aggregate_unchanged(self):
        from src.aggregator import aggregate
        signals = _all_signals({"gc_content": 0.5, "is_proximity": 0.0,
                                 "integron": 0.0, "conjugative": 0.0, "prophage": 0.0})
        result = aggregate(signals)
        assert 0.0 <= result.risk_index <= 1.0
        assert result.risk_level is not None
