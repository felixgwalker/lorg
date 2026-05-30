import math
from dataclasses import dataclass
from typing import List

from .pwm_scanner import scan_pwms


@dataclass
class ConservationResult:
    """Results of multi-species conservation scoring for a single element."""
    mean_pairwise_identity: float        # Mean % identity across all species pairs (0–100)
    conserved_positions_fraction: float  # Fraction of positions with ≥90 % identity across all seqs
    entropy_conservation_score: float    # Mean (1 - H/H_max) across alignment columns (0–1)
    n_sequences: int                     # Number of sequences included in the analysis
    n_positions: int                     # Alignment length used


def score_conservation(seqs_across_species: List[str]) -> ConservationResult:
    """Compute conservation metrics for a set of aligned sequences from N species.

    All sequences should represent the same regulatory element extracted from
    orthologous loci; they may differ in length (shorter seqs are ignored at
    positions beyond their end).

    Parameters
    ----------
    seqs_across_species : list of str
        Aligned DNA sequences, one per species.  Non-ACGT characters and gaps
        are treated as unknowns and excluded from identity counts.

    Returns
    -------
    ConservationResult dataclass with three complementary metrics:

    mean_pairwise_identity
        Average % identical bases across all N*(N-1)/2 species pairs, computed
        over the overlapping region of each pair.  Range 0–100.

    conserved_positions_fraction
        Fraction of alignment columns in which ≥90 % of sequences share the
        same base (strict majority criterion).  Range 0–1.

    entropy_conservation_score
        For each column the Shannon entropy H = -sum(p * log2(p)) over the
        four nucleotide frequencies.  H_max = log2(4) = 2.  Conservation at
        that column = 1 - H/H_max.  The reported value is the mean across all
        columns that have at least one resolved base.  Range 0–1.
    """
    seqs = [s.upper() for s in seqs_across_species if s]
    n = len(seqs)

    if n == 0:
        return ConservationResult(
            mean_pairwise_identity=0.0,
            conserved_positions_fraction=0.0,
            entropy_conservation_score=0.0,
            n_sequences=0,
            n_positions=0,
        )

    if n == 1:
        return ConservationResult(
            mean_pairwise_identity=100.0,
            conserved_positions_fraction=1.0,
            entropy_conservation_score=1.0,
            n_sequences=1,
            n_positions=len(seqs[0]),
        )

    max_len = max(len(s) for s in seqs)
    H_max = math.log2(4)  # 2.0

    # ------------------------------------------------------------------ #
    # 1. Mean pairwise identity
    # ------------------------------------------------------------------ #
    pair_identities = []
    for i in range(n):
        for j in range(i + 1, n):
            s1, s2 = seqs[i], seqs[j]
            overlap = min(len(s1), len(s2))
            if overlap == 0:
                pair_identities.append(0.0)
                continue
            valid = 0
            matches = 0
            for k in range(overlap):
                b1, b2 = s1[k], s2[k]
                if b1 in "ACGT" and b2 in "ACGT":
                    valid += 1
                    if b1 == b2:
                        matches += 1
            if valid > 0:
                pair_identities.append(100.0 * matches / valid)
            else:
                pair_identities.append(0.0)

    mean_pairwise_identity = sum(pair_identities) / len(pair_identities)

    # ------------------------------------------------------------------ #
    # 2. Per-column analysis: conserved_positions_fraction & entropy score
    # ------------------------------------------------------------------ #
    conserved_col_count = 0
    entropy_scores = []

    for pos in range(max_len):
        bases_at_pos = []
        for s in seqs:
            if pos < len(s) and s[pos] in "ACGT":
                bases_at_pos.append(s[pos])

        if not bases_at_pos:
            continue

        total = len(bases_at_pos)

        # Conserved-positions metric: ≥90 % of present sequences agree
        counts = {"A": 0, "C": 0, "G": 0, "T": 0}
        for b in bases_at_pos:
            counts[b] += 1
        max_count = max(counts.values())
        if max_count / total >= 0.90:
            conserved_col_count += 1

        # Shannon entropy conservation
        H = 0.0
        for cnt in counts.values():
            if cnt > 0:
                p = cnt / total
                H -= p * math.log2(p)
        if H_max > 0:
            col_conservation = 1.0 - H / H_max
        else:
            col_conservation = 1.0
        entropy_scores.append(col_conservation)

    n_valid_positions = len(entropy_scores)
    conserved_positions_fraction = (
        conserved_col_count / n_valid_positions if n_valid_positions > 0 else 0.0
    )
    entropy_conservation_score = (
        sum(entropy_scores) / n_valid_positions if n_valid_positions > 0 else 0.0
    )

    return ConservationResult(
        mean_pairwise_identity=round(mean_pairwise_identity, 4),
        conserved_positions_fraction=round(conserved_positions_fraction, 4),
        entropy_conservation_score=round(entropy_conservation_score, 4),
        n_sequences=n,
        n_positions=n_valid_positions,
    )


def pairwise_identity(seq1, seq2):
    if not seq1 or not seq2:
        return 0.0
    length = min(len(seq1), len(seq2))
    matches = sum(1 for a, b in zip(seq1[:length], seq2[:length]) if a == b)
    return matches / length


def score_element_conservation(element_id, species_seqs, species_list, ref_species):
    ref_seq = species_seqs.get(ref_species, {}).get(element_id, "")
    ref_pwm = scan_pwms(ref_seq)

    per_species_identity = {}
    per_species_motifs = {}
    all_motif_names = list(ref_pwm.keys())

    for sp in species_list:
        seq = species_seqs.get(sp, {}).get(element_id, "")
        identity = pairwise_identity(ref_seq, seq)
        per_species_identity[sp] = round(identity, 4)

        pwm_result = scan_pwms(seq)
        motif_found = {m: pwm_result[m]["found"] for m in all_motif_names}
        per_species_motifs[sp] = motif_found

    ref_motifs_found = {m: ref_pwm[m]["found"] for m in all_motif_names}
    motifs_in_all = {
        m: all(per_species_motifs[sp].get(m, False) for sp in species_list)
        for m in all_motif_names
    }
    motif_retention = sum(1 for found in motifs_in_all.values() if found) / max(len(all_motif_names), 1)

    non_ref = [sp for sp in species_list if sp != ref_species]
    mean_identity = (
        sum(per_species_identity[sp] for sp in non_ref) / len(non_ref)
        if non_ref else 1.0
    )

    combined_score = 0.6 * mean_identity + 0.4 * motif_retention

    if combined_score > 0.7:
        classification = "conserved"
    elif combined_score >= 0.4:
        classification = "partially_conserved"
    else:
        classification = "diverged"

    return {
        "element_id": element_id,
        "mean_sequence_identity": round(mean_identity, 4),
        "motif_retention": round(motif_retention, 4),
        "combined_score": round(combined_score, 4),
        "classification": classification,
        "per_species_identity": per_species_identity,
        "per_species_motifs": per_species_motifs,
        "ref_motifs": ref_motifs_found,
    }
