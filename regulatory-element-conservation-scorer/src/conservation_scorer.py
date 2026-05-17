from .pwm_scanner import scan_pwms


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
