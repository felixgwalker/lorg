"""QC checks for HDR template arms."""

import re
from collections import Counter


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def check_gc(seq: str, label: str) -> dict:
    gc = gc_content(seq)
    flag = gc < 0.40 or gc > 0.60
    return {
        "check": f"gc_{label}",
        "value": round(gc, 4),
        "flag": flag,
        "note": f"GC={gc:.1%} {'(LOW)' if gc < 0.40 else '(HIGH)' if gc > 0.60 else '(OK)'}",
    }


def check_homopolymer(seq: str, label: str, max_run: int = 5) -> dict:
    longest = 0
    current = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current += 1
            longest = max(longest, current)
        else:
            current = 1
    flag = longest > max_run
    return {
        "check": f"homopolymer_{label}",
        "value": longest,
        "flag": flag,
        "note": f"Longest homopolymer run: {longest} {'(FLAGGED)' if flag else '(OK)'}",
    }


def check_repetitive(seq: str, label: str, k: int = 4) -> dict:
    if len(seq) < k:
        return {"check": f"repeat_{label}", "value": 0.0, "flag": False, "note": "Sequence too short"}
    kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    counts = Counter(kmers)
    most_common_frac = counts.most_common(1)[0][1] / len(kmers)
    flag = most_common_frac > 0.50
    return {
        "check": f"repeat_{label}",
        "value": round(most_common_frac, 4),
        "flag": flag,
        "note": f"Most common {k}-mer fraction: {most_common_frac:.1%} {'(FLAGGED)' if flag else '(OK)'}",
    }


def _reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def check_hairpin(seq: str, label: str, min_stem: int = 4) -> dict:
    """Count potential hairpin-forming complementary stretches."""
    rc = _reverse_complement(seq)
    n = len(seq)
    hairpin_count = 0
    for i in range(n - min_stem):
        for j in range(i + min_stem + 3, n):
            stem = seq[i:i + min_stem]
            window = rc[n - j - min_stem:n - j]
            if stem == window:
                hairpin_count += 1
    flag = hairpin_count > 3
    return {
        "check": f"hairpin_{label}",
        "value": hairpin_count,
        "flag": flag,
        "note": f"Potential hairpin sites: {hairpin_count} {'(RISK)' if flag else '(OK)'}",
    }


def run_all_checks(left_arm: str, right_arm: str) -> list[dict]:
    checks = []
    checks.append(check_gc(left_arm, "left_arm"))
    checks.append(check_gc(right_arm, "right_arm"))
    checks.append(check_homopolymer(left_arm, "left_arm"))
    checks.append(check_homopolymer(right_arm, "right_arm"))
    checks.append(check_repetitive(left_arm, "left_arm"))
    checks.append(check_repetitive(right_arm, "right_arm"))
    checks.append(check_hairpin(left_arm, "left_arm"))
    checks.append(check_hairpin(right_arm, "right_arm"))
    return checks


def arm_length_variants(ref_seq: str, cut_pos: int, lengths: list[int] = None) -> list[dict]:
    """Test multiple arm lengths and report GC content for each."""
    if lengths is None:
        lengths = [50, 75, 100, 150, 200]
    rows = []
    for arm_len in lengths:
        left_start = max(0, cut_pos - arm_len)
        right_end = min(len(ref_seq), cut_pos + arm_len)
        left_arm = ref_seq[left_start:cut_pos]
        right_arm = ref_seq[cut_pos:right_end]
        rows.append({
            "arm_length": arm_len,
            "left_arm_actual_len": len(left_arm),
            "right_arm_actual_len": len(right_arm),
            "left_arm_gc": round(gc_content(left_arm), 4),
            "right_arm_gc": round(gc_content(right_arm), 4),
            "left_gc_flag": gc_content(left_arm) < 0.40 or gc_content(left_arm) > 0.60,
            "right_gc_flag": gc_content(right_arm) < 0.40 or gc_content(right_arm) > 0.60,
        })
    return rows
