"""QC checks for HDR template arms."""

import re
from collections import Counter
from dataclasses import dataclass, field


# Synthesis length limits (bp).
_SSODN_MAX_LEN = 200   # single-stranded oligonucleotide
_PLASMID_MAX_LEN = 5000  # plasmid / dsDNA donor


@dataclass
class QCResult:
    """Result of a full HDR template quality-control check.

    Attributes
    ----------
    passed : bool
        ``True`` when no hard-fail conditions were detected.  Warnings
        (``warnings`` list non-empty) do not cause a fail on their own.
    warnings : list[str]
        Human-readable warning/error messages, one per issue found.
    checks : list[dict]
        Detailed per-check results from :func:`run_all_checks`.
    """
    passed: bool
    warnings: list[str] = field(default_factory=list)
    checks: list[dict] = field(default_factory=list)


def check_template(template: dict, edit_spec: dict, homology_arm_len: int) -> QCResult:
    """Run comprehensive QC on an HDR donor template.

    Parameters
    ----------
    template : dict
        Template dict as returned by :func:`~src.edit_applicator.build_template`,
        containing keys ``left_arm``, ``edit_seq``, ``right_arm``,
        ``full_template``, ``left_arm_len``, ``right_arm_len``.
    edit_spec : dict
        Edit specification dict (same format as accepted by
        :func:`~src.edit_applicator.apply_edit`), with keys
        ``type``, ``ref_bases``, and ``alt_bases``.
    homology_arm_len : int
        Requested (nominal) arm length used when building the template.

    Returns
    -------
    QCResult
        Dataclass with ``passed``, ``warnings``, and ``checks``.
    """
    warnings: list[str] = []
    hard_fail = False

    left_arm: str = template.get("left_arm", "")
    right_arm: str = template.get("right_arm", "")
    edit_seq: str = template.get("edit_seq", "")
    full_template: str = template.get("full_template", "")

    # ------------------------------------------------------------------
    # 1. Homology arm length
    # ------------------------------------------------------------------
    min_arm = 30
    if len(left_arm) < min_arm:
        msg = (
            f"Left homology arm is only {len(left_arm)} bp "
            f"(minimum required: {min_arm} bp)"
        )
        warnings.append(msg)
        hard_fail = True
    if len(right_arm) < min_arm:
        msg = (
            f"Right homology arm is only {len(right_arm)} bp "
            f"(minimum required: {min_arm} bp)"
        )
        warnings.append(msg)
        hard_fail = True

    # ------------------------------------------------------------------
    # 2. Edit incorporated correctly
    # ------------------------------------------------------------------
    alt_bases = str(edit_spec.get("alt_bases", "")).upper()
    if alt_bases and edit_seq.upper() != alt_bases:
        warnings.append(
            f"Edit sequence in template ('{edit_seq}') does not match "
            f"expected alt_bases ('{alt_bases}')"
        )
        hard_fail = True

    # Verify that the full_template actually contains the edit sequence.
    if alt_bases and alt_bases not in full_template.upper():
        warnings.append(
            f"Alt bases '{alt_bases}' not found anywhere in full template sequence"
        )
        hard_fail = True

    # ------------------------------------------------------------------
    # 3. Internal PAM sites in payload (the edit region only)
    # ------------------------------------------------------------------
    # A PAM site inside the payload means Cas9 could cut the donor.
    payload = edit_seq.upper()
    pam_in_payload: list[int] = []
    for i in range(len(payload) - 2):
        if payload[i + 1] == "G" and payload[i + 2] == "G":
            pam_in_payload.append(i)
    if pam_in_payload:
        warnings.append(
            f"NGG PAM site(s) detected in payload/edit sequence at positions: "
            f"{pam_in_payload} — Cas9 may re-cut after HDR"
        )

    # ------------------------------------------------------------------
    # 4. Template length within synthesis limits
    # ------------------------------------------------------------------
    tlen = len(full_template)
    if tlen <= _SSODN_MAX_LEN:
        synth_class = "ssODN"
    elif tlen <= _PLASMID_MAX_LEN:
        synth_class = "plasmid/dsDNA"
    else:
        synth_class = None
        warnings.append(
            f"Template length {tlen} bp exceeds maximum synthesis limit "
            f"({_PLASMID_MAX_LEN} bp for plasmid donors)"
        )
        hard_fail = True

    if synth_class:
        # Informational — not a warning, but add to checks via a note.
        pass

    if tlen > _SSODN_MAX_LEN and tlen <= _PLASMID_MAX_LEN:
        warnings.append(
            f"Template length {tlen} bp exceeds ssODN limit ({_SSODN_MAX_LEN} bp); "
            "suitable for plasmid/dsDNA synthesis only"
        )

    # ------------------------------------------------------------------
    # 5. Standard sequence-quality checks on both arms
    # ------------------------------------------------------------------
    checks = run_all_checks(left_arm, right_arm)
    for chk in checks:
        if chk["flag"]:
            warnings.append(f"QC check '{chk['check']}': {chk['note']}")

    return QCResult(
        passed=not hard_fail,
        warnings=warnings,
        checks=checks,
    )


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
