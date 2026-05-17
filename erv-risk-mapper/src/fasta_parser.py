"""Parse FASTA genome into chromosome sequences."""

import random
from pathlib import Path


def parse_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    chrom, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if chrom is not None:
                    seqs[chrom] = "".join(parts).upper()
                chrom = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if chrom is not None:
        seqs[chrom] = "".join(parts).upper()
    return seqs


def make_demo_genome(seed: int = 42) -> dict[str, str]:
    """Generate synthetic 100kb genome (5 chroms x 20kb) with 8 planted ERV-like elements."""
    rng = random.Random(seed)
    bases = "ACGT"

    # ERV family signatures (hardcoded k-mers used by erv_detector)
    family_signatures = {
        "HERV-K": ["TGAAAGAC", "GTCTTCA", "AATAAA", "GAGNNNNNCTG", "CCCGGG"],
        "HERV-H": ["TGAAAGAC", "AATAAA", "GCNNNNGC", "TATTAG", "CCCGGG"],
        "HERV-W": ["TGAAAGAC", "GTCTTCA", "AATAAA", "GCCCGG", "TTTAAA"],
        "HERV-E": ["TGAAAGAC", "AATAAA", "CTGCAG", "ACATGT", "CCCGGG"],
        "HERV-L": ["TGAAAGAC", "AATAAA", "GGATCC", "GAGCTC", "GTCTTCA"],
    }

    # Planted ERV elements: (chrom_idx, pos, family, ltr_type)
    planted = [
        (0, 1000, "HERV-K", "full"),
        (0, 12000, "HERV-H", "solo"),
        (1, 500, "HERV-W", "full"),
        (1, 8000, "HERV-K", "partial"),
        (2, 3000, "HERV-E", "full"),
        (2, 15000, "HERV-L", "solo"),
        (3, 6000, "HERV-W", "full"),
        (4, 2000, "HERV-H", "partial"),
    ]

    chroms = {}
    for ci in range(5):
        seq = list(rng.choice(bases) for _ in range(20000))
        chroms[f"chr{ci + 1}"] = seq

    chrom_names = [f"chr{i + 1}" for i in range(5)]

    for chrom_idx, pos, family, ltr_type in planted:
        chrom = chrom_names[chrom_idx]
        seq = chroms[chrom]
        sigs = family_signatures[family]

        # Plant at least 3 signature k-mers near this position
        for j, sig in enumerate(sigs[:4]):
            plant_pos = pos + j * 120
            clean_sig = sig.replace("N", "A").replace("G", "G")
            for k, b in enumerate(clean_sig):
                if plant_pos + k < len(seq):
                    seq[plant_pos + k] = b

        # Plant LTR motifs
        if ltr_type in ("full", "partial"):
            ltr_5 = list("TGAAAGAC")
            for k, b in enumerate(ltr_5):
                if pos + k < len(seq):
                    seq[pos + k] = b

        if ltr_type == "full":
            ltr_3 = list("GTCTTCA")
            end_pos = pos + 4000
            for k, b in enumerate(ltr_3):
                if end_pos + k < len(seq):
                    seq[end_pos + k] = b

        # Plant poly-A signal
        polyA = list("AATAAA")
        pa_pos = pos + 200
        for k, b in enumerate(polyA):
            if pa_pos + k < len(seq):
                seq[pa_pos + k] = b

        # Plant a long ORF (for high-risk elements: full LTR)
        if ltr_type == "full":
            orf_pos = pos + 400
            # Simple ORF: ATG + 510 bases + TAA
            orf = list("ATG") + [rng.choice(bases) for _ in range(510)] + list("TAA")
            for k, b in enumerate(orf):
                if orf_pos + k < len(seq):
                    seq[orf_pos + k] = b

        chroms[chrom] = seq

    return {k: "".join(v) for k, v in chroms.items()}
