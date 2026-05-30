"""Sliding window ROH detection from per-individual genotype arrays."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src.vcf_parser import GenotypeData


@dataclass
class ROHSegment:
    """A single Run of Homozygosity segment."""
    individual_id: str
    chrom: str
    start_pos: int
    end_pos: int
    n_snps: int
    length_bp: int
    length_class: str   # "short" | "medium" | "long"
    mean_homozygosity: float


# Length class boundaries (base pairs)
# short  : <  100 kb   — ancient signal (hundreds of generations ago)
# medium : 100 kb – 1 Mb — moderate inbreeding (tens of generations)
# long   : >  1 Mb     — recent inbreeding (last few generations)
SHORT_MAX_BP = 100_000       # upper boundary for "short"
MEDIUM_MAX_BP = 1_000_000    # upper boundary for "medium"


def detect_roh(
    gdata: GenotypeData,
    window_size: int = 50,
    step: int = 1,
    homo_threshold: float = 0.95,
    min_snps: int = 10,
) -> list[ROHSegment]:
    """Detect ROH segments using a sliding window approach.

    The algorithm marks each SNP position as "in an ROH" if the window
    centred on it (or starting at it) has a homozygosity fraction >=
    *homo_threshold*.  Contiguous runs of marked positions are then merged
    into ROHSegment objects.  Segments with fewer than *min_snps* SNPs are
    discarded.

    Length classes assigned to each segment:
      ``"short"``  — < 100 kb  (ancient bottleneck signal)
      ``"medium"`` — 100 kb – 1 Mb  (moderate / historical inbreeding)
      ``"long"``   — > 1 Mb  (recent inbreeding, last few generations)
    """
    positions = np.array(gdata.positions, dtype=np.int64)
    homo = np.array(gdata.is_homozygous, dtype=np.float64)
    n = len(positions)

    if n < window_size:
        return []

    in_roh = np.zeros(n, dtype=bool)
    for start in range(0, n - window_size + 1, step):
        end = start + window_size
        frac = homo[start:end].mean()
        if frac >= homo_threshold:
            in_roh[start:end] = True

    # Merge contiguous flagged positions into segments
    segments: list[ROHSegment] = []
    i = 0
    while i < n:
        if in_roh[i]:
            j = i
            while j < n and in_roh[j]:
                j += 1
            seg_snps = j - i
            if seg_snps >= min_snps:
                start_pos = int(positions[i])
                end_pos = int(positions[j - 1])
                length_bp = end_pos - start_pos
                mean_homo = float(homo[i:j].mean())
                lclass = _length_class(length_bp)
                segments.append(ROHSegment(
                    individual_id=gdata.individual_id,
                    chrom=gdata.chrom,
                    start_pos=start_pos,
                    end_pos=end_pos,
                    n_snps=seg_snps,
                    length_bp=length_bp,
                    length_class=lclass,
                    mean_homozygosity=round(mean_homo, 4),
                ))
            i = j
        else:
            i += 1
    return segments


def _length_class(length_bp: int) -> str:
    """Classify an ROH by physical length into short / medium / long."""
    if length_bp >= MEDIUM_MAX_BP:
        return "long"
    if length_bp >= SHORT_MAX_BP:
        return "medium"
    return "short"
