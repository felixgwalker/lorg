from dataclasses import dataclass
import numpy as np


@dataclass
class FrohResult:
    froh: float
    classification: str   # "low" | "moderate" | "high" | "very high"
    total_roh_bp: int
    genome_size_bp: float


def _classify_froh(froh: float) -> str:
    if froh < 0.05:
        return "low"
    elif froh < 0.125:
        return "moderate"
    elif froh < 0.25:
        return "high"
    else:
        return "very high"


def calculate_froh(
    roh_segments,
    genome_size_bp: float = 2.7e9,
) -> FrohResult:
    """Calculate FROH = total ROH length / genome_size_bp.

    Parameters
    ----------
    roh_segments : iterable of (start, end) int pairs (in bp) representing
                   ROH segments for one individual, or an iterable of lengths.
    genome_size_bp : reference genome size in base-pairs (default 2.7 Gb).

    Returns
    -------
    FrohResult with froh value, classification, total ROH bp, and genome size.
    """
    total_roh_bp = 0
    for seg in roh_segments:
        if hasattr(seg, "__len__") and len(seg) == 2:
            start, end = seg
            total_roh_bp += abs(int(end) - int(start))
        else:
            total_roh_bp += int(seg)

    froh = total_roh_bp / genome_size_bp if genome_size_bp > 0 else 0.0
    froh = float(np.clip(froh, 0.0, 1.0))

    return FrohResult(
        froh=froh,
        classification=_classify_froh(froh),
        total_roh_bp=total_roh_bp,
        genome_size_bp=genome_size_bp,
    )


def detect_roh_sliding_window(
    is_homozygous_row: np.ndarray,
    positions: list,
    window_size: int = 50,
    min_hom_fraction: float = 0.85,
    min_roh_bp: int = 100_000,
) -> list:
    """Detect ROH segments using a sliding window over homozygous sites.

    Parameters
    ----------
    is_homozygous_row : 1-D bool array of length n_sites (True = homozygous).
    positions : list of genomic positions (bp) of length n_sites.
    window_size : number of SNPs per sliding window.
    min_hom_fraction : minimum fraction of homozygous SNPs in window to call ROH.
    min_roh_bp : minimum physical length (bp) for a segment to be reported.

    Returns
    -------
    List of (start_bp, end_bp) tuples for each detected ROH segment.
    """
    n = len(is_homozygous_row)
    if n < window_size:
        return []

    in_roh = np.zeros(n, dtype=bool)
    half = window_size // 2

    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        frac = np.mean(is_homozygous_row[lo:hi])
        if frac >= min_hom_fraction:
            in_roh[i] = True

    # Collect contiguous ROH stretches
    segments = []
    in_seg = False
    seg_start_idx = 0
    for i in range(n):
        if in_roh[i] and not in_seg:
            in_seg = True
            seg_start_idx = i
        elif not in_roh[i] and in_seg:
            in_seg = False
            start_bp = positions[seg_start_idx]
            end_bp = positions[i - 1]
            if abs(end_bp - start_bp) >= min_roh_bp:
                segments.append((start_bp, end_bp))
    if in_seg:
        start_bp = positions[seg_start_idx]
        end_bp = positions[n - 1]
        if abs(end_bp - start_bp) >= min_roh_bp:
            segments.append((start_bp, end_bp))

    return segments


def compute_froh(genotype_matrix, min_run_length=10):
    """Compute FROH from a binary genotype matrix (SNP-count based).

    genotype_matrix: 2D array shape (n_individuals, n_sites), values 0=het, 1=hom
    Returns array of FROH per individual.
    """
    n_indiv, n_sites = genotype_matrix.shape
    froh_values = []
    for i in range(n_indiv):
        gt = genotype_matrix[i]
        total_hom_run = 0
        run = 0
        for j in range(n_sites):
            if gt[j] == 1:
                run += 1
            else:
                if run >= min_run_length:
                    total_hom_run += run
                run = 0
        if run >= min_run_length:
            total_hom_run += run
        froh = total_hom_run / n_sites if n_sites > 0 else 0.0
        froh_values.append(froh)
    return np.array(froh_values)


def froh_from_population_data(population_data, n_snps=500, rng=None):
    """Derive FROH estimates from population data dicts containing F_initial.

    Uses F_initial as approximate FROH (they are related measures of inbreeding).
    """
    if rng is None:
        rng = np.random.default_rng(42)

    for indiv in population_data:
        F = indiv["F_initial"]
        noise = rng.normal(0, 0.02)
        froh = float(np.clip(F + noise, 0.0, 1.0))
        indiv["FROH"] = round(froh, 4)

    return population_data
