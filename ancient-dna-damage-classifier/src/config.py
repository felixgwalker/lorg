"""
Configuration constants for the Ancient DNA Damage Classifier.

All numeric thresholds, defaults, and classification parameters live here.
Update this file to adjust sensitivity without touching algorithm code.
"""

# ── BAM filtering defaults ────────────────────────────────────────────────────
DEFAULT_MIN_MAPQ: int = 30
DEFAULT_MIN_LENGTH: int = 30        # bp; reads shorter than this are excluded
DEFAULT_N_TERMINAL: int = 25        # positions profiled at each terminus

# ── Damage model defaults ─────────────────────────────────────────────────────
BACKGROUND_ERROR_RATE: float = 0.001    # sequencing/mapping error rate for modern reads
DEFAULT_PRIOR_ANCIENT: float = 0.9     # library-level prior P(ancient)

# ── Classification thresholds ─────────────────────────────────────────────────
AUTH_THRESHOLD: float = 0.85    # P(ancient|data) >= AUTH_THRESHOLD -> "authentic"
CONT_THRESHOLD: float = 0.15    # P(ancient|data) <= CONT_THRESHOLD -> "contaminated"
                                 # between these bounds -> "ambiguous"

# ── Decay model quality grading ───────────────────────────────────────────────
# Tuple: (minimum_r_squared, minimum_amplitude, grade_label)
# Evaluated in order; first matching threshold wins.
GRADE_THRESHOLDS: list[tuple[float, float, str]] = [
    (0.90, 0.10, "strong"),
    (0.70, 0.05, "moderate"),
    (0.40, 0.02, "weak"),
    (0.00, 0.00, "absent"),
]

# ── Geometric decay model bounds ─────────────────────────────────────────────
# f(x) = amplitude * (1 - rate)^x + background
DECAY_PARAM_BOUNDS: tuple[list[float], list[float]] = (
    [0.0, 0.0, 0.0],    # lower: amplitude, rate, background
    [1.0, 1.0, 0.5],    # upper: amplitude, rate, background
)
DECAY_PARAM_P0: list[float] = [0.15, 0.3, 0.001]    # initial guess

# ── Output filenames ──────────────────────────────────────────────────────────
OUTFILE_DAMAGE_CSV: str = "damage_frequencies.csv"
OUTFILE_PLOT_PNG: str = "damage_profile.png"
OUTFILE_PLOT_SVG: str = "damage_profile.svg"
OUTFILE_READS_TSV: str = "read_classifications.tsv"
OUTFILE_SUMMARY_JSON: str = "summary_report.json"
OUTFILE_SUMMARY_TXT: str = "summary_report.txt"

# ── Minimum fraction of reads with usable MD tags ─────────────────────────────
MIN_MD_TAG_FRACTION: float = 0.10    # error if fewer than this fraction have MD tags
