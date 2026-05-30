"""Configuration for Prime Edit Design Assistant."""

PAM_PATTERNS: dict[str, str] = {
    "SpCas9": "NGG",
    "SaCas9": "NNGRRT",
    "Cas9-NG": "NG",
    "SpRY": "NRN",
    "Cas12a": "TTTV",
}

PBS_MIN = 8
PBS_MAX = 15
PBS_OPTIMAL_MIN = 10
PBS_OPTIMAL_MAX = 13
PBS_GC_OPTIMAL_MIN = 0.30
PBS_GC_OPTIMAL_MAX = 0.70

RT_MIN = 10
RT_MAX = 16
RT_OPTIMAL_MIN = 10
RT_OPTIMAL_MAX = 14

NICK_WINDOW_MIN = 40
NICK_WINDOW_MAX = 90

SPACER_LENGTH = 20
