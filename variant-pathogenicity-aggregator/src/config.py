"""Configuration for Variant Pathogenicity Aggregator."""

ACMG_CLASSES = ["pathogenic", "likely_pathogenic", "VUS", "likely_benign", "benign"]

PATHOGENIC_THRESHOLD = 6
LIKELY_PATHOGENIC_THRESHOLD = 4
BENIGN_THRESHOLD = -6
LIKELY_BENIGN_THRESHOLD = -4

CRITERION_WEIGHTS: dict[str, int] = {
    "PVS1": 8, "PS1": 4, "PS2": 4, "PS3": 4, "PS4": 4,
    "PM1": 2, "PM2": 2, "PM3": 2, "PM4": 2, "PM5": 2, "PM6": 2,
    "PP1": 1, "PP2": 1, "PP3": 1, "PP4": 1, "PP5": 1,
    "BA1": -8,
    "BS1": -4, "BS2": -4, "BS3": -4, "BS4": -4,
    "BP1": -1, "BP2": -1, "BP3": -1, "BP4": -1,
    "BP5": -1, "BP6": -1, "BP7": -1,
}

MIN_CLINVAR_STARS = 1
