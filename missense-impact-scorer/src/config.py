"""Configuration for Missense Impact Scorer."""

BENIGN_THRESHOLD = 0.3
PATHOGENIC_THRESHOLD = 0.7

CONSERVATION_WEIGHT = 0.40
BLOSUM_WEIGHT = 0.30
PHYSICOCHEMICAL_WEIGHT = 0.30

PHYLOP_CUTOFF = 2.5
GERP_CUTOFF = 2.0

BLOSUM62_MATRIX_URL = "https://www.ncbi.nlm.nih.gov/Class/StructMeth/Alignment/blosum62.txt"

PHYSICOCHEMICAL_GROUPS: dict[str, str] = {
    "A": "nonpolar", "V": "nonpolar", "I": "nonpolar", "L": "nonpolar",
    "M": "nonpolar", "F": "aromatic", "W": "aromatic", "P": "nonpolar",
    "G": "nonpolar", "S": "polar", "T": "polar", "C": "polar",
    "Y": "aromatic", "H": "positive", "D": "negative", "E": "negative",
    "N": "polar", "Q": "polar", "K": "positive", "R": "positive",
}
