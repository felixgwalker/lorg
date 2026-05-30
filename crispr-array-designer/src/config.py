"""Configuration for CRISPR Array Designer."""

CAS_SYSTEM_PAMS: dict[str, str] = {
    "Cas12a": "TTTV",
    "Cas9": "NGG",
    "Cas12b": "ATTN",
}

DIRECT_REPEATS: dict[str, str] = {
    "Cas12a": "AATTTCTACTGTTGTAGAT",
    "Cas9": "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
}

SPACER_LENGTH: dict[str, int] = {
    "Cas12a": 23,
    "Cas9": 20,
    "Cas12b": 20,
}

MIN_SPACER_GC = 0.25
MAX_SPACER_GC = 0.75
MAX_ARRAY_SPACERS = 8
