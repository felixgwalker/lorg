"""Configuration for Transcription Factor Site Scanner."""

DEFAULT_P_VALUE_THRESHOLD = 1e-4
MIN_IC_CONTENT = 8.0
PSEUDOCOUNT = 0.1

BACKGROUND_FREQ = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

SCAN_BOTH_STRANDS = True
FDR_CORRECTION = True
