"""Configuration for Repeat Element Classifier."""

REPEAT_MASKER_CLASSES = ["LINE", "SINE", "LTR", "DNA", "Satellite", "Simple_repeat", "Low_complexity"]

MIN_ELEMENT_LENGTH = 50
MAX_DIVERGENCE = 50.0

LANDSCAPE_DIVERGENCE_BINS = list(range(0, 55, 5))
