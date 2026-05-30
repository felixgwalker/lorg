"""Configuration for PAM Flexibility Predictor."""

CAS_VARIANT_PAM_TABLE: dict[str, dict] = {
    "SpCas9":    {"pam": "NGG",    "position": "3prime", "length": 3, "size_aa": 1368},
    "SaCas9":    {"pam": "NNGRRT", "position": "3prime", "length": 6, "size_aa": 1053},
    "Cas9-NG":   {"pam": "NG",     "position": "3prime", "length": 2, "size_aa": 1368},
    "SpRY":      {"pam": "NRN",    "position": "3prime", "length": 3, "size_aa": 1368},
    "AsCas12a":  {"pam": "TTTV",   "position": "5prime", "length": 4, "size_aa": 1307},
    "LbCas12a":  {"pam": "TTTN",   "position": "5prime", "length": 4, "size_aa": 1228},
    "Cas12b":    {"pam": "ATTN",   "position": "5prime", "length": 4, "size_aa":  1090},
    "CasX":      {"pam": "TTCN",   "position": "5prime", "length": 4, "size_aa":  986},
}

PAM_DENSITY_WINDOW = 200
