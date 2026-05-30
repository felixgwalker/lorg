"""Configuration for Cas Variant Selector."""

CAS_VARIANT_DATABASE: dict[str, dict] = {
    "SpCas9":    {"pam": "NGG",    "size_aa": 1368, "editing_types": ["knockout", "base-edit", "prime-edit", "activation"], "aav_compatible": False},
    "SaCas9":    {"pam": "NNGRRT","size_aa": 1053, "editing_types": ["knockout", "base-edit"], "aav_compatible": True},
    "Cas9-NG":   {"pam": "NG",    "size_aa": 1368, "editing_types": ["knockout", "base-edit", "prime-edit"], "aav_compatible": False},
    "SpRY":      {"pam": "NRN",   "size_aa": 1368, "editing_types": ["knockout", "base-edit", "prime-edit"], "aav_compatible": False},
    "AsCas12a":  {"pam": "TTTV",  "size_aa": 1307, "editing_types": ["knockout"], "aav_compatible": False},
    "LbCas12a":  {"pam": "TTTN",  "size_aa": 1228, "editing_types": ["knockout"], "aav_compatible": True},
    "Cas12b":    {"pam": "ATTN",  "size_aa": 1090, "editing_types": ["knockout"], "aav_compatible": True},
    "CasX":      {"pam": "TTCN",  "size_aa":  986, "editing_types": ["knockout"], "aav_compatible": True},
    "ABE8e":     {"pam": "NGG",   "size_aa": 1450, "editing_types": ["base-edit"], "aav_compatible": False},
    "CBE4max":   {"pam": "NGG",   "size_aa": 1400, "editing_types": ["base-edit"], "aav_compatible": False},
    "PE2":       {"pam": "NGG",   "size_aa": 1700, "editing_types": ["prime-edit"], "aav_compatible": False},
}

AAV_SIZE_LIMIT_AA = 1200
