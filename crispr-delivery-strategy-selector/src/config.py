"""Configuration for CRISPR Delivery Strategy Selector."""

DELIVERY_STRATEGIES: dict[str, dict] = {
    "plasmid": {
        "description": "Plasmid transfection (lipofection or electroporation)",
        "max_payload_kb": 10,
        "transient": True,
        "integration_risk": "low",
        "primary_cell_score": 0.30,
        "dividing_cell_score": 0.80,
    },
    "rnp": {
        "description": "Ribonucleoprotein electroporation",
        "max_payload_kb": None,
        "transient": True,
        "integration_risk": "none",
        "primary_cell_score": 0.75,
        "dividing_cell_score": 0.80,
    },
    "lnp": {
        "description": "Lipid nanoparticle (mRNA + sgRNA)",
        "max_payload_kb": 5,
        "transient": True,
        "integration_risk": "none",
        "primary_cell_score": 0.60,
        "dividing_cell_score": 0.65,
    },
    "aav": {
        "description": "Adeno-associated virus",
        "max_payload_kb": 4.7,
        "transient": True,
        "integration_risk": "very_low",
        "primary_cell_score": 0.85,
        "dividing_cell_score": 0.50,
    },
    "lentivirus": {
        "description": "Lentiviral transduction",
        "max_payload_kb": 8,
        "transient": False,
        "integration_risk": "high",
        "primary_cell_score": 0.70,
        "dividing_cell_score": 0.75,
    },
}
