"""Configuration for Biosafety Risk Assessor."""

SELECT_AGENT_DB = "select_agents.fa"
VIRULENCE_FACTOR_DB = "vfdb.fa"
ANTIBIOTIC_RESISTANCE_DB = "card.fa"
TOXIN_DB = "toxprot.fa"

BLAST_EVALUE = 1e-5
MIN_IDENTITY = 50.0
MIN_COVERAGE = 50.0

HIGH_RISK_SCORE_THRESHOLD = 0.7
MEDIUM_RISK_SCORE_THRESHOLD = 0.3
