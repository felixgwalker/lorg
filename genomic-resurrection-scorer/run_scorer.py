#!/usr/bin/env python3
"""
Top-level entry point for the Genomic Resurrection Scorer pipeline.

Run from the project root (genomic-resurrection-scorer/):

    python run_scorer.py data/thylacine_case_study/metrics.json

    # Write report to the frontend public folder so the Next.js dashboard
    # can display it automatically:
    python run_scorer.py data/thylacine_case_study/metrics.json \\
        --output src/frontend/public/reports/report.json
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from pipeline.scorer import main  # noqa: E402

sys.exit(main())
