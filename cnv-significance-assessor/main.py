#!/usr/bin/env python3
"""
CNV Significance Assessor — entry point.

Run from the project root (cnv-significance-assessor/):

    python main.py variants.bed annotation.gff3 --output-dir results/
    python main.py --help
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from pipeline import main  # noqa: E402

sys.exit(main())
