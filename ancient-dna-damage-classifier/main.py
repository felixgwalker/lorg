#!/usr/bin/env python3
"""
Ancient DNA Damage Classifier — entry point.

Run from the project root (ancient-dna-damage-classifier/):

    python main.py reads.bam --output-dir results/
    python main.py --help
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from pipeline import main  # noqa: E402

sys.exit(main())
