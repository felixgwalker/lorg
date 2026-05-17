#!/usr/bin/env python3
"""
CRISPR Base Editor Window Visualiser — entry point.

Run from the project root (crispr-base-editor-window-visualiser/):

    python main.py GUIDE TARGET EDITOR [options]
    python main.py --help

Examples:
    python main.py GCACTGACCTGAGTTCAGTG GCACTGACCTGAGTTCAGTGNGG ABE8e
    python main.py GCACTGACCTGAGTTCAGTG GCACTGACCTGAGTTCAGTGNGG BE4max -o results/my_guide
    python main.py GCACTGACCTGAGTTCAGTG GCACTGACCTGAGTTCAGTGNGG custom \\
        --editor-class CBE --window-start 4 --window-end 7
    python main.py --list-editors
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from pipeline import main  # noqa: E402

sys.exit(main())
