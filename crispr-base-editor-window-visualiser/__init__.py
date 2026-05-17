"""
CRISPR Base Editor Window Visualiser — public API.

Programmatic usage:

    from crispr_base_editor_window_visualiser.src.config import EDITOR_PROFILES
    from crispr_base_editor_window_visualiser.src.analyser import analyse_sequence
    from crispr_base_editor_window_visualiser.src.pipeline import run_pipeline

    editor = EDITOR_PROFILES["ABE8e"]
    result = run_pipeline(
        guide_rna  = "GCACTGACCTGAGTTCAGTG",
        target_dna = "GCACTGACCTGAGTTCAGTGNGG",
        editor     = editor,
        output_prefix = "output/my_guide",
    )
"""

from src.config import (  # noqa: F401
    EDITOR_PROFILES,
    BaseEditorProfile,
    build_custom_profile,
)
from src.analyser import analyse_sequence, PositionInfo  # noqa: F401
from src.pipeline import run_pipeline                    # noqa: F401
