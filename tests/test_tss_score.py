from pathlib import Path
import importlib.util

import numpy as np


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "tss_score.py"
SPEC = importlib.util.spec_from_file_location("tss_score", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def test_enrichment_uses_distal_edge_background_and_profile_maximum():
    # Forty 10 bp bins: the first and last 100 bp define a background of 2.
    profile = np.full(40, 2.0)
    profile[20] = 14.0
    assert MOD.enrichment(profile, flank_bins=10) == 7.0


def test_zero_background_is_reported_as_zero_not_infinite():
    profile = np.zeros(40)
    profile[20] = 10.0
    assert MOD.enrichment(profile, flank_bins=10) == 0.0
