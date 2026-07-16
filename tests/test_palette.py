"""Palette invariants shared by shell commands and cross-sample figures."""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from palette import (  # noqa: E402
    CATEG,
    colours_for_samples,
    condition_colours,
    shell_colours,
)

THEME_R = ROOT / "workflow" / "scripts" / "atac_theme.R"


def test_python_and_r_palettes_are_identical():
    """atac_theme.R is the palette the R figures use. It must not drift."""
    text = THEME_R.read_text()
    block = re.search(r"ATAC_CATEGORICAL <- c\((.*?)\)", text, re.S)
    assert block, "ATAC_CATEGORICAL not found in atac_theme.R"
    r_hues = re.findall(r"#[0-9a-fA-F]{6}", block.group(1))
    assert r_hues == CATEG


def test_hue_follows_condition_not_row_order():
    """Reordering the sample sheet must not repaint a group."""
    cond_of = {"a1": "alpha", "a2": "alpha", "b1": "beta", "b2": "beta"}
    forward = colours_for_samples(["a1", "a2", "b1", "b2"], cond_of)
    shuffled = colours_for_samples(["b2", "a1", "b1", "a2"], cond_of)
    assert forward == [CATEG[0], CATEG[0], CATEG[1], CATEG[1]]
    assert shuffled == [CATEG[1], CATEG[0], CATEG[1], CATEG[0]]


def test_replicates_of_a_condition_share_a_hue():
    cond_of = {"x_rep1": "x", "x_rep2": "x", "y_rep1": "y"}
    hues = colours_for_samples(["x_rep1", "x_rep2", "y_rep1"], cond_of)
    assert hues[0] == hues[1] != hues[2]


def test_dropping_a_condition_does_not_repaint_the_survivors():
    """A filter that changes the group count must leave the others alone."""
    assert condition_colours(["alpha", "beta"])["alpha"] == CATEG[0]
    assert condition_colours(["alpha"])["alpha"] == CATEG[0]


def test_condition_order_is_independent_of_input_order():
    assert condition_colours(["z", "a"]) == condition_colours(["a", "z"])


def test_too_many_conditions_refuses_rather_than_inventing_hues():
    with pytest.raises(ValueError, match="only 12 categorical hues"):
        condition_colours([f"c{i}" for i in range(13)])


def test_shell_colours_are_quoted_so_the_shell_cannot_eat_them():
    """A bare hash starts a shell comment and truncates the command line."""
    cond_of = {"s1": "a", "s2": "b"}
    out = shell_colours(["s1", "s2"], cond_of)
    assert out == f"'{CATEG[0]}' '{CATEG[1]}'"
    assert not re.search(r"(^|\s)#", out), "an unquoted # would open a shell comment"


def test_shell_colours_survive_a_real_shell():
    """Round-trip through bash: the quoted colours must arrive as separate argv."""
    import subprocess

    cond_of = {"s1": "a", "s2": "b"}
    cmd = f"printf '%s\\n' {shell_colours(['s1', 's2'], cond_of)}"
    out = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, check=True)
    assert out.stdout.split() == [CATEG[0], CATEG[1]]


def test_unquoted_hex_really_does_get_eaten_by_the_shell():
    """Guards the premise. If this ever fails, the quoting above is cargo cult."""
    import subprocess

    cmd = f"printf '%s\\n' {CATEG[0]} {CATEG[1]}"
    out = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "", "expected the shell to swallow bare hex as a comment"
