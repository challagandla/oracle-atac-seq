#!/usr/bin/env python3
"""The pipeline's colour palette, in one place.

Every Python figure and every deepTools call takes its colours from here. The R
scripts read the same values from `atac_theme.R`; `tests/test_palette.py` fails
if the two ever drift apart.

Why this file exists: a hue has to mean the same thing in every panel of the
report. A hue is bound to a *condition*, never to a sample's position in the
sample sheet, so adding a replicate or renaming a sample cannot repaint the
groups. Conditions are sorted before assignment, which makes the mapping a pure
function of the condition names -- the rule `atac_condition_colours()` already
follows in atac_theme.R.

The categorical hues are handed out in fixed order and never cycled. They are
colour-blind-safe: across the first four slots the worst all-pairs separation is
dE 24.2 under protanopia, against a floor of 12.
"""
from __future__ import annotations

import shlex

# Colour-blind-safe categorical hues, fixed slot order (mirrors ATAC_CATEGORICAL
# in atac_theme.R). Slots are handed out in order, never recycled.
CATEG = ["#2a78d6", "#1baf7a", "#eda100", "#008300", "#4a3aa7", "#e34948",
         "#e87ba4", "#eb6834", "#6d4b9f", "#00a2b3", "#8c8c00", "#a6611a"]

# Sequential ramp for magnitude (one hue, light -> dark), by deepTools' name.
SEQUENTIAL_CMAP = "Blues"


def condition_colours(conditions):
    """Map each distinct condition to a hue. Sorted, so row order cannot matter.

    Returns {condition: hex}. Beyond 12 conditions a categorical legend stops
    being readable anyway, so this refuses rather than inventing hues that have
    never been checked for colour-blind separation.
    """
    uniq = sorted(set(conditions))
    if len(uniq) > len(CATEG):
        raise ValueError(
            f"{len(uniq)} conditions but only {len(CATEG)} categorical hues. "
            "Collapse the rare groups, or disable categorical report overlays."
        )
    return {c: CATEG[i] for i, c in enumerate(uniq)}


def colours_for_samples(samples, cond_of):
    """One hue per sample, taken from that sample's condition.

    Replicates of a condition deliberately share a hue: on an overlay plot the
    replicate lines should superimpose, and that is the check a reader makes.
    `cond_of` maps sample -> condition.
    """
    colour = condition_colours(cond_of[s] for s in samples)
    return [colour[cond_of[s]] for s in samples]


def shell_colours(samples, cond_of):
    """`colours_for_samples`, quoted for a POSIX shell.

    A bare `#2a78d6` on a command line starts a comment: the shell drops it and
    the rest of the line, so the flag it belonged to loses its argument and any
    trailing redirect vanishes with it. Anything interpolated into a `shell:`
    block has to come through here.
    """
    return " ".join(shlex.quote(c) for c in colours_for_samples(samples, cond_of))
