"""Library-complexity fragment identity invariants."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pysam
import pytest


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "library_complexity.py"
SPEC = importlib.util.spec_from_file_location("library_complexity", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def read1(*, start, mate_start, template_length, reverse=False):
    read = pysam.AlignedSegment()
    read.flag = 83 if reverse else 99
    read.reference_id = 0
    read.reference_start = start
    read.next_reference_id = 0
    read.next_reference_start = mate_start
    read.template_length = template_length
    read.mapping_quality = 60
    return read


def test_fragment_key_uses_template_bounds_not_reverse_alignment_start():
    # Both records describe chr1:100-250. The reverse read's aligned left edge
    # differs after trimming, but its fragment 5' boundary does not.
    first = read1(start=200, mate_start=100, template_length=-150, reverse=True)
    trimmed = read1(start=205, mate_start=100, template_length=-150, reverse=True)
    assert MOD.fragment_key(first) == (100, 250, True)
    assert MOD.fragment_key(trimmed) == MOD.fragment_key(first)


def test_fragment_key_rejects_cross_reference_or_zero_length_pairs():
    cross = read1(start=100, mate_start=200, template_length=150)
    cross.next_reference_id = 1
    zero = read1(start=100, mate_start=100, template_length=0)
    assert MOD.fragment_key(cross) is None
    assert MOD.fragment_key(zero) is None


def test_complexity_metrics_retain_their_auditable_integer_components():
    nrf, pbc1, pbc2 = MOD.complexity_metrics(
        total=10, distinct=7, singletons=5, doubletons=1
    )
    assert nrf == pytest.approx(0.7)
    assert pbc1 == pytest.approx(5 / 7)
    assert pbc2 == pytest.approx(5.0)

    source = SCRIPT.read_text()
    assert "M1_singletons\\tM2_doubletons\\t" in source
    assert "NRF\\tPBC1\\tPBC2" in source
