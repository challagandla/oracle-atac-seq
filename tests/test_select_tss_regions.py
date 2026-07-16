import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "select_tss_regions.py"
SPEC = importlib.util.spec_from_file_location("select_tss_regions", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def write_regions(path, count=10):
    path.write_text(
        "".join(f"chr{index + 1}\t{index * 10}\t{index * 10 + 1}\n"
                for index in range(count)),
        encoding="utf-8",
    )


def test_bounded_selection_is_deterministic_and_spans_the_full_bed(tmp_path):
    source = tmp_path / "all.bed"
    output = tmp_path / "selected.bed"
    write_regions(source)

    assert MOD.select_regions(source, output, 3) == (3, 10)
    assert output.read_text(encoding="utf-8").splitlines() == [
        "chr2\t10\t11",
        "chr6\t50\t51",
        "chr9\t80\t81",
    ]


def test_zero_retains_every_region(tmp_path):
    source = tmp_path / "all.bed"
    output = tmp_path / "selected.bed"
    write_regions(source, count=4)

    assert MOD.select_regions(source, output, 0) == (4, 4)
    assert output.read_text(encoding="utf-8") == source.read_text(encoding="utf-8")


@pytest.mark.parametrize("maximum", [-1, -10])
def test_negative_limits_fail(tmp_path, maximum):
    source = tmp_path / "all.bed"
    write_regions(source)

    with pytest.raises(ValueError, match="zero or a positive"):
        MOD.select_regions(source, tmp_path / "selected.bed", maximum)
