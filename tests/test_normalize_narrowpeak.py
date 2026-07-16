import importlib.util
import stat
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "normalize_narrowpeak.py"
SPEC = importlib.util.spec_from_file_location("normalize_narrowpeak", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def record(*, score="1170", summit="25"):
    return f"chr1\t100\t200\tpeak_1\t{score}\t.\t42.5\t12\t9\t{summit}\n"


def test_normalizes_only_the_ucsc_display_score(tmp_path):
    path = tmp_path / "peaks.narrowPeak"
    path.write_text(record(), encoding="utf-8")
    path.chmod(0o640)

    assert MOD.normalize_file(path, path) == 1
    fields = path.read_text(encoding="utf-8").rstrip().split("\t")
    assert fields == [
        "chr1", "100", "200", "peak_1", "1000", ".",
        "42.5", "12", "9", "25",
    ]
    assert stat.S_IMODE(path.stat().st_mode) == 0o640


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("chr1\t100\t200\tpeak\t1\t.\t2\t3\t4\n", "exactly 10"),
        (record(score="not-a-number"), "score must be an integer"),
        (record().replace("\t100\t200\t", "\t100.0\t200\t"), "chromStart must"),
        (record().replace("\t100\t200\t", "\t1e2\t200\t"), "chromStart must"),
        (record().replace("\t100\t200\t", "\t 100\t200\t"), "chromStart must"),
        (record(summit="100"), "lies outside"),
        (record().replace("\t42.5\t", "\t-0.1\t"), "signalValue must"),
        (record().replace("\t12\t9\t", "\t-2\t9\t"), "pValue must"),
        ("", "no records"),
    ],
)
def test_rejects_malformed_or_empty_files_without_replacing_input(
    tmp_path, text, message
):
    path = tmp_path / "peaks.narrowPeak"
    path.write_text(text, encoding="utf-8")

    with pytest.raises(ValueError, match=message):
        MOD.normalize_file(path, path)
    assert path.read_text(encoding="utf-8") == text
