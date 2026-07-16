from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "gtf_to_tss.py"
SPEC = importlib.util.spec_from_file_location("gtf_to_tss", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def gene(chrom: str, start: int, end: int, strand: str, gene_id: str) -> str:
    return (
        f'{chrom}\ttest\tgene\t{start}\t{end}\t.\t{strand}\t.\t'
        f'gene_id "{gene_id}";\n'
    )


def test_rejects_gtf_without_gene_features(tmp_path):
    chrom = tmp_path / "chrom.sizes"
    gtf = tmp_path / "genes.gtf"
    chrom.write_text("chr1\t1000\n")
    gtf.write_text(
        'chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1";\n'
    )

    with pytest.raises(ValueError, match="no gene features"):
        MOD.extract_tss(gtf, chrom)


def test_harmonizes_chr_prefix_and_mitochondrial_aliases(tmp_path):
    chrom = tmp_path / "chrom.sizes"
    gtf = tmp_path / "genes.gtf"
    chrom.write_text("chr1\t1000\nchrM\t100\n")
    gtf.write_text(
        gene("1", 101, 200, "+", "plus")
        + gene("MT", 20, 80, "-", "minus")
    )

    records, stats = MOD.extract_tss(gtf, chrom)

    assert records == [
        ("chr1", 100, 101, "plus", ".", "+"),
        ("chrM", 79, 80, "minus", ".", "-"),
    ]
    assert stats["renamed"] == 2
    assert stats["retained_fraction"] == 1


def test_rejects_severe_partial_contig_mismatch(tmp_path):
    chrom = tmp_path / "chrom.sizes"
    gtf = tmp_path / "genes.gtf"
    chrom.write_text("chr1\t1000\n")
    mapped = "".join(
        gene("chr1", 10 + index * 20, 20 + index * 20, "+", f"g{index}")
        for index in range(9)
    )
    gtf.write_text(mapped + gene("other", 10, 20, "+", "missing"))

    with pytest.raises(ValueError, match="refusing biased TSS QC"):
        MOD.extract_tss(gtf, chrom)


def test_accepts_exactly_ninety_five_percent_with_warning(tmp_path):
    chrom = tmp_path / "chrom.sizes"
    gtf = tmp_path / "genes.gtf"
    chrom.write_text("chr1\t1000\n")
    mapped = "".join(
        gene("chr1", 10 + index * 20, 20 + index * 20, "+", f"g{index}")
        for index in range(19)
    )
    gtf.write_text(mapped + gene("other", 10, 20, "+", "missing"))

    records, stats = MOD.extract_tss(gtf, chrom)

    assert len(records) == 19
    assert stats["retained_fraction"] == 0.95
    assert stats["warn"] == 1


@pytest.mark.parametrize(
    "strand,start,end",
    [("+", 1001, 1100), ("-", 900, 1001)],
)
def test_rejects_out_of_bounds_tss(tmp_path, strand, start, end):
    chrom = tmp_path / "chrom.sizes"
    gtf = tmp_path / "genes.gtf"
    chrom.write_text("chr1\t1000\n")
    gtf.write_text(gene("chr1", start, end, strand, "outside"))

    with pytest.raises(ValueError, match="0/1 gene TSSs"):
        MOD.extract_tss(gtf, chrom)
