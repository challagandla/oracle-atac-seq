# Third-party tools, data & licenses

This pipeline's own code is **MIT** (see `LICENSE`). It **bundles no third-party source code and no
reference data** — tools are installed via conda/bioconda and invoked as separate processes, and
reference genomes/annotations/blacklists are downloaded by the user. (Files under `.test/` are tiny
synthetic fixtures, e.g. a 9-byte `blacklist.bed`, not real data.) Invoking a separate program is not
a derivative work, so the MIT license is unaffected by the licenses of the tools it orchestrates.

## ⚠️ Restriction that affects commercial use

| Dependency | Used by | Terms |
|---|---|---|
| **HOMER** | motif enrichment (`workflow/rules/motif.smk`, `workflow/envs/motif.yaml`) | HOMER is **freeware for academic / non-profit use, is not open-source, and may not be redistributed**; commercial use requires contacting the author (C. Benner, Salk/UCSD). The HOMER genome packages installed via `configureHomer.pl` are derived from **UCSC** data (academic/non-profit; commercial needs a UCSC license). |

HOMER is conda-installed and invoked (not bundled), so this repository redistributes nothing. But the
**motif step is academic/non-commercial**. For commercial use, obtain HOMER/UCSC permissions or skip
the HOMER motif step (peak calling, differential accessibility, and JASPAR-based footprinting do not
depend on HOMER).

## Tools (installed via conda; invoked, not bundled)

| Tool | License | Role |
|---|---|---|
| Bowtie2 / BWA | GPL-3.0 | alignment |
| samtools | MIT | |
| MACS3 | BSD-3-Clause | peak calling |
| Genrich | MIT | peak calling |
| deepTools | GPL-3.0 | coverage/QC |
| Picard | MIT | dedup/metrics |
| **HOMER** | **academic/non-profit; not redistributable** (see above) | motif enrichment |
| sra-tools | Public Domain (US Gov) | SRA download |
| JASPAR motif matrices | CC0 | footprinting input |

If footprinting is configured with a **MEME-format** motif database instead of JASPAR, note the **MEME
Suite is free for non-commercial use only** (commercial use requires a license).

**R/Bioconductor:** DESeq2 (LGPL ≥3), ChIPseeker/annotation packages (Artistic-2.0/GPL), ggplot2 (MIT).
Invoked within R; not redistributed.

**On the GPL tools:** called as independent executables (Snakemake rules) — mere aggregation, not
linking — so they impose no copyleft obligation on this MIT pipeline.

## Reference data (downloaded at run time; not redistributed)
Genome/annotation from GENCODE/Ensembl (open); ENCODE blacklist from the Boyle Lab (open, cite). None
are bundled in this repository.

## Bottom line
No code-incorporation conflict, no bundled code/data. The one operative constraint is **HOMER**
(academic/non-profit; not redistributable; commercial use needs author/UCSC permission). Everything
else is freely usable with citation.
