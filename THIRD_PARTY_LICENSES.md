# Third-party software and data

The ORACLE workflow code in this repository is licensed under the
[MIT License](LICENSE). The workflow installs and invokes separate
third-party packages through Conda/Bioconda; it does not relicense those
packages. Reference genomes, annotations, blacklists, and motif databases also
retain their source-specific terms.

This file is a practical inventory, not legal advice. Users distributing an
environment, container, reference bundle, or derived database should review
the exact package build and upstream terms themselves.

## Main command-line packages

The license labels below follow the metadata of the declared Conda/Bioconda
packages used by the workflow at the time this file was reviewed. Environment
specifications pin or constrain compatibility-sensitive tools; they are not
immutable lock files.

| Package | Workflow role | Package metadata license |
|---|---|---|
| fastp | read trimming and QC | MIT |
| FastQC | raw-read QC | GPL-3.0-or-later |
| MultiQC | combined report | GPL-3.0-or-later |
| Bowtie2 | paired-end alignment | GPL-3.0-or-later |
| SAMtools | BAM processing | MIT |
| Picard | duplicate marking and metrics | MIT |
| deepTools | Tn5 shifting, coverage, and QC plots | MIT |
| MACS3 | peak calling | BSD-3-Clause |
| Genrich | optional condition-level peak cross-check | MIT |
| Subread/featureCounts | peak quantification | GPL-3.0-only |
| SRA Toolkit | SRA download and conversion | Public Domain in package metadata; review bundled notices |
| HOMER 4.11 | motif enrichment | GNU GPL v3 in the [Bioconda package metadata](https://anaconda.org/bioconda/homer) |

HOMER is installed from the pinned Bioconda package and receives the configured
FASTA directly. The workflow does not install or redistribute HOMER genome
packages. If a different HOMER build or auxiliary dataset is used, review that
build's metadata and data terms separately.

The optional TOBIAS environment and all transitive dependencies retain their
own licenses. Inspect the resolved Conda package metadata for the exact build
used in a run.

## R, Bioconductor, and Python packages

DESeq2, ashr, ChIPseeker, clusterProfiler, chromVAR, JASPAR2020, genome
annotation packages, pandas, pysam, plotting libraries, Snakemake, and their
dependencies are installed as separate packages. Their licenses are not
uniform. The authoritative record for a particular run is the resolved package
metadata plus each upstream project's license.

## JASPAR and other motif databases

[JASPAR CORE](https://jaspar.elixir.no/about/) data are provided under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) and require
attribution. This applies to JASPAR motif data independently of the license of
the R package or software used to read them.

TOBIAS accepts user-supplied motif databases. A file's syntax does not determine
its license: a motif file in MEME format does **not** acquire the MEME Suite
software license merely because it uses that format. Cite and follow the terms
of the database from which the motifs were obtained. The workflow does not run
the MEME Suite.

## Reference genomes, annotations, and blacklists

Reference data are downloaded or supplied by the user and are not covered by
this repository's MIT license. Record and review the terms for the exact
versions used, including:

- Ensembl or another FASTA/GTF provider;
- the Boyle Lab/ENCODE blacklist source;
- TxDb, OrgDb, BSgenome, and related annotation data;
- JASPAR or another motif collection;
- GO and KEGG resources used for enrichment.

Do not assume that a data resource is unrestricted because it can be downloaded
automatically. Preserve required attribution and any version-specific notices
with the analysis provenance.
