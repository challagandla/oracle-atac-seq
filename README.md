# ORACLE ATAC-seq

[![CI](https://github.com/challagandla/oracle-atac-seq/actions/workflows/ci.yml/badge.svg)](https://github.com/challagandla/oracle-atac-seq/actions/workflows/ci.yml)

ORACLE ATAC-seq is a Snakemake workflow for reproducible analysis of
paired-end bulk ATAC-seq data. It starts from local FASTQ files or SRA run
accessions and produces filtered alignments, accessibility tracks, peaks,
quality-control reports, a consensus count matrix, differential-accessibility
results, annotations, motifs, and optional transcription-factor analyses.

The workflow is **ENCODE-informed**, not a drop-in implementation of the
official ENCODE ATAC-seq pipeline. In particular, it creates fixed-width peaks
and requires within-condition replicate support; it does **not** implement IDR.
Use the [official ENCODE pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline)
when exact ENCODE processing or submission compliance is required.

## Scope

- Paired-end bulk ATAC-seq
- Linux workstations, servers, and WSL2
- Local FASTQs, SRA run accessions, or a mixture of both
- Human, mouse, rat, or carefully configured custom references
- Two-group and covariate-aware differential designs

Single-cell ATAC-seq, single-end libraries, IDR, spike-in normalization, and
clinical interpretation are outside the current scope.

## What the workflow does

```text
FASTQ or SRA
  -> FastQC and fastp
  -> Bowtie2 paired-end alignment
  -> proper-pair, MAPQ, mitochondrial, duplicate, and blacklist filtering
  -> fragment-level and Tn5 cut-site CPM bigWigs
  -> per-sample MACS3 peaks
  -> fixed-width, within-condition replicate-supported consensus peaks
  -> per-library Rsubread featureCounts and checked peak-by-sample matrix merge
  -> DESeq2 differential accessibility
  -> ChIPseeker annotation and clusterProfiler enrichment
  -> HOMER motif enrichment
  -> optional chromVAR and TOBIAS
  -> MultiQC plus publication-oriented QC and differential figures
```

Genrich can also pool replicates within each condition as a secondary peak
calling cross-check. The primary consensus/count matrix is built from
per-sample MACS3 peaks.

## Quick start

### 1. Clone and install

```bash
git clone https://github.com/challagandla/oracle-atac-seq.git
cd oracle-atac-seq
bash setup.sh
```

The setup script creates the small Snakemake runner environment and the
rule-specific Conda environments. Shell activation is not required. Verify an
existing installation with:

```bash
bash setup.sh --check
```

### 2. Create project configuration copies

The tracked files are reusable templates. Copy them before entering project
metadata or paths; never put project records into the tracked templates:

```bash
cp config/samples.tsv config/samples.project.tsv
cp config/project.example.yaml config/project_overrides.yaml
```

Both copies are ignored by Git. Edit `config/samples.project.tsv`. Its first six
columns are required; `include` is optional, and additional design columns such
as `batch` or `donor` are allowed. A false `include` value keeps a library
documented but excludes it from the analysis.

Both `sample` and `condition` must match
`^[A-Za-z0-9][A-Za-z0-9._-]*$`: start with a letter or digit, followed only by
letters, digits, dots, underscores, or hyphens.

Local read paths must end in `.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz`.

| sample | condition | replicate | fq1 | fq2 | sra | include |
|---|---|---:|---|---|---|---|
| control_r1 | control | 1 | data/control_r1_R1.fq.gz | data/control_r1_R2.fq.gz | | true |
| control_r2 | control | 2 | data/control_r2_R1.fq.gz | data/control_r2_R2.fq.gz | | true |
| treatment_r1 | treatment | 1 | data/treatment_r1_R1.fq.gz | data/treatment_r1_R2.fq.gz | | true |
| treatment_r2 | treatment | 2 | data/treatment_r2_R1.fq.gz | data/treatment_r2_R2.fq.gz | | true |

For SRA input, leave `fq1` and `fq2` empty and provide an SRR, ERR, or DRR run
accession in `sra`. Each row must provide local paired FASTQs or one SRA run,
never both.

### 3. Configure the analysis

Edit `config/project_overrides.yaml`, set its sample-sheet path, and review at
least:

```yaml
samples: "config/samples.project.tsv"
```

Replace `my_project` in that copied overlay's five output/reference paths so
each analysis starts in an isolated namespace.

- `genome.build`, reference paths, and blacklist assembly
- `diffacc.design` and `diffacc.contrast`
- `peaks.consensus_min_replicates` and `peaks.consensus_peak_width`
- optional annotation, enrichment, motif, chromVAR, and footprinting switches
- CPU settings under `resources`

Positive differential log2 fold change means greater accessibility in the
contrast numerator:

```yaml
diffacc:
  design: "~condition"
  contrast: ["condition", "treatment", "control"]
```

Use biological replicates. The workflow requires at least two included
libraries per condition for differential analysis; three or more are strongly
preferred. Preflight also requires every included condition to have at least
`peaks.consensus_min_replicates` libraries so each condition can contribute
condition-specific consensus peaks.

For `genome.build: "custom"`, provide local FASTA and GTF paths or a complete
Ensembl species/release/assembly tuple. Also configure the effective and MACS
genome sizes. At least 95% of GTF gene TSSs must map within the configured
FASTA after conservative name reconciliation. Enrichment requires a matching taxonomic ID and OrgDb; chromVAR
requires that taxonomic ID and a matching BSgenome. Add the named R packages to
`workflow/envs/r.yaml`; annotation can use the configured TxDb or GTF fallback.

### 4. Preview, review QC, then run

```bash
bash run.sh --dry-run \
  --configfile config/config.yaml config/project_overrides.yaml
bash run.sh --cores 8 qc_review \
  --configfile config/config.yaml config/project_overrides.yaml
```

Open `results/my_project/qc/qc_review_report.html`, record technical inclusion
or exclusion decisions without consulting differential outcomes, and freeze
the sample sheet. Then repeat the dry-run and launch the complete analysis:

```bash
bash run.sh --dry-run \
  --configfile config/config.yaml config/project_overrides.yaml
bash run.sh --cores 8 \
  --configfile config/config.yaml config/project_overrides.yaml
```

The dedicated `qc_review` target has an exact manifest that excludes consensus
counts, DESeq2, enrichment, motif, footprint, and chromVAR results. The final
`multiqc_report.html` is a combined post-analysis report, not the outcome-blind
inclusion gate.

Run one target when a full run is not needed:

```bash
bash run.sh --cores 8 \
  results/my_project/qc/multiqc_report.html \
  --configfile config/config.yaml config/project_overrides.yaml
bash run.sh --cores 8 \
  results/my_project/consensus/consensus_peaks.bed \
  --configfile config/config.yaml config/project_overrides.yaml
```

The workflow is restartable. Re-running the same command schedules only
missing, incomplete, or outdated work. Keep this base-then-overlay order for
dry-runs, production runs, targets, restarts, and summaries.

## Start with quality control

Open `results/my_project/qc/qc_review_report.html` before building or
interpreting differential results. Review every library for:

- read quality and adapter content
- mapping and usable-fragment yield
- duplicate rate and library complexity
- mitochondrial fraction
- fragment-size periodicity
- TSS enrichment
- FRiP
- outcome-blind genome-bin PCA, sample correlation, and replicate similarity

Do not exclude a sample solely because it weakens a desired contrast. Record a
technical reason, preserve the original row, set its `include` cell to `false`
in `config/samples.project.tsv`, and rerun the dry-run before analysis.

## Main outputs

| Path | Meaning |
|---|---|
| `results/my_project/qc/qc_review_report.html` | Outcome-blind library QC report used before inclusion is frozen |
| `results/my_project/qc/multiqc_report.html` | Final combined QC and analysis-figure report |
| `results/my_project/qc/summary/qc_decisions_mqc.tsv` | Transparent per-library QC metrics and review flags |
| `results/my_project/filtered/` | Analysis-ready BAM files |
| `results/my_project/coverage/` | Fragment-level and Tn5 cut-site CPM bigWig tracks |
| `results/my_project/peaks/macs3/` | Per-sample narrowPeak files |
| `results/my_project/peaks/genrich/` | Optional per-condition Genrich peaks |
| `results/my_project/consensus/consensus_peaks.bed` | Replicate-supported consensus universe |
| `results/my_project/counts/per_sample/` | Independently restartable featureCounts tables and assignment summaries |
| `results/my_project/counts/consensus_counts.tsv` | Peak-by-sample fragment counts |
| `results/my_project/diffacc/` | DESeq2 table, finite-p-value tested peak set, normalized counts, and figures |
| `results/my_project/annotation/` | Genomic annotation of consensus peaks |
| `results/my_project/enrichment/` | GO and optional KEGG over-representation results |
| `results/my_project/motif/` | HOMER results for more-open and less-open peaks |
| `results/my_project/chromvar/` | Optional motif deviation results |
| `results/my_project/footprint/` | Optional TOBIAS outputs from condition-merged final filtered BAMs |
| `results/my_project/figures/` | Vector fragment-size, TSS, and outcome-blind sample-similarity figures |
| `results/my_project/provenance/` | Effective config, sample copy, raw-input/reference/environment hashes, and run manifest |
| `results/my_project/provenance/project_record/` | Project-managed base/overlay/sample copies, hashes, software check, and run summary; see the tutorial |
| `logs/my_project/` | Per-rule diagnostic logs |

These paths match the copied project overlay. Substitute the project name if
you changed it; custom root settings replace the corresponding prefixes.

## Standards boundary

ORACLE follows widely used bulk ATAC-seq practices: paired-end mapping,
high-confidence filtering, duplicate and blacklist removal, Tn5-aware
processing, ATAC-specific QC, replicate-aware peaks, count-based differential
analysis, and reproducible per-rule software environments.

It should still be reviewed for each study. Reference FASTA, GTF, chromosome
names, blacklist, and annotation packages must describe the same assembly.
There is no official Boyle Lab mm39 blacklist; do not apply an mm10 blacklist
to mm39 coordinates without a documented, validated liftover. Thresholds are
guides rather than substitutes for experimental context.

The consensus method is intentionally transparent: peaks are converted to
fixed-width summit-centered intervals, retained when supported by the
configured number of biological replicates **within at least one condition**,
then combined across conditions. This preserves condition-specific peaks, but
it is not the ENCODE IDR procedure and does not calculate rescue or
self-consistency ratios.

## Documentation

- [Complete beginner-to-advanced tutorial](TUTORIAL.md)
- [Contributing guide](CONTRIBUTING.md)
- [Third-party licenses and data terms](THIRD_PARTY_LICENSES.md)
- [Citation metadata](CITATION.cff)

## License and citation

The workflow code is available under the [MIT License](LICENSE). Tools and
reference data retain their own licenses and terms; see
[THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md).

If you use the workflow, cite the software repository and the methods relevant
to the stages you report. The tutorial provides a curated reference list.
