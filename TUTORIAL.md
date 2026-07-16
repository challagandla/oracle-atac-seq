# ORACLE ATAC-seq: complete beginner-to-advanced tutorial

This tutorial explains the biology, study design, installation, configuration,
execution, quality control, interpretation, and limitations of the ORACLE
paired-end bulk ATAC-seq workflow. Read the design and QC sections before
starting a production study; a successful command is not by itself evidence
that the resulting biology is reliable.

## Contents

1. [What ATAC-seq measures](#1-what-atac-seq-measures)
2. [What this workflow is, and is not](#2-what-this-workflow-is-and-is-not)
3. [Design the experiment before sequencing](#3-design-the-experiment-before-sequencing)
4. [Prepare the inputs](#4-prepare-the-inputs)
5. [Install the workflow](#5-install-the-workflow)
6. [Choose and verify the reference genome](#6-choose-and-verify-the-reference-genome)
7. [Configure a run](#7-configure-a-run)
8. [Preview and run](#8-preview-and-run)
9. [What happens at each stage](#9-what-happens-at-each-stage)
10. [Quality control: make a decision before biology](#10-quality-control-make-a-decision-before-biology)
11. [Understand the main outputs](#11-understand-the-main-outputs)
12. [Differential accessibility done carefully](#12-differential-accessibility-done-carefully)
13. [Motifs, chromVAR, and TOBIAS](#13-motifs-chromvar-and-tobias)
14. [Restart, recover, and update safely](#14-restart-recover-and-update-safely)
15. [Provenance and reproducibility](#15-provenance-and-reproducibility)
16. [Important limitations](#16-important-limitations)
17. [Troubleshooting](#17-troubleshooting)
18. [Methods and citations](#18-methods-and-citations)

## Quick QC glossary

- **FRiP**: fraction of final usable fragments that overlap the library's
  MACS3 peaks.
- **TSS enrichment**: enrichment of Tn5 insertions around transcription start
  sites relative to flanking background.
- **NRF**: non-redundant fraction, or distinct pre-deduplication fragment
  locations divided by all fragments considered.
- **PBC1**: PCR bottleneck coefficient 1, or locations observed once divided
  by all distinct fragment locations.
- **PBC2**: PCR bottleneck coefficient 2, or locations observed once divided
  by locations observed exactly twice.

The complexity table retains the integer M1 singleton and M2 doubleton counts,
so NRF, PBC1, and PBC2 can be recalculated after temporary pre-deduplication
BAMs are removed.

## 1. What ATAC-seq measures

ATAC-seq uses the Tn5 transposase to insert sequencing adapters preferentially
into accessible DNA. After paired-end sequencing, the genomic positions and
lengths of the resulting fragments provide complementary information:

- short fragments are enriched in nucleosome-free regulatory regions;
- fragments near one nucleosome in length report chromatin organization;
- local fragment pileups identify accessible peaks;
- changes in fragment counts across conditions identify differential
  accessibility;
- sequence motifs and footprint-like patterns can suggest candidate
  transcription-factor regulators.

ATAC-seq measures accessibility, not transcription-factor binding, enhancer
activity, or gene expression directly. Motif enrichment, nearest-gene
annotation, chromVAR deviations, and footprints are hypotheses that should be
integrated with orthogonal evidence.

## 2. What this workflow is, and is not

ORACLE is an **ENCODE-informed** Snakemake workflow for paired-end bulk
ATAC-seq. It uses commonly adopted tools and QC concepts, including Bowtie2,
duplicate and blacklist filtering, ATAC-specific Tn5 handling, FRiP, TSS
enrichment, library-complexity metrics, replicate-aware peak construction,
DESeq2, and MultiQC.

It is not an exact reimplementation of the official ENCODE pipeline and should
not be described as fully ENCODE-standard. The most important distinction is
peak reproducibility:

1. MACS3 calls peaks independently in each included biological replicate.
2. Peaks are converted to fixed-width windows around their summits.
3. A window must be supported by the configured number of replicates within at
   least one condition.
4. Supported windows are combined across conditions to form the counting
   universe.

This preserves condition-specific accessibility. It is transparent and useful
for differential counting, but it is **not IDR**. The workflow does not
calculate IDR scores, rescue ratios, self-consistency ratios, or pooled
pseudoreplicate diagnostics. Use the
[official ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline)
when exact ENCODE processing or submission requirements apply.

## 3. Design the experiment before sequencing

### Biological replication

The differential workflow requires at least two included biological
replicates per condition. Two replicates allow a model to run but provide weak
variance estimation and limited protection against an outlier. Three or more
independent biological replicates per condition are strongly preferred.

A biological replicate is an independently obtained biological unit. A second
sequencing lane from the same library is a technical replicate, not a new
biological replicate. Combine technical lanes for the same library before
listing that library once in the sample sheet, and retain the lane-to-library
mapping in the study metadata.

### Balance and confounding

Balance library preparation, sequencing lane, operator, processing date, sex,
age, donor, and other important covariates across conditions whenever
possible. A model cannot separate treatment from batch if every control was
prepared in one batch and every treatment sample in another.

Add measured covariates as sample-sheet columns and include justified terms in
the DESeq2 design, for example:

```yaml
diffacc:
  design: "~batch + condition"
  contrast: ["condition", "treatment", "control"]
```

For paired or repeated-donor designs, a formula such as `~donor + condition`
can be appropriate only when donor and condition are not confounded and the
design matrix is full rank. Decide the model before looking at differential
results.

### Depth and library preparation

Required depth depends on organism, cell number, signal quality, and the
question. A common planning range for human bulk ATAC-seq is tens of millions
of paired-end reads per library, with deeper sequencing for weak or complex
samples. Pilot data should be assessed by usable unique fragments and peak
saturation rather than raw reads alone.

Use matched read length and broadly similar sequencing depth across the study.
Record nuclei preparation, transposition conditions, PCR cycles, library kit,
sequencer, read length, and any deviations from the protocol.

### Controls

Standard ATAC-seq does not use an input-DNA control like ChIP-seq. Biological
controls, batch balance, and independent replication are still essential.
Include positive process controls when the experiment or sample type is new.

## 4. Prepare the inputs

### Local paired FASTQs

Each library needs an R1 and R2 FASTQ. Paths must end in `.fastq`, `.fq`,
`.fastq.gz`, or `.fq.gz`; gzip compression is recommended for large raw files.
Confirm that files exist and compressed files are intact:

```bash
test -s data/control_r1_R1.fastq.gz
test -s data/control_r1_R2.fastq.gz
gzip -t data/control_r1_R1.fastq.gz
gzip -t data/control_r1_R2.fastq.gz
```

FASTQ paths are interpreted from the repository root unless absolute paths are
used. Avoid spaces in paths and sample identifiers. Keep raw data read-only
and backed up outside the repository.

### SRA input

Use run accessions beginning with SRR, ERR, or DRR. A study accession, GEO
series, BioProject, or sample accession is not a sequencing run and cannot be
used directly. Verify that every selected run is paired-end and belongs to the
intended assay, organism, assembly context, condition, and biological
replicate.

The workflow downloads SRA runs with the SRA Toolkit. Production runs should
leave `sra.max_spots` at zero; a spot cap creates a subsample suitable only for
plumbing tests and will often be too shallow for peak calling.

### Create project copies from the tracked templates

`config/samples.tsv` and `config/config.yaml` are reusable tracked templates.
Do not enter project data or paths in them. Create the ignored working copies
once, then edit only those copies:

```bash
cp config/samples.tsv config/samples.project.tsv
cp config/project.example.yaml config/project_overrides.yaml
```

Set this value in `config/project_overrides.yaml`:

```yaml
samples: "config/samples.project.tsv"
```

The copied overlay already puts generated files under a `my_project`
namespace. Replace that token with a short name for this analysis before the
first dry run; do not reuse an earlier project's output or reference roots.

The sample sheet is tab-separated. Its first six columns are required;
`include` is optional, and extra design columns may follow:

| Column | Meaning |
|---|---|
| `sample` | Unique safe library identifier |
| `condition` | Safe biological-group identifier used for replicate support and usually differential testing |
| `replicate` | Replicate identifier unique within a condition |
| `fq1` | Local R1 FASTQ path, blank for SRA mode |
| `fq2` | Local R2 FASTQ path, blank for SRA mode |
| `sra` | SRR/ERR/DRR run, blank for local mode |
| `include` | Optional strict true/false QC gate; blank means included |

Both `sample` and `condition` must match the exact contract
`^[A-Za-z0-9][A-Za-z0-9._-]*$`: the first character is a letter or digit, and
every later character is a letter, digit, dot, underscore, or hyphen. Spaces,
slashes, and shell punctuation are rejected.

Example:

```text
sample	condition	replicate	fq1	fq2	sra	include	batch
control_r1	control	1	data/control_r1_R1.fastq.gz	data/control_r1_R2.fastq.gz		true	A
control_r2	control	2	data/control_r2_R1.fastq.gz	data/control_r2_R2.fastq.gz		true	B
treatment_r1	treatment	1	data/treatment_r1_R1.fastq.gz	data/treatment_r1_R2.fastq.gz		true	A
treatment_r2	treatment	2	data/treatment_r2_R1.fastq.gz	data/treatment_r2_R2.fastq.gz		true	B
```

Extra columns such as `batch`, `donor`, `sex`, or `timepoint` may be used in a
DESeq2 design. Keep them simple categorical values and avoid missing values in
columns used by the model.

Each row must provide exactly one input mode:

- both `fq1` and `fq2`, with `sra` blank; or
- one `sra` run, with both FASTQ columns blank.

The workflow validates identifiers, file existence, input mode, duplicate
sample names, replicate labels, contrast levels, replicate counts, and feature
dependencies before scheduling expensive work.

### Using the inclusion gate

Add `include` before a reviewed run when a library has a documented technical
failure. Preserve excluded rows rather than deleting them. Good reasons include
a failed library preparation, sample swap, severe contamination, or a
predefined QC failure. A weak biological effect is not a technical exclusion
criterion.

Keep an exclusion reason in the study metadata or an extra sample-sheet column.
Re-run the dry-run after changing inclusion because it changes replicate
support, the consensus universe, and the statistical design.

## 5. Install the workflow

### Requirements

Plan for:

- a 64-bit Linux server/workstation or WSL2;
- a reliable network connection for Conda packages and references;
- at least 20 GB free for software environments, plus references, FASTQs,
  intermediate BAMs, and results;
- enough RAM for alignment, FastQC, sorting, and R analysis;
- substantially more temporary and final storage than the compressed FASTQs.

The installer also supports macOS, but the continuous integration target is
Linux and some bioinformatics packages may vary by platform.

### Installation

From a clean checkout:

```bash
git clone https://github.com/challagandla/oracle-atac-seq.git
cd oracle-atac-seq
bash setup.sh
```

The setup script:

1. finds an existing Conda installation or installs a pinned, checksummed
   Miniforge release;
2. creates the `oracle-atac-runner` controller environment;
3. creates rule-specific Conda environments from `workflow/envs/*.yaml`;
4. checks the main executables and imports.

No environment activation is required. Check an installation without changing
it:

```bash
bash setup.sh --check
```

For a smaller initial setup, install only the runner and allow the first run to
create rule environments as needed:

```bash
bash setup.sh --runner-only
```

Do not manually merge all rule dependencies into one environment. Per-rule
environments reduce R, Python, Java, and native-library conflicts.

## 6. Choose and verify the reference genome

The reference FASTA, GTF, chromosome sizes, blacklist, annotation packages,
and motif background must use the same assembly and compatible chromosome
names. Assembly mismatch can remove valid reads, misplace annotations, or
produce plausible-looking but incorrect results.

`config/config.yaml` provides human, mouse, `mouse_mm10`, rat, and custom build
options. Review the resolved assembly rather than relying on the species name.
For local references, set explicit FASTA, GTF, and blacklist paths.
Provide the FASTA uncompressed (`.fa`); decompress an ordinary `.fa.gz` before
configuration so SAMtools and Bowtie2 can index it consistently.

Important blacklist rules:

- never use a blacklist from a different assembly by filename similarity;
- there is no official Boyle Lab mm39 blacklist;
- do not apply mm10 coordinates to mm39 without a documented, validated
  liftover and post-liftover checks;
- when blacklist removal is enabled, the resolved BED must contain at least
  one valid interval; an empty file is an error, not a no-op;
- for an assembly without a suitable blacklist, document that limitation
  instead of silently substituting another assembly.

For a custom genome, choose one complete reference route:

- set local, uncompressed FASTA and GTF paths; or
- leave those paths blank and provide the complete Ensembl species, release,
  and assembly tuple used to construct the download URLs.

Both routes require realistic effective and MACS genome sizes. Enrichment needs
a positive `taxid` plus `orgdb`; chromVAR needs that taxonomy identity plus a
matching `bsgenome` for sequence and motif selection. Annotation can use an
optional `txdb` or build one from the GTF. The named OrgDb, TxDb, and BSgenome
packages must exist in `workflow/envs/r.yaml` and describe the same assembly;
configuration strings do not install them. TSS extraction reconciles only
conservative chromosome aliases, rejects out-of-bounds genes, and requires at
least 95% of GTF gene TSSs to map inside the alignment FASTA; a failure here is
a reference mismatch to fix, not a QC threshold to weaken.

## 7. Configure a run

Edit `config/project_overrides.yaml` rather than the tracked base template or
rule code. The base configuration is loaded first and the project overlay
second, so values in the overlay take precedence. The most important sections
are below.

### Output locations

```yaml
results_dir: "results/my_project"
raw_dir: "results/my_project/fastq"
processed_dir: "results/my_project"
logs_dir: "logs/my_project"
reference_dir: "resources/reference_my_project"
```

Replace `my_project` once in all five paths. Separate raw-download, processed,
log, and reference roots prevent stale files from another analysis entering a
report or being mistaken for the current assembly. Keep the matching config
and sample sheet with that analysis.

### Filtering

The defaults retain high-quality proper pairs, remove mitochondrial alignments
and PCR duplicates, apply an assembly-matched blacklist when available, and
perform the Tn5 shift. Review rather than casually weakening these controls:

```yaml
filtering:
  min_mapq: 30
  remove_mito: true
  mitochondrial_contigs: ["chrM", "MT", "chrMT", "Mito", "M"]
  remove_blacklist: true
  keep_proper_pairs: true
  remove_duplicates: true
  tn5_shift: true
```

`keep_proper_pairs` and `tn5_shift` are required workflow invariants. The first
keeps peak calling, counting, FRiP, library-complexity, and usable-fragment QC
on one fragment universe; the second makes cut-site tracks and insertion-based
QC biologically meaningful. Preflight rejects `false` for either setting.

For a custom reference, replace or extend `mitochondrial_contigs` with the
exact mitochondrial sequence accession from that FASTA (for example,
`NC_012920.1`). The same list controls filtering and the reported mitochondrial
fraction, so those two stages cannot silently disagree. When removal is enabled,
the workflow fails if none of those names is present in the BAM header; this
protects against a silently wrong zero fraction. Set `remove_mito: false` only
when the reference genuinely has no mitochondrial sequence.

Very large custom references need Bowtie2's large index files. Set
`alignment.large_index: true` before the first run; the workflow then requests
and tracks the six `.bt2l` files instead of `.bt2` files.

Duplicate removal is appropriate for the default bulk workflow, but very low
input assays need special care because real biological fragments and PCR
duplicates can be difficult to distinguish.

### Peak calling and consensus

```yaml
peaks:
  run_genrich: true
  macs3_qvalue: 0.01
  consensus_min_replicates: 2
  consensus_peak_width: 500
```

MACS3 is run in paired-end BAMPE mode. Do not add single-end
`--shift`/`--extsize` options to `macs3_extra`; the input validator rejects
that incompatible combination.

`consensus_min_replicates` is evaluated within each condition, not across the
entire study. A value of two means that a condition-specific peak can survive
when supported by two replicates from that condition. It does not mean that a
peak must appear in both conditions. During preflight, every included condition
must contain at least this many distinct biological replicates; otherwise that
condition could never contribute a condition-specific peak and the run stops
before computation.

The fixed width prevents broad merged intervals from dominating counts and
makes intervals comparable, but the chosen width is an analysis parameter.
Inspect summit placement and peak architecture in a genome browser.

### Differential accessibility

```yaml
diffacc:
  enabled: true
  method: "DESeq2"
  design: "~condition"
  contrast: ["condition", "treatment", "control"]
  fdr: 0.05
  lfc_threshold: 0
```

The contrast is `[factor, numerator, denominator]`. Positive log2 fold change
means more accessible in the numerator. Every design variable must be a
sample-sheet column.

Significance, direction, and `lfc_threshold` use DESeq2's unshrunken
`log2FoldChange` together with `padj`. The ashr-shrunken value is used for
effect ranking and MA/volcano visualization, not to decide membership in the
up/down BED files. HOMER and functional enrichment inherit those raw-effect
up/down calls.

### Optional stages

```yaml
annotation:
  enabled: true

functional_enrichment:
  enabled: true
  kegg: false

motif:
  enabled: true

chromvar:
  enabled: false

footprinting:
  enabled: false
  motif_db: ""
```

Dependencies are enforced:

- motif enrichment requires differential up/down peak sets;
- functional enrichment runs GO analysis and requires differential results and
  annotation;
- KEGG is optional and opt-in because it calls a live service during analysis;
- TOBIAS requires a user-supplied motif database;
- non-human annotation and chromVAR require matching R genome packages.

Enable optional stages only when their biological assumptions and reference
requirements are satisfied.

## 8. Preview and run

### Preflight

Always inspect the planned DAG before computing:

```bash
bash run.sh --dry-run \
  --configfile config/config.yaml config/project_overrides.yaml
```

A dry-run checks configuration and input logic but does not prove that remote
files are downloadable, every tool can process the full data, or the libraries
have adequate biological quality.

### Outcome-blind QC pass

Build the dedicated library-review report before differential analysis:

```bash
bash run.sh --cores 8 qc_review \
  --configfile config/config.yaml config/project_overrides.yaml
```

Open `results/my_project/qc/qc_review_report.html` and follow Section 10. Its
exact manifest includes per-library reads, alignments, peaks, signal, FRiP,
complexity, TSS, genome-bin correlation/PCA, and replicate-level QC, but
excludes consensus counts, DESeq2,
enrichment, motif, footprint, and chromVAR results. Record technical decisions,
set any rejected row's `include` cell to `false`, freeze the sample sheet, and repeat
the dry-run.

Genome-bin correlation and PCA require at least two included libraries. A
single-library project can be used only for technical QC, not replication or
differential analysis. Put this complete overlay in the project configuration;
preflight intentionally rejects the default multi-library analysis settings:

```yaml
qc:
  sample_similarity: false
peaks:
  consensus_min_replicates: 1
diffacc:
  enabled: false
functional_enrichment:
  enabled: false
motif:
  enabled: false
footprinting:
  enabled: false
chromvar:
  enabled: false
```

The single-library report still provides raw-read, alignment, complexity,
fragment-size, peak, FRiP, TSS, and track-level QC. It cannot establish
reproducibility or support condition-level biological claims.

Categorical QC/report overlays support at most 12 conditions. For a larger
design, consolidate display groups or set both `qc.sample_similarity: false`
and `report.enabled: false`; the statistical tables remain available.

### Full run after the inclusion set is frozen

```bash
bash run.sh --cores 8 \
  --configfile config/config.yaml config/project_overrides.yaml
```

Choose cores based on available CPU, RAM, and storage throughput. More parallel
jobs can make performance worse when several FastQC, sorting, or alignment jobs
compete for memory and disk.

### Specific targets

```bash
# Rebuild only the outcome-blind QC review report.
bash run.sh --cores 8 qc_review \
  --configfile config/config.yaml config/project_overrides.yaml

# Stop after the consensus universe.
bash run.sh --cores 8 \
  results/my_project/consensus/consensus_peaks.bed \
  --configfile config/config.yaml config/project_overrides.yaml

# Build the differential results and their prerequisites.
bash run.sh --cores 8 \
  results/my_project/diffacc/differential_accessibility.tsv \
  --configfile config/config.yaml config/project_overrides.yaml
```

### Configuration overlays

Always load the tracked base first and the ignored project overlay second:

```bash
bash run.sh --dry-run \
  --configfile config/config.yaml config/project_overrides.yaml
bash run.sh --cores 8 \
  --configfile config/config.yaml config/project_overrides.yaml
```

Use the identical config-file order for dry-runs, production runs, restarts,
and summaries.

## 9. What happens at each stage

### Raw-read QC and trimming

FastQC reports base quality, adapter content, sequence duplication, and other
read-level features. fastp trims adapters and low-quality sequence. ATAC-seq
libraries often show non-random sequence content because Tn5 integration is
not random; interpret warnings in assay context.

### Alignment and filtering

Bowtie2 aligns paired reads while allowing long inserts that include
nucleosomal fragments. The workflow sorts alignments, keeps configured proper
pairs and mapping quality, removes mitochondrial alignments when enabled,
marks/removes duplicates, harmonizes blacklist chromosome names, filters
blacklisted regions, and indexes the final BAM.

Track attrition across these steps. A high raw mapping rate can still end in
few usable nuclear fragments because of duplicates, mitochondrial DNA,
multimapping, or improper pairs.

### Tn5 shift and coverage

Tn5 integration produces a strand-dependent offset between read alignment and
cut site. The shifted BAM supports cut-site-aware analyses. deepTools creates
CPM-normalized bigWig tracks for genome browsers.

The workflow writes a fragment-level CPM bigWig for browser visualization and
a one-base Tn5 cut-site CPM bigWig for cut-site-resolution QC and inspection.
CPM tracks are useful for visualization, not a substitute for count-based
DESeq2 normalization. Global chromatin changes can make equal-total scaling
misleading.

### Peak calling

MACS3 calls paired-end peaks independently for every sample. A library with no
peaks is considered a failure rather than allowed to generate empty downstream
tables. When enabled, Genrich pools the filtered replicates of each condition
in ATAC mode and produces a secondary condition-level peak set.

Genrich is a cross-check, not the source of the default differential count
matrix. Large disagreement between callers or replicates should trigger
inspection of library depth, backgrounds, chromosome naming, and browser
tracks.

### Consensus and counting

The workflow normalizes per-sample peaks to fixed-width, summit-centered
intervals and requires within-condition replicate support. It then uses
deterministic non-maximum suppression to choose one non-overlapping
representative at each supported locus across conditions. Rsubread's
featureCounts counts paired fragments in this shared universe. Each library is
counted in an independent process with isolated scratch space, then a validation step merges
the results. The merge refuses mismatched peak rows, sample labels, non-integer
counts, zero assigned libraries, or a disagreement between the matrix column
sum and the featureCounts assignment summary. This design makes failed samples
individually restartable and prevents a partial cohort matrix from publishing.

The count universe depends on included samples. Excluding a library or changing
the consensus parameters can change every downstream row, so comparisons
between runs must record those settings.

### Differential and tertiary analysis

DESeq2 models peak counts using the configured design and contrast. ChIPseeker
annotates peaks; clusterProfiler tests genes assigned to differential peaks;
HOMER tests sequence motif enrichment; chromVAR estimates motif-associated
deviations; TOBIAS performs bias-corrected footprint analysis.

Each later stage adds assumptions. Interpret results in layers rather than
treating agreement among dependent analyses as independent validation.

## 10. Quality control: make a decision before biology

Open `results/my_project/qc/qc_review_report.html` first, then inspect the
underlying files and tracks. Compare libraries within the same experiment;
absolute thresholds are guides, not universal pass/fail rules. After the final
run, `results/my_project/qc/multiqc_report.html` adds outcome figures for a
combined report, but it is not the inclusion-decision report.

Common bulk ATAC-seq guideposts include:

| Metric | Useful guide | Interpretation |
|---|---|---|
| Mapping and retained fraction | High and broadly balanced | Follow attrition from raw to proper pairs to final usable fragments |
| Usable fragments | `<1,000,000` fail; `>=1,000,000` and `<10,000,000` review; `>=10,000,000` no flag | Required depth depends on organism, cell type, signal, and question |
| Mitochondrial fraction | `>0.50` fail; `>0.30` and `<=0.50` review; `<=0.30` no flag | High values often indicate poor nuclei preparation or excess mitochondrial carryover |
| Per-sample peaks | `<1,000` fail; `>=1,000` and `<20,000` review; `>=20,000` no flag | Interpret with depth, genome size, FRiP, and browser signal |
| FRiP | `<0.20` fail; `>=0.20` and `<0.30` review; `>=0.30` no flag | Same-library MACS3 peaks make this an optimistic within-library QC metric; it is not directly comparable to independent or official ENCODE FRiP |
| TSS enrichment | `<5` fail; `>=5` and `<7` review; `>=7` no flag for the human preset | Implementation and gene annotation change the scale |
| NRF | `<0.50` fail; `>=0.50` and `<0.80` review; `>=0.80` no flag | Lower values indicate redundancy; values above 0.9 are commonly considered strong |
| PBC1 | Display only; no automated status boundary | Values above 0.9 are a useful complexity guide |
| PBC2 | Display only; no automated status boundary | Values above 3 are a useful complexity guide; interpret with depth and NRF/PBC1 |
| Fragment sizes | Nucleosome-free enrichment plus mono-/di-nucleosome periodicity | Missing periodicity can indicate poor nuclei or high background |
| Replicate similarity | Same-condition samples should be coherent | Inspect PCA, correlation, coverage, and peak overlap together |

Read-1 records represent fragments in these flow-through metrics. The
mitochondrial fraction is raw aligned proper-pair fragments with either mate on
a configured mitochondrial contig divided by all raw aligned proper-pair
fragments, before MAPQ and duplicate filtering. The retained fraction is final
usable proper-pair fragments divided by that same raw proper-pair denominator.
These definitions matter when comparing values from another pipeline.

These are approximate ENCODE-informed bands. The workflow places transparent
`pass`, `review`, or `fail` flags in
`results/my_project/qc/summary/qc_decisions_mqc.tsv` and MultiQC, but it never
excludes a library automatically. The inequalities above are exact: an equality belongs
to the band shown with `>=` or `<=`. PBC1 and PBC2 are reported for review but
do not change `QC_status`. TSS flagging is applied only to the human preset
because the scale is annotation- and organism-dependent. The workflow's TSS
score is not guaranteed to be numerically identical to the official ENCODE
pipeline. For a formal submission, calculate and report the official metrics
with the official implementation.

ORACLE calculates each library's FRiP against MACS3 peaks called from that same
library. This is useful for internal QC but benefits from peak-selection on the
data being scored. Compare it only within a consistently processed study, and
do not present it as an independent validation statistic or as numerically
equivalent to official ENCODE FRiP.

### A practical QC order

1. Confirm sample identity and expected read pairing.
2. Review FastQC and fastp for catastrophic read or adapter problems.
3. Compare total, mapped, proper-pair, mitochondrial, duplicate, and final
   usable-fragment counts.
4. Inspect NRF, PBC1, and PBC2 for bottlenecks.
5. Examine fragment periodicity and TSS enrichment.
6. Review FRiP and the number/shape of per-sample peaks.
7. Inspect bigWigs at known open loci and background regions.
8. Compare biological replicates in PCA, correlation, fingerprint, and peak
   support.
9. Decide and document exclusions without consulting desired differential
   outcomes.
10. Freeze the included sample set before final differential interpretation.

No single metric rescues a poor library. Conversely, one warning does not
automatically invalidate a library if the remaining evidence is coherent and
the reason is understood.

## 11. Understand the main outputs

Paths below use the copied overlay's `my_project` namespace. Substitute the
name chosen for your analysis.

### QC reports and unsupervised figures

- `results/my_project/qc/qc_review_report.html`: outcome-blind report used to
  freeze library inclusion.
- `results/my_project/qc/multiqc_report.html`: final combined report, including
  post-analysis figures.
- `results/my_project/qc/summary/qc_decisions_mqc.tsv`: auditable library
  metrics and review flags.
- `results/my_project/figures/qc_sample_pca.pdf` and
  `qc_sample_correlation.pdf`: vector, genome-bin sample-similarity figures
  computed without differential outcomes.

### Alignments and tracks

- `results/my_project/filtered/<sample>.filtered.bam`: final analysis-ready paired
  alignments.
- `results/my_project/shifted/<sample>.shifted.bam`: Tn5-shifted alignments used by
  cut-site-aware stages.
- `results/my_project/coverage/<sample>.cpm.bw`: CPM-normalized paired-fragment browser
  track.
- `results/my_project/coverage/<sample>.cutsites.cpm.bw`: one-base CPM-normalized Tn5
  insertion track.

Load bigWigs and peaks into IGV or the UCSC Genome Browser with the exact same
assembly. Check known promoters, expected cell-type loci, broad background,
blacklisted regions, and replicate concordance.

### Peaks and consensus

- `results/my_project/peaks/macs3/<sample>_peaks.narrowPeak`: per-sample MACS3 calls.
- `results/my_project/peaks/genrich/<condition>.narrowPeak`: optional pooled Genrich calls.
- `results/my_project/consensus/consensus_peaks.bed`: fixed-width,
  within-condition-supported peak universe.
- `results/my_project/consensus/consensus_peaks.saf`: featureCounts representation of that
  universe.

The consensus BED is suitable for quantification and many downstream tools,
but it is not an IDR peak set.

### Counts and differential results

- `results/my_project/counts/consensus_counts.tsv`: raw fragment counts.
- `results/my_project/counts/per_sample/<sample>.counts.tsv`: independently
  computed inputs retained for restartability and audit.
- `results/my_project/counts/consensus_counts.tsv.summary`: per-library
  assignment categories checked against the merged count matrix.
- `results/my_project/diffacc/differential_accessibility.tsv`: statistical results.
- `results/my_project/diffacc/normalized_counts.tsv`: DESeq2-normalized counts.
- `results/my_project/diffacc/tested_peaks.bed`: peaks with finite raw DESeq2 p-values;
  this is the opportunity universe used by enrichment and HOMER background.
- `results/my_project/diffacc/PCA_plot.pdf` and correlation/distance heatmaps: sample
  structure.
- `results/my_project/diffacc/MA_plot.pdf`, `volcano_plot.pdf`, and
  `differential_peaks_heatmap.pdf`: effect and significance summaries.

The differential table includes raw log2 fold change, shrunken effect size, raw
p-value, and adjusted p-value. Use:

- `padj` plus the absolute raw `log2FoldChange` threshold for significance and
  up/down membership;
- `lfcShrink` for stable effect-size ranking and visualization, not for the
  significance cutoff;
- normalized counts and browser tracks to check whether a result is driven by
  one sample;
- genomic coordinates as the primary peak identity.

An `NA` adjusted p-value commonly reflects low information or independent
filtering; it is not equivalent to a significant or unchanged peak.

### Annotation and enrichment

`results/my_project/annotation/consensus_peaks.annotated.tsv` reports genomic context and
nearby genes. A nearest gene is not necessarily the regulated gene,
particularly for distal enhancers.

GO enrichment, and optional KEGG enrichment when enabled, inherit the
peak-to-gene mapping, finite-p-value tested gene universe, matching OrgDb,
database version, and differential threshold. Mapping coverage is checked
against actual OrgDb Entrez keys before testing. Treat terms in
`results/my_project/enrichment/` as summaries of candidates, not proof of a pathway
mechanism.

## 12. Differential accessibility done carefully

Before interpreting p-values:

1. verify that numerator and denominator are in the intended order;
2. check that PCA and sample correlations agree with recorded metadata;
3. confirm that important covariates are modeled and not confounded;
4. inspect the count distribution and sample-level values for top peaks;
5. consider whether global accessibility shifts challenge the normalization
   assumptions;
6. report the number of tested peaks, FDR, effect-size rule, and replicate
   counts.

The default contrast:

```yaml
contrast: ["condition", "treatment", "control"]
```

means treatment versus control. Positive values are more accessible in the
treatment condition; negative values are more accessible in control.

Do not interpret an arbitrary p-value cutoff without effect size. Predefine a
biologically meaningful absolute log2-fold-change threshold when appropriate,
and show complete effect distributions. Very small studies can produce unstable
dispersion and effect estimates even when the software returns a table.

## 13. Motifs, chromVAR, and TOBIAS

### HOMER motif enrichment

HOMER compares sequence motifs in more-open and less-open differential peak
sets against the other finite-p-value DESeq2-tested peaks. This measured-locus
background reduces opportunity bias relative to arbitrary genomic sequence.
Review:

- the number and width of tested peaks;
- GC and sequence-composition matching;
- known and de novo motif results;
- related TFs sharing nearly identical DNA-binding motifs;
- whether the TF is expressed in the relevant cells.

A motif result supports a candidate TF family, not direct binding or causal
regulation.

### chromVAR

chromVAR calculates motif-associated accessibility deviations while using
GC/depth-matched background peaks. In this bulk workflow, samples rather than
single cells are the units. With few samples, motif variability rankings and
group patterns can be unstable; at least two included libraries are required.

Use chromVAR to prioritize motif families and compare patterns with DESeq2,
expression, and known biology. The output table preserves stable JASPAR motif
IDs as well as display names, and a seeded serial bootstrap backend makes the
default run reproducible. Do not treat a high deviation score as proof that one
specific family member is active.

### TOBIAS

TOBIAS corrects Tn5 sequence bias, calculates footprint scores, and runs
BINDetect with a user-provided motif database. Enable it only after checking:

- the motif database's assembly independence, provenance, version, and license;
- adequate depth and signal quality;
- the consensus peak universe;
- condition labels and the final filtered, unshifted BAMs that are merged within
  each condition for ATACorrect;
- sufficient compute and storage.

Footprint depth is affected by cleavage bias, local accessibility, TF residence
time, nucleosomes, and sequencing depth. A footprint is not equivalent to a
ChIP-seq binding event. Validate important candidates with orthogonal data or
experiments.

ORACLE merges included replicate BAMs within each condition before ATACorrect
and BINDetect. The resulting condition-pooled comparison is descriptive: it
does not model biological-replicate variance and must not be reported as
replicate-aware statistical evidence of differential TF binding.

JASPAR CORE data are available under CC BY 4.0 and require attribution. A
MEME-format file does not inherit the MEME Suite software license merely
because of its format; the actual motif database source and terms control.

## 14. Restart, recover, and update safely

Snakemake records completed outputs. After an interruption, use the same
command and configuration:

```bash
bash run.sh --cores 8 \
  --configfile config/config.yaml config/project_overrides.yaml
```

`run.sh` always enables Snakemake's incomplete-output recovery, so that normal
command automatically reschedules any job left incomplete by the interruption.

If the working directory is locked, first confirm that no workflow process is
still running. Only then:

```bash
bash run.sh --unlock \
  --configfile config/config.yaml config/project_overrides.yaml
```

Useful diagnostics:

```bash
# Preview what would run.
bash run.sh --dry-run \
  --configfile config/config.yaml config/project_overrides.yaml

# Summarize tracked output state.
bash run.sh --summary \
  --configfile config/config.yaml config/project_overrides.yaml

# Re-run one rule after a deliberate parameter/code review.
bash run.sh --cores 8 --forcerun macs3_callpeak \
  --configfile config/config.yaml config/project_overrides.yaml

# Allow more time for files on network storage to appear.
bash run.sh --cores 8 --latency-wait 60 \
  --configfile config/config.yaml config/project_overrides.yaml
```

Do not delete the complete results tree to solve one failed rule. Read the
rule's file under `logs/my_project/`, identify the cause, remove only a confirmed corrupt
output when necessary, and rerun. Keep the same config-file order on every
restart.

After pulling workflow updates, perform a new dry-run. Scientific logic,
environment definitions, or output formats may have changed; do not silently
mix products from different workflow commits in one analysis.

## 15. Provenance and reproducibility

The full default target writes the workflow-generated record under
`results/my_project/provenance/`:

- `effective_config.yaml`: the merged configuration used by Snakemake;
- `samples.tsv`: a copy of the configured sample sheet, including rows excluded
  by the `include` gate;
- `software_environments.sha256.tsv`: environment-specification hashes;
- `raw_inputs.sha256.tsv`: sample/mate, resolved path, byte count, and SHA-256
  for every selected FASTQ actually analyzed;
- `run_manifest.json`: workflow commit/dirty state, selected samples,
  conditions, ordered config sources, Snakemake version, raw-input identities,
  and reference checksums.

These files are a strong baseline, but they do not replace project-level
records. For every final analysis, also retain:

- the workflow Git commit;
- the exact config and sample sheet;
- the workflow-generated FASTQ checksums and any repository/archive accession;
- reference FASTA/GTF/blacklist source, version, assembly, and checksums;
- the dry-run and software-check output;
- the run command, start/end dates, and compute platform;
- inclusion/exclusion decisions and reasons;
- the MultiQC report, logs, and final result tables;
- versions and sources of motif and pathway databases.

Keep the additional project record under
`results/my_project/provenance/project_record/`. The following captures the tracked base,
ordered overlay, configured sample sheet, their hashes, the software check, and
a summary evaluated with the same base-then-overlay order as the run:

```bash
record="results/my_project/provenance/project_record"
mkdir -p "$record"
git rev-parse HEAD > "$record/workflow_commit.txt"
cp config/config.yaml config/project_overrides.yaml \
  config/samples.project.tsv "$record/"
if command -v sha256sum >/dev/null 2>&1; then
  sha256sum config/config.yaml config/project_overrides.yaml \
    config/samples.project.tsv > "$record/configuration.sha256"
else
  shasum -a 256 config/config.yaml config/project_overrides.yaml \
    config/samples.project.tsv > "$record/configuration.sha256"
fi
bash setup.sh --check > "$record/software_check.txt" 2>&1
bash run.sh --summary \
  --configfile config/config.yaml config/project_overrides.yaml \
  > "$record/output_summary.tsv"
```

Store this provenance with the analysis, not in a public workflow commit if it
contains local paths, accessions, sample metadata, or project details. Do not
commit FASTQs, BAMs, results, logs, credentials, protected metadata, or
study-specific analysis files to the public repository.

Conda YAML files improve reproducibility, but solver state, operating system,
remote reference content, and database releases can still change. For
long-lived regulated or collaborative work, archive resolved package exports
or containers and reference checksums in the project record.

## 16. Important limitations

- Paired-end bulk ATAC-seq only; no single-cell or single-end workflow.
- Fixed-width within-condition replicate support, not IDR.
- No pooled-pseudoreplicate rescue/self-consistency analysis.
- No spike-in normalization or explicit global-accessibility model.
- No automatic sample identity, sex, contamination, or genotype concordance
  check.
- No automatic technical-lane merging.
- Non-human annotation/chromVAR packages require manual verification.
- Some assemblies lack a validated official blacklist.
- Nearest-gene and pathway analyses do not establish regulatory targets.
- Motif enrichment, chromVAR, and TOBIAS do not prove TF binding.
- CPM bigWigs are visualization tracks, not differential measurements.
- Conda environments are versioned specifications, not immutable containers.
- Results require expert review and are not intended for clinical decisions.

## 17. Troubleshooting

### The dry-run reports invalid inputs

Read the complete error list. Common causes are missing columns, duplicated
sample IDs, only one mate, a non-run SRA accession, absent FASTQs, invalid
include values, missing contrast levels, or fewer than two included replicates
per condition.

### Reference or blacklist chromosomes do not match

Confirm the FASTA contig names, GTF seqnames, BAM header, peak files, and
blacklist. `chr1` and `1` can often be harmonized, but assembly mismatch cannot
be repaired by renaming chromosomes.

### MACS3 reports no peaks

Inspect usable fragment count, mapping/filtering attrition, TSS enrichment,
FRiP context, background tracks, genome size, assembly, and MACS3 log. Do not
lower the q-value or replicate requirement until the library quality and
reference are understood. A shallow test subsample may simply be incapable of
supporting peak calling.

### MultiQC looks incomplete or contains unexpected samples

Confirm that `results_dir`, `processed_dir`, and `logs_dir` belong to this run,
then inspect the listed samples and source paths. Use isolated output roots for
separate studies.

### R annotation or chromVAR fails for a non-human genome

Verify TxDb, OrgDb, BSgenome, assembly, and chromosome naming. Add reviewed
matching packages to the appropriate environment YAML, reinstall the affected
environment, dry-run, and document the change.

### TOBIAS will not start

`footprinting.enabled` requires an existing motif database path. Check its
format, permissions, source, version, and license. Then confirm that the final
filtered, unshifted BAMs used for condition merging, the FASTA, and consensus
peaks use the same assembly.

### A job is killed or the machine becomes unresponsive

Reduce `--cores`, check the rule log and system memory, and avoid running many
FastQC, alignment, or sorting tasks simultaneously. On shared systems, use a
reviewed scheduler profile with realistic per-rule resources.

## 18. Methods and citations

If results are published, cite this repository and the tools/stages actually
used. The environment/specification hashes identify what the workflow declared;
the final methods must also record resolved versions and cite every enabled
stage. Core references and project records are:

1. Buenrostro JD, Giresi PG, Zaba LC, Chang HY, Greenleaf WJ. Transposition of
   native chromatin for fast and sensitive epigenomic profiling of open
   chromatin, DNA-binding proteins and nucleosome position.
   *Nature Methods* (2013). [doi:10.1038/nmeth.2688](https://doi.org/10.1038/nmeth.2688)
2. Yan F, Powell DR, Curtis DJ, Wong NC. From reads to insight: a hitchhiker's
   guide to ATAC-seq data analysis. *Genome Biology* (2020).
   [doi:10.1186/s13059-020-1929-3](https://doi.org/10.1186/s13059-020-1929-3)
3. ENCODE Project. [ATAC-seq data standards](https://www.encodeproject.org/atac-seq/)
   and [official ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline).
4. Li Q, Brown JB, Huang H, Bickel PJ. Measuring reproducibility of high-throughput
   experiments. *Annals of Applied Statistics* (2011).
   [doi:10.1214/11-AOAS466](https://doi.org/10.1214/11-AOAS466)
5. Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ
   preprocessor. *Bioinformatics* (2018).
   [doi:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)
6. Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2.
   *Nature Methods* (2012).
   [doi:10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923)
7. Amemiya HM, Kundaje A, Boyle AP. The ENCODE Blacklist: identification of
   problematic regions of the genome. *Scientific Reports* (2019).
   [doi:10.1038/s41598-019-45839-z](https://doi.org/10.1038/s41598-019-45839-z)
8. Zhang Y et al. Model-based analysis of ChIP-Seq (MACS).
   *Genome Biology* (2008).
   [doi:10.1186/gb-2008-9-9-r137](https://doi.org/10.1186/gb-2008-9-9-r137)
9. Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose
   program for assigning sequence reads to genomic features.
   *Bioinformatics* (2014).
   [doi:10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656)
10. Love MI, Huber W, Anders S. Moderated estimation of fold change and
    dispersion for RNA-seq data with DESeq2. *Genome Biology* (2014).
    [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)
11. Yu G, Wang LG, He QY. ChIPseeker: an R/Bioconductor package for ChIP peak
    annotation, comparison and visualization. *Bioinformatics* (2015).
    [doi:10.1093/bioinformatics/btv145](https://doi.org/10.1093/bioinformatics/btv145)
12. Heinz S et al. Simple combinations of lineage-determining transcription
    factors prime cis-regulatory elements required for macrophage and B cell
    identities. *Molecular Cell* (2010).
    [doi:10.1016/j.molcel.2010.05.004](https://doi.org/10.1016/j.molcel.2010.05.004)
13. Schep AN et al. chromVAR: inferring transcription-factor-associated
    accessibility from single-cell epigenomic data. *Nature Methods* (2017).
    [doi:10.1038/nmeth.4401](https://doi.org/10.1038/nmeth.4401)
14. Bentsen M et al. ATAC-seq footprinting unravels kinetics of transcription
    factor binding during zygotic genome activation. *Nature Communications*
    (2020). [doi:10.1038/s41467-020-18035-1](https://doi.org/10.1038/s41467-020-18035-1)
15. Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis
    results for multiple tools and samples in a single report.
    *Bioinformatics* (2016).
    [doi:10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)
16. Mölder F et al. Sustainable data analysis with Snakemake.
    *F1000Research* (2021).
    [doi:10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2)
17. Andrews S. [FastQC: a quality-control tool for high-throughput sequence data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
18. Li H et al. The Sequence Alignment/Map format and SAMtools.
    *Bioinformatics* (2009).
    [doi:10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352)
19. Broad Institute. [Picard toolkit documentation](https://broadinstitute.github.io/picard/).
20. Ramírez F et al. deepTools2: a next generation web server for deep-sequencing
    data analysis. *Nucleic Acids Research* (2016).
    [doi:10.1093/nar/gkw257](https://doi.org/10.1093/nar/gkw257)
21. Stephens M. False discovery rates: a new deal. *Biostatistics* (2017).
    [doi:10.1093/biostatistics/kxw041](https://doi.org/10.1093/biostatistics/kxw041)
22. Wu T et al. clusterProfiler 4.0: a universal enrichment tool for interpreting
    omics data. *The Innovation* (2021).
    [doi:10.1016/j.xinn.2021.100141](https://doi.org/10.1016/j.xinn.2021.100141)
23. Gaspar JM. [Genrich ATAC-seq peak caller](https://github.com/jsh58/Genrich).
24. NCBI. [SRA Toolkit](https://github.com/ncbi/sra-tools) and
    [usage documentation](https://github.com/ncbi/sra-tools/wiki).
25. Fornes O et al. JASPAR 2020: update of the open-access database of
    transcription factor binding profiles. *Nucleic Acids Research* (2020).
    [doi:10.1093/nar/gkz1001](https://doi.org/10.1093/nar/gkz1001)
26. Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing
    genomic features. *Bioinformatics* (2010).
    [doi:10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)
