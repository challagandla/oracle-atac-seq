# ATAC-seq Snakemake Pipeline — Full-Length (Bulk) Paired-End

[![CI](https://github.com/challagandla/oracle-atacseq/actions/workflows/ci.yml/badge.svg)](https://github.com/challagandla/oracle-atacseq/actions/workflows/ci.yml)

A reproducible, genome-agnostic Snakemake workflow for bulk ATAC-seq, built on
current best practices from the most-cited methods literature. It takes raw
paired-end FASTQs (or SRA accessions) all the way to differential chromatin
accessibility, peak annotation, transcription-factor motif enrichment, and
TF footprinting — with a single MultiQC report tying the QC together.

Everything runs on Ubuntu/Linux with **Snakemake + conda (Python & R)**. Each
step pins its own conda environment, so the only thing you install by hand is
Snakemake itself.

> **Scope.** This is for *bulk / full-length* ATAC-seq (the standard Buenrostro
> assay), not single-cell ATAC. For scATAC use ArchR/Signac or the `scvi-tools`
> PeakVI route.

---

## Table of contents

1. [What the pipeline does](#1-what-the-pipeline-does)
2. [Method choices and why](#2-method-choices-and-why-the-citations)
3. [Repository layout](#3-repository-layout)
4. [Installation](#4-installation)
5. [Quick start](#5-quick-start)
6. [Configuring your run](#6-configuring-your-run)
7. [Inputs: local FASTQs or SRA](#7-inputs-local-fastqs-or-sra)
8. [Genomes: human, mouse, rat, custom](#8-genomes-human-mouse-rat-custom)
9. [Running on a workstation or cluster](#9-running-on-a-workstation-or-cluster)
10. [Understanding the outputs](#10-understanding-the-outputs)
11. [Quality-control checklist](#11-quality-control-checklist-encode-thresholds)
12. [Step-by-step tutorial (worked example)](#12-step-by-step-tutorial-worked-example)
13. [Troubleshooting](#13-troubleshooting)
14. [References](#14-references)
15. [Downstream, tertiary analysis & recommendations](#15-downstream-tertiary-analysis--scientific-recommendations)

---

## 1. What the pipeline does

```
FASTQ (local or SRA)
   │  fastp            adapter/quality trimming (Nextera/Tn5 auto-detect)
   ▼
Bowtie2  -X 2000       paired-end alignment, keep nucleosomal fragments
   ▼
Filtering              MAPQ≥30 · proper pairs · drop chrM · MarkDuplicates ·
   │                   remove ENCODE blacklist regions
   ▼
Tn5 shift  (+4/−5)     base-pair-accurate cut sites (alignmentSieve --ATACshift)
   ├─────────────► deepTools bamCoverage ─► CPM bigWig (IGV/UCSC tracks)
   ▼
Peak calling           MACS3 (per sample)  +  Genrich -j (per condition)
   ▼
Consensus peaks        merged, reproducible peak set (≥N samples)
   ▼
featureCounts          peak × sample fragment-count matrix
   ▼
DESeq2                 differential accessibility (ashr-shrunken LFC)
   │                   + publication figures: volcano · MA · PCA + scree ·
   │                     sample correlation/distance · DA-peak heatmap
   ├─► ChIPseeker      genomic annotation + nearest gene
   ├─► clusterProfiler GO / KEGG enrichment of genes near up/down peaks
   ├─► HOMER           motif enrichment on up/down peaks
   ├─► TOBIAS          TF footprinting per condition (optional)
   └─► chromVAR        motif deviations + variability + TF heatmap (optional)

QC completeness        library complexity (NRF/PBC1/PBC2) · deepTools
   │                   fingerprint · fragment-size/nucleosome overlay ·
   │                   TSS metagene heatmap + numeric enrichment score
   ▼
MultiQC                one HTML report: FastQC, fastp, Bowtie2, Picard, fragment
                       sizes, TSS score, FRiP, complexity, fingerprint,
                       featureCounts + embedded DA/QC figures
```

## 2. Method choices and why (the citations)

These selections track the most-cited ATAC-seq methods papers and the ENCODE
ATAC-seq data standards. Full citations in [§14](#14-references).

| Step | Tool | Rationale |
|------|------|-----------|
| Trimming | **fastp** | Fast, auto-detects Nextera/Tn5 adapters for PE reads; single JSON for MultiQC (Chen et al. 2018). |
| Alignment | **Bowtie2** `-X 2000` | ENCODE-standard ATAC aligner; the large insert ceiling retains mono-/di-nucleosomal fragments (Langmead & Salzberg 2012). |
| Duplicate removal | **Picard MarkDuplicates** | ENCODE-standard; ATAC libraries are PCR-amplified and duplication-prone. |
| Blacklist | **ENCODE blacklist v2** | Removing hyper-signal artifact regions is *the* single most impactful QC step (Amemiya et al. 2019). |
| Tn5 shift | **deepTools `alignmentSieve --ATACshift`** | Canonical +4/−5 bp shift recentres reads on the Tn5 cut site (Buenrostro et al. 2013). |
| Peaks (primary) | **MACS3** `--nomodel --shift -75 --extsize 150` | De-facto ATAC peak caller; treats each read end as a cut and builds a 150 bp pileup (Zhang et al. 2008; MACS3). |
| Peaks (cross-check) | **Genrich** `-j` (ATAC mode) | Pools replicates, models the Tn5 cut, emits one reproducible peak list per group (Gaspar 2018). |
| Quantification | **featureCounts** | Fast fragment counting over the consensus peak set (Liao et al. 2014). |
| Differential | **DESeq2** | Robust negative-binomial model widely used for count-based accessibility (Love et al. 2014). |
| Annotation | **ChIPseeker** | Standard peak→feature/gene annotation in Bioconductor (Yu et al. 2015). |
| Functional enrichment | **clusterProfiler** | GO/KEGG over-representation of genes near differential peaks (Wu et al. 2021). |
| Motifs | **HOMER** | De-novo + known TF motif enrichment on differential peaks (Heinz et al. 2010). |
| Footprinting | **TOBIAS** | Bias-corrected, differential TF footprinting from bulk ATAC (Bentsen et al. 2020). |
| Motif deviations | **chromVAR** | Bias-corrected per-sample motif accessibility variability (Schep et al. 2017). |
| Library complexity | **NRF / PBC1 / PBC2** | ENCODE core ATAC QC for PCR bottlenecking / redundancy. |
| Signal:noise | **deepTools plotFingerprint** | Cumulative read-enrichment; separates focused from background libraries. |
| TSS enrichment | **deepTools computeMatrix** | Numeric TSS-enrichment score + metagene heatmap; the headline ATAC S/N metric. |
| Aggregate QC | **MultiQC** | One report across all tools (Ewels et al. 2016). |
| Workflow engine | **Snakemake** | Reproducible, conda-isolated, scalable DAG execution (Mölder et al. 2021). |

The overall design follows the widely-cited end-to-end guide *"From reads to
insight: a hitchhiker's guide to ATAC-seq data analysis"* (Yan et al., *Genome
Biology* 2020), and the ENCODE ATAC-seq processing standards.

## 3. Repository layout

```
atacseq-snakemake/
├── README.md                 ← this tutorial
├── environment.yaml          ← creates the `atacseq-smk` env (Snakemake)
├── config/
│   ├── config.yaml           ← MAIN knobs: genome, filtering, peaks, design
│   └── samples.tsv           ← your samples (one row per replicate)
├── profiles/default/         ← run profile (cores, conda, retries)
├── workflow/
│   ├── Snakefile             ← includes all rule modules + target rule
│   ├── rules/                ← one .smk per stage (refs, trim, align, …)
│   ├── envs/                 ← per-stage conda environments (pinned)
│   └── scripts/              ← Python + R helpers
└── .test/                    ← tiny fixtures for a DAG smoke-test
```

## 4. Installation

You need **conda** (Miniforge/Mambaforge recommended) on Linux/Ubuntu or WSL2.

```bash
# 1. Get the code
cd atacseq-snakemake

# 2. Create the controller environment (Snakemake itself)
conda env create -f environment.yaml
conda activate atacseq-smk

# 3. Everything else is installed automatically per-rule on first run,
#    because we pass --use-conda. No manual tool installs needed.
```

Two tool families need a **one-time external data install** (large genome
packages, so not auto-installed):

```bash
# HOMER genome (for motif enrichment) — pick your build:
#   after the motif env is created on first run, or in a manual homer env:
configureHomer.pl -install hg38      # or mm39 / mm10 / rn7

# Bioconductor TxDb/OrgDb/BSgenome (for ChIPseeker + chromVAR) — human example:
#   add these to workflow/envs/r.yaml or install into that env:
#   bioconductor-txdb.hsapiens.ucsc.hg38.knowngene
#   bioconductor-org.hs.eg.db
#   bioconductor-bsgenome.hsapiens.ucsc.hg38
```

If you skip those, the core pipeline (through peaks, consensus, counts, DESeq2)
still runs; only annotation/motif/footprint/chromVAR need them.

## 5. Quick start

```bash
conda activate atacseq-smk

# 1. Edit config/samples.tsv  → list your FASTQs (or SRA accessions)
# 2. Edit config/config.yaml  → set genome.build and the diffacc contrast
# 3. Preview the plan (no compute):
snakemake -n

# 4. Run it (16 cores, auto-build conda envs):
snakemake --use-conda --cores 16

# …or just the QC report, or just up to peaks:
snakemake --use-conda --cores 16 results/qc/multiqc_report.html
snakemake --use-conda --cores 16 results/peaks/macs3/ctrl_rep1_peaks.narrowPeak
```

## 6. Configuring your run

Everything lives in `config/config.yaml`. The most important knobs:

- **`genome.build`** — `human` | `mouse` | `mouse_mm10` | `rat` | `custom`.
  Leave `fasta`/`gtf`/`blacklist` blank to auto-download from Ensembl, or set
  local paths to skip downloads.
- **`filtering`** — MAPQ cutoff, mito removal, dedup, blacklist, Tn5 shift. The
  defaults are ENCODE-standard; you rarely need to change them.
- **`peaks`** — MACS3 q-value/flags, whether to also run Genrich, and
  `consensus_min_overlap` (how many samples must share a peak to keep it).
- **`diffacc`** — the model. Set `design`, and `contrast: [factor, A, B]` to
  test **A vs B** (log2FC > 0 means more open in A). The factor and any
  covariates must be **column names in `samples.tsv`**.
- **`motif` / `footprinting` / `chromvar` / `annotation`** — feature switches.
  `footprinting` and `chromvar` are **off by default** (slow / need extra data).

## 7. Inputs: local FASTQs or SRA

`config/samples.tsv` is tab-separated, one row per sample (per replicate):

```
sample      condition   replicate   fq1                         fq2                         sra
ctrl_rep1   control     1           data/fastq/c1_R1.fq.gz      data/fastq/c1_R2.fq.gz
ctrl_rep2   control     2           data/fastq/c2_R1.fq.gz      data/fastq/c2_R2.fq.gz
treat_rep1  treatment   1           data/fastq/t1_R1.fq.gz      data/fastq/t1_R2.fq.gz
treat_rep2  treatment   2           data/fastq/t2_R1.fq.gz      data/fastq/t2_R2.fq.gz
# Or pull from SRA (leave fq1/fq2 blank, fill sra):
gm12878_r1  control     1                                                                   SRR891268
```

Mode is auto-detected per row: filled `fq1` → local; only `sra` → download via
the SRA Toolkit. You can mix both in one sheet. Add extra columns (e.g.
`donor`, `batch`) and reference them in the DESeq2 `design` to model covariates.

## 8. Genomes: human, mouse, rat, custom

Built-in presets ship the correct effective genome size, MACS `gsize`, ENCODE
blacklist URL, and the R annotation package names:

| build | assembly | blacklist | notes |
|-------|----------|-----------|-------|
| `human` | GRCh38 / hg38 | ENCODE v2 hg38 | default |
| `mouse` | GRCm39 / mm39 | ENCODE v2 (mm10 lifted) | newest mouse |
| `mouse_mm10` | GRCm38 / mm10 | ENCODE v2 mm10 | legacy mm10 annotations |
| `rat` | mRatBN7.2 / rn7 | *(none official)* | see note below |
| `custom` | your own | your own | fill `genome.custom.*` |

**Rat note.** There is no official ENCODE blacklist for rn7. The pipeline
detects this and simply skips blacklist filtering (a warning is logged). If you
have a community blacklist BED, set `genome.blacklist` to its path and it will
be used. For motif/annotation in rat, install `org.Rn.eg.db`,
`TxDb.Rnorvegicus.UCSC.rn7.refGene`, and the matching BSgenome.

## 9. Running on a workstation or cluster

**Workstation** — use the bundled profile (cores/conda/retries preset):

```bash
snakemake --profile profiles/default
```

**Cluster (SLURM)** — copy the profile and add an executor. With Snakemake ≥8
use `snakemake-executor-plugin-slurm`; on 7.x use a generic cluster profile:

```bash
# 7.x example
snakemake --use-conda --jobs 100 \
  --cluster "sbatch -c {threads} --mem={resources.mem_mb} -t 04:00:00" \
  --default-resources mem_mb=16000
```

Per-rule threads come from `config.resources`; override any rule with
`--set-threads bowtie2_align=24`.

## 10. Understanding the outputs

Everything lands under `results/`:

```
results/
├── trimmed/                 fastp-trimmed FASTQs
├── aligned/                 raw sorted BAMs (+ .bai)
├── filtered/                *.filtered.bam  ← analysis-ready, deduped, no-blacklist
├── shifted/                 *.shifted.bam   ← Tn5-shifted (used for peaks/signal)
├── coverage/                *.cpm.bw        ← genome-browser tracks
├── peaks/
│   ├── macs3/               <sample>_peaks.narrowPeak  (per sample)
│   └── genrich/             <condition>.narrowPeak     (per group)
├── consensus/               consensus_peaks.bed / .saf
├── counts/                  consensus_counts.tsv (peak × sample matrix)
├── diffacc/                 differential_accessibility.tsv (+ ashr-shrunken LFC),
│                            diffacc_summary.tsv, normalized_counts.tsv,
│                            up_peaks.bed / down_peaks.bed, and figures:
│                            volcano_plot · MA_plot · PCA_plot + scree_plot ·
│                            sample_correlation_heatmap · sample_distance_heatmap ·
│                            differential_peaks_heatmap  (PDF + PNG)
├── annotation/              consensus_peaks.annotated.tsv, feature_distribution.pdf
├── enrichment/              enrichment_dotplots.pdf + GO/KEGG result tables
├── motif/{up,down}/         HOMER homerResults.html (+ knownResults)
├── footprint/               TOBIAS corrected tracks + BINDetect per condition
├── chromvar/                chromvar_deviations.tsv, chromvar_variability.tsv/.pdf,
│                            chromvar_deviation_heatmap.pdf
├── figures/                 fragment_size_distribution, tss_enrichment_heatmap,
│                            tss_enrichment_profile  (aggregate, publication-ready)
└── qc/
    ├── complexity/          NRF / PBC1 / PBC2 per sample
    ├── fingerprint/         deepTools fingerprint + JS-distance metrics
    ├── tss/                 per-sample + aggregate TSS metagene, numeric score
    └── multiqc_report.html  ← START HERE
```

The single most useful file is **`results/qc/multiqc_report.html`** — open it
first; it now embeds the headline DA figures (volcano, PCA, correlation) and the
QC figures (TSS metagene, fragment overlay) alongside the standard tool metrics.
The key biological result is
**`results/diffacc/differential_accessibility.tsv`** (columns: peak coords,
`log2FoldChange`, ashr-shrunken `lfcShrink`, `padj`, …), with the up/down BEDs
feeding HOMER and clusterProfiler. All vector PDFs under `results/figures/` and
`results/diffacc/` are drawn on a shared colour-blind-safe theme for direct use
in figures.

## 11. Quality-control checklist (ENCODE thresholds)

Check these in the MultiQC report before trusting downstream results:

- **Alignment rate** > 80–95 % to the nuclear genome.
- **Mitochondrial fraction** — high (>20–40 %) is common but wasteful; very
  high suggests poor nuclei prep. (Removed by the pipeline regardless.)
- **Fragment-size distribution** — must show the periodic ATAC pattern: a large
  sub-147 bp nucleosome-free peak, then mono-/di-nucleosome bumps (~200/~400 bp).
  A flat distribution = failed transposition.
- **TSS enrichment score** — a single number per library in the MultiQC
  "TSS enrichment score" table (and `results/qc/tss/tss_enrichment_mqc.tsv`).
  ENCODE flags **<5 (hg38) as poor, >7 as ideal**. The aggregate metagene
  heatmap/profile is in `results/figures/tss_enrichment_*`.
- **FRiP** (fraction of reads in peaks) — ENCODE recommends **> 0.2–0.3**.
  See `results/qc/frip/`.
- **Library complexity** — the MultiQC "Library complexity" table reports
  **NRF** (non-redundant fraction, ideal >0.9), **PBC1** (>0.9) and **PBC2**
  (>3). Low values flag PCR bottlenecking / over-amplified libraries — a
  complement to Picard's duplication rate.
- **Fingerprint** — `results/qc/fingerprint/` (and MultiQC): a strongly
  bowed cumulative curve = good signal concentration; a near-diagonal curve =
  little enrichment over background.
- **Replicate PCA / correlation** (`results/diffacc/PCA_plot.pdf`,
  `sample_correlation_heatmap.pdf`) — replicates should cluster by condition.

> **Worked QC example (this pipeline on Buenrostro GSE47753).** Running the
> included 13-library dataset makes the value of these metrics concrete: the
> 500-cell GM12878 and later CD4 time-points score FRiP ≈ 0.001, TSS enrichment
> < 2 and a flat fragment distribution, whereas the 50k-cell GM12878 and
> `cd4_day2_rep2` reach TSS enrichment 4–15 with clear nucleosome periodicity.
> **These metrics agree with each other** — use them together to decide which
> libraries to keep before interpreting differential results.

## 12. Step-by-step tutorial (worked example)

A complete public-data run (human GM12878 ATAC), end to end.

**Step 1 — Set up the sample sheet** to pull two replicates from SRA. Edit
`config/samples.tsv`:

```
sample        condition   replicate   fq1   fq2   sra
gm12878_rep1  control     1                       SRR891268
gm12878_rep2  control     2                       SRR891269
```

(For a *differential* test you need two conditions; here we just demo QC +
peaks. Add a `treatment` group to exercise DESeq2.)

**Step 2 — Point at the human genome.** In `config/config.yaml`:

```yaml
genome:
  build: "human"     # auto-downloads GRCh38 + hg38 blacklist on first run
```

**Step 3 — Preview the DAG** (no computation, catches config errors):

```bash
snakemake -n
```

You should see jobs for `sra_download`, `download_genome_fasta`,
`bowtie2_build`, `fastp`, `bowtie2_align`, the filtering chain,
`macs3_callpeak`, `consensus_peaks`, `featurecounts`, and `multiqc`.

**Step 4 — Run the QC-and-peaks core** on 16 cores:

```bash
snakemake --use-conda --cores 16 \
  results/qc/multiqc_report.html \
  results/consensus/consensus_peaks.bed
```

First run spends time downloading the genome, building the Bowtie2 index, and
creating conda envs; subsequent runs reuse all of it.

**Step 5 — Inspect QC.** Open `results/qc/multiqc_report.html`. Confirm the
fragment-size periodicity, TSS enrichment, and FRiP pass (§11). Load a
`results/coverage/*.cpm.bw` track in IGV next to a housekeeping gene — you
should see a sharp promoter peak.

**Step 6 — Differential accessibility.** With a two-condition sheet, set the
contrast in `config.yaml`:

```yaml
diffacc:
  enabled: true
  design: "~condition"
  contrast: ["condition", "treatment", "control"]
```

then:

```bash
snakemake --use-conda --cores 16 results/diffacc/differential_accessibility.tsv
```

Sort that TSV by `padj`; positive `log2FoldChange` = more open in *treatment*.

**Step 7 — Motifs in the differential peaks** (needs the HOMER genome
installed, §4):

```bash
snakemake --use-conda --cores 16 \
  results/motif/up/homerResults.html \
  results/motif/down/homerResults.html
```

Open the HOMER HTML to see enriched known + de-novo TF motifs driving the
opened/closed regions.

**Step 8 (optional) — TF footprinting** across conditions with TOBIAS. Provide
a motif database and switch it on:

```yaml
footprinting:
  enabled: true
  motif_db: "resources/motifs/JASPAR2022_CORE_vertebrates.meme"
```

```bash
snakemake --use-conda --cores 16
```

BINDetect output ranks TFs by differential binding between your conditions.

## 13. Troubleshooting

- **`conda: command not found` during `--use-conda`** — activate the controller
  env (`conda activate atacseq-smk`) and ensure conda/mamba is on `PATH`.
- **Genome download fails** — Ensembl FTP can rate-limit; set `genome.fasta`
  and `genome.gtf` to local copies, or re-run (rules retry).
- **`Genrich` errors on name sorting** — Genrich requires name-sorted BAMs; the
  pipeline handles this (`namesort_for_genrich`). If you feed your own BAMs,
  keep them name-sorted.
- **ChIPseeker can't find a TxDb** — it falls back to building one from your
  GTF automatically; for gene *symbols* install the matching `OrgDb`.
- **chromVAR/HOMER skipped** — they need the BSgenome / HOMER genome packages
  (§4). The rest of the pipeline is unaffected.
- **Few or no consensus peaks** — lower `peaks.consensus_min_overlap`, check
  FRiP/TSS QC; low-complexity libraries yield few reproducible peaks.
- **MACS3 wants to re-call peaks after upgrading** — the default `macs3_extra`
  no longer passes `-B` (it wrote multi-GB pileup bedGraphs nothing consumes;
  the narrowPeak output is unchanged). Snakemake sees the changed params and
  re-runs MACS3. To keep existing peaks on a resumed run, add
  `--rerun-triggers mtime`. Re-enable pileups by putting `-B` back in config.
- **Re-run after editing config** — `snakemake --use-conda --cores 16 -R $(snakemake --list-params-changes)` or simply target the affected outputs.

## 14. References

1. Buenrostro JD, Giresi PG, Zaba LC, Chang HY, Greenleaf WJ. *Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position.* **Nature Methods** 2013;10:1213–1218. (Original ATAC-seq; Tn5 +4/−5 shift.)
2. Yan F, Powell DR, Curtis DJ, Wong NC. *From reads to insight: a hitchhiker's guide to ATAC-seq data analysis.* **Genome Biology** 2020;21:22. (End-to-end best-practice guide this pipeline follows.) PMID 32014034.
3. ENCODE Project Consortium. *ATAC-seq Data Standards and Processing Pipeline.* ENCODE Portal (encodeproject.org). (QC thresholds, filtering standards.)
4. Langmead B, Salzberg SL. *Fast gapped-read alignment with Bowtie 2.* **Nature Methods** 2012;9:357–359.
5. Chen S, Zhou Y, Chen Y, Gu J. *fastp: an ultra-fast all-in-one FASTQ preprocessor.* **Bioinformatics** 2018;34:i884–i890.
6. Amemiya HM, Kundaje A, Boyle AP. *The ENCODE Blacklist: Identification of Problematic Regions of the Genome.* **Scientific Reports** 2019;9:9354.
7. Zhang Y, Liu T, Meyer CA, et al. *Model-based Analysis of ChIP-Seq (MACS).* **Genome Biology** 2008;9:R137. (MACS; MACS3 is the current release.)
8. Gaspar JM. *Genrich: detecting sites of genomic enrichment.* 2018. github.com/jsh58/Genrich. (ATAC mode `-j`.)
9. Liao Y, Smyth GK, Shi W. *featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features.* **Bioinformatics** 2014;30:923–930.
10. Love MI, Huber W, Anders S. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* **Genome Biology** 2014;15:550.
11. Yu G, Wang LG, He QY. *ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization.* **Bioinformatics** 2015;31:2382–2383.
12. Heinz S, Benner C, Spann N, et al. *Simple combinations of lineage-determining transcription factors prime cis-regulatory elements (HOMER).* **Molecular Cell** 2010;38:576–589.
13. Bentsen M, Goymann P, Schultheis H, et al. *ATAC-seq footprinting unravels kinetics of transcription factor binding during zygotic genome activation (TOBIAS).* **Nature Communications** 2020;11:4267. PMID 32848148.
14. Schep AN, Wu B, Buenrostro JD, Greenleaf WJ. *chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data.* **Nature Methods** 2017;14:975–978.
15. Ramírez F, Ryan DP, Grüning B, et al. *deepTools2: a next generation web server for deep-sequencing data analysis.* **Nucleic Acids Research** 2016;44:W160–W165. (bamCoverage, alignmentSieve `--ATACshift`, computeMatrix.)
16. Ewels P, Magnusson M, Lundin S, Käller M. *MultiQC: summarize analysis results for multiple tools and samples in a single report.* **Bioinformatics** 2016;32:3047–3048.
17. Mölder F, Jablonski KP, Letcher B, et al. *Sustainable data analysis with Snakemake.* **F1000Research** 2021;10:33.
18. Li H, Handsaker B, Wysoker A, et al. *The Sequence Alignment/Map format and SAMtools.* **Bioinformatics** 2009;25:2078–2079.
19. Wu T, Hu E, Xu S, et al. *clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.* **The Innovation** 2021;2:100141. (GO/KEGG over-representation.)

## 15. Downstream, tertiary analysis & scientific recommendations

Everything below runs from the same config and the same real data
(`ATAC-PRJNA207663_GSE47753`). Each stage is a switch in `config.yaml`.

**What now ships end-to-end**

| Layer | Output | Enable |
|-------|--------|--------|
| Differential figures | volcano, MA (ashr-shrunken), PCA + scree, sample correlation & distance heatmaps, top-DA-peak z-score heatmap | `diffacc.enabled` |
| Functional enrichment | clusterProfiler GO (BP/MF) + KEGG dot-plots + tables of genes near up/down peaks | `functional_enrichment.enabled` (needs `annotation.enabled`) |
| Motif deviations | chromVAR variability plot + top-variable-TF deviation heatmap | `chromvar.enabled` |
| ENCODE QC | NRF/PBC1/PBC2, deepTools fingerprint, fragment/nucleosome overlay, TSS metagene heatmap + numeric score | `qc.library_complexity` / `qc.fingerprint` / `report.enabled` |

All figures use one colour-blind-safe theme (`workflow/scripts/atac_theme.R`),
export **vector PDF + PNG**, and the headline ones are folded into MultiQC.

**Scientific recommendations (be critical of your own data)**

1. **Gate samples on QC before differential testing.** Combine FRiP (>0.2),
   TSS enrichment (>5–7), NRF/PBC and a clean nucleosome ladder — they agree.
   In the bundled dataset the 500-cell and late CD4 libraries fail all four and
   should be dropped or interpreted with heavy caution.
2. **Model known batches.** If `samples.tsv` carries a batch/donor column, set
   `diffacc.design: "~batch + condition"` — the DESeq2 script honours any
   formula whose terms exist in the sheet. Shrunken LFCs (`lfcShrink`, ashr) are
   used for ranking/MA/volcano so weakly-supported peaks don't dominate.
3. **Consensus strategy.** The default keeps peaks seen in ≥2 samples across the
   whole cohort. For many heterogeneous groups, consider building consensus
   *within condition* first (or an IDR step) to avoid one large group dominating.
4. **Signal normalisation.** bigWigs are CPM; for cross-sample browser
   comparisons `--normalizeUsing RPGC` (1× genome coverage) is often preferable —
   change it in `coverage.smk`.
5. **Genrich vs MACS3.** Both are produced; treat Genrich (per-condition,
   replicate-aware) as the reproducible cross-check on the per-sample MACS3 set.
6. **Rat (rn7)** has no ENCODE blacklist — supply one via `genome.blacklist` or
   expect a few artefact peaks.

**Natural extensions** (not yet wired): per-locus genome-browser panels
(pyGenomeTracks) for top DA peaks, IDR-based reproducible peaks, and
motif→footprint integration (enable `footprinting` with a JASPAR/MEME db).

---

*Built for reproducible bulk ATAC-seq on Ubuntu/Linux with Snakemake, conda,
Python and R. Verified by Snakemake DAG dry-run (95-job full workflow, all
stages enabled) and by running the new figure/QC scripts on real Buenrostro
GSE47753 data.*

---

## License & usage

The pipeline's own code is **MIT** (see [LICENSE](LICENSE)). It bundles no third-party code or data;
tools are conda-installed and invoked, so the MIT license is unaffected by the (incl. GPL) licenses of
those tools. Full breakdown: [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md).

> ⚠️ **The HOMER motif step is academic / non-profit only** — HOMER is freeware, not open-source, and
> not redistributable; commercial use requires the author's permission (and its genome packages are
> UCSC-derived). For commercial use, obtain permission or skip the HOMER motif step; peak calling,
> differential accessibility, and JASPAR footprinting do not depend on it.
