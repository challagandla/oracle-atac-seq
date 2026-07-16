# =============================================================================
# common.smk — sample-sheet parsing, genome presets, helper functions
# =============================================================================
import os
import re
import shlex
import sys
import pandas as pd

# -----------------------------------------------------------------------------
# Output directory layout
# -----------------------------------------------------------------------------
RESULTS = config.get("results_dir", "results")
RAW = config.get("raw_dir") or f"{RESULTS}/fastq"
PROCESSED = config.get("processed_dir") or RESULTS
LOGS = config.get("logs_dir", "logs")
REF = config.get("reference_dir", "resources/reference")

# -----------------------------------------------------------------------------
# Sample sheet
# -----------------------------------------------------------------------------
samples_all = (
    pd.read_csv(
        config["samples"],
        sep="\t",
        dtype=str,
        comment="#",
    )
)
samples_all = samples_all.fillna("")

# Validate the biological and configuration contract before Snakemake schedules
# any work.  The optional include column is the deliberate QC-review gate:
# excluded libraries remain documented in the sheet but do not enter the DAG.
sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
from validate_inputs import validate_config_and_samples  # noqa: E402

samples = validate_config_and_samples(config, samples_all)
samples = samples.set_index("sample", drop=False).sort_index()

SAMPLES = list(samples["sample"])
MITOCHONDRIAL_CONTIGS = tuple(config["filtering"]["mitochondrial_contigs"])


def mitochondrial_contig_args(option):
    """Return one safely quoted command-line option per configured contig."""
    return " ".join(
        f"{option} {shlex.quote(contig)}" for contig in MITOCHONDRIAL_CONTIGS
    )


def is_sra(sample):
    """A sample is SRA-mode if it has no local fq1 but does have an accession."""
    row = samples.loc[sample]
    return (row.get("fq1", "") == "") and (row.get("sra", "") != "")


def sra_accession(sample):
    return samples.loc[sample, "sra"]


def raw_fastqs(sample):
    """Return (R1, R2) FASTQ paths for a sample, whether local or downloaded."""
    if is_sra(sample):
        return (
            f"{RAW}/{sample}_R1.fastq.gz",
            f"{RAW}/{sample}_R2.fastq.gz",
        )
    row = samples.loc[sample]
    return (row["fq1"], row["fq2"])


SCRIPTS = "workflow/scripts"


def script_inputs(*names):
    """The scripts a rule runs, declared as inputs.

    Snakemake tracks code changes for the `script:` and `run:` directives, but a
    `shell:` line that happens to call `Rscript workflow/scripts/foo.R` is opaque
    to it: editing foo.R leaves every output it produced looking up to date, and
    the next run quietly keeps the stale figures. Listing the script here is what
    makes `--rerun-triggers mtime` notice.

    Pass every file the rule actually executes, including the ones the entry
    point sources or imports (atac_theme.R, palette.py).
    """
    return [f"{SCRIPTS}/{n}" for n in names]


def conditions():
    return sorted(set(samples["condition"]))


def samples_in_condition(cond):
    return list(samples.loc[samples["condition"] == cond, "sample"])


# -----------------------------------------------------------------------------
# Wildcard constraints — pin wildcards to known values so Snakemake never
# mis-splits sample names that contain underscores (e.g. "ctrl_rep1_R1").
# -----------------------------------------------------------------------------
wildcard_constraints:
    sample="|".join(re.escape(s) for s in SAMPLES) if SAMPLES else "$^",
    cond="|".join(re.escape(c) for c in sorted(set(samples["condition"]))) or "$^",
    read="R1|R2",
    direction="up|down",


# -----------------------------------------------------------------------------
# Genome presets
# -----------------------------------------------------------------------------
# effective_genome_size: deepTools mappable size for 2x ~50-100bp reads.
# macs_gsize: MACS3 -g shorthand or integer.
# ENCODE blacklist v2 (Amemiya et al. 2019) URLs.
GENOME_PRESETS = {
    "human": {
        "ensembl_species": "homo_sapiens",
        "ensembl_assembly": "GRCh38",
        "ensembl_release": 111,
        "taxid": 9606,
        "effective_genome_size": 2913022398,
        "macs_gsize": "hs",
        "blacklist_url": "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz",
        "blacklist_md5": "393688b4f06c9ce26165d47433dd8c37",
        "txdb": "TxDb.Hsapiens.UCSC.hg38.knownGene",
        "orgdb": "org.Hs.eg.db",
        "homer_genome": "hg38",
        "bsgenome": "BSgenome.Hsapiens.UCSC.hg38",
    },
    "mouse": {  # GRCm39 / mm39
        "ensembl_species": "mus_musculus",
        "ensembl_assembly": "GRCm39",
        "ensembl_release": 111,
        "taxid": 10090,
        "effective_genome_size": 2654621783,
        "macs_gsize": "mm",
        # Boyle Lab publishes mm10, not mm39. Chromosome renaming is not a
        # coordinate lift, so mm39 requires an assembly-matched user BED.
        "blacklist_url": "",
        "blacklist_md5": "",
        "txdb": "TxDb.Mmusculus.UCSC.mm39.refGene",
        "orgdb": "org.Mm.eg.db",
        "homer_genome": "mm39",
        "bsgenome": "BSgenome.Mmusculus.UCSC.mm39",
    },
    "mouse_mm10": {
        "ensembl_species": "mus_musculus",
        "ensembl_assembly": "GRCm38",
        "ensembl_release": 102,
        "taxid": 10090,
        "effective_genome_size": 2652783500,
        "macs_gsize": "mm",
        "blacklist_url": "https://www.encodeproject.org/files/ENCFF543DDX/@@download/ENCFF543DDX.bed.gz",
        "blacklist_md5": "4d2d98597e8301eddf7bc7b6818e2142",
        "txdb": "TxDb.Mmusculus.UCSC.mm10.knownGene",
        "orgdb": "org.Mm.eg.db",
        "homer_genome": "mm10",
        "bsgenome": "BSgenome.Mmusculus.UCSC.mm10",
    },
    "rat": {  # mRatBN7.2 / rn7
        "ensembl_species": "rattus_norvegicus",
        "ensembl_assembly": "mRatBN7.2",
        "ensembl_release": 111,
        "taxid": 10116,
        "effective_genome_size": 2626580772,
        "macs_gsize": "2.6e9",
        "blacklist_url": "",  # No official ENCODE blacklist for rn7; see README.
        "blacklist_md5": "",
        "txdb": "TxDb.Rnorvegicus.UCSC.rn7.refGene",
        "orgdb": "org.Rn.eg.db",
        "homer_genome": "rn7",
        "bsgenome": "BSgenome.Rnorvegicus.UCSC.rn7",
    },
}


def genome_cfg():
    build = config["genome"]["build"]
    if build == "custom":
        c = dict(config["genome"]["custom"])
        c["blacklist_url"] = ""
        return c
    if build not in GENOME_PRESETS:
        raise ValueError(
            f"Unknown genome.build '{build}'. "
            f"Choose one of {list(GENOME_PRESETS)} or 'custom'."
        )
    return GENOME_PRESETS[build]


GENOME = genome_cfg()


def macs_gsize():
    return GENOME["macs_gsize"]


def effective_genome_size():
    return GENOME["effective_genome_size"]


# Reference file targets (produced by refs.smk or supplied by the user)
def genome_fasta():
    return config["genome"]["fasta"] or f"{REF}/genome.fa"


def genome_gtf():
    return config["genome"]["gtf"] or f"{REF}/genes.gtf"


def blacklist_source():
    """The blacklist as supplied: user-provided path, or the downloaded ENCODE BED.

    This is the file *before* chromosome names are reconciled with the genome.
    Nothing downstream should read it directly.
    """
    if config["genome"]["blacklist"]:
        return config["genome"]["blacklist"]
    if GENOME.get("blacklist_url"):
        return f"{REF}/blacklist.raw.bed"
    return ""


def blacklist_bed():
    """The blacklist actually used for filtering: named the way the genome is.

    A user-supplied blacklist goes through the same harmonisation as a downloaded
    one -- supplying your own file is not a promise that it uses your genome's
    chromosome names.
    """
    if not config["filtering"].get("remove_blacklist", True):
        return ""
    if not blacklist_source():
        return ""   # no blacklist available -> filtering step is skipped
    # A distinct name: a user may legitimately point genome.blacklist at
    # {reference_dir}/blacklist.bed, and a rule whose input is its own output
    # is a cycle.
    return f"{REF}/blacklist.harmonized.bed"


def bowtie2_index_prefix():
    return f"{REF}/bowtie2/genome"


def bowtie2_index_files():
    suffix = "bt2l" if config.get("alignment", {}).get("large_index", False) else "bt2"
    parts = ["1", "2", "3", "4", "rev.1", "rev.2"]
    return [f"{bowtie2_index_prefix()}.{part}.{suffix}" for part in parts]


def bowtie2_build_mode():
    return "--large-index" if config.get("alignment", {}).get("large_index", False) else ""


# -----------------------------------------------------------------------------
# Final-output collector — respects the feature switches in config.yaml
# -----------------------------------------------------------------------------
def collect_final_outputs():
    out = []

    # Always: per-sample filtered BAMs, peaks, bigWigs, outcome-blind QC review,
    # and the final combined MultiQC report.
    out += expand(f"{PROCESSED}/filtered/{{s}}.filtered.bam", s=SAMPLES)
    out += expand(f"{PROCESSED}/peaks/macs3/{{s}}_peaks.narrowPeak", s=SAMPLES)
    out += expand(f"{PROCESSED}/coverage/{{s}}.cpm.bw", s=SAMPLES)
    out += expand(f"{PROCESSED}/coverage/{{s}}.cutsites.cpm.bw", s=SAMPLES)
    out += [f"{RESULTS}/qc/qc_review_report.html"]
    out += [f"{RESULTS}/qc/multiqc_report.html"]
    out += [f"{RESULTS}/consensus/consensus_peaks.bed"]
    out += [f"{RESULTS}/counts/consensus_counts.tsv"]

    if config["peaks"].get("run_genrich", False):
        out += expand(f"{PROCESSED}/peaks/genrich/{{c}}.narrowPeak", c=conditions())

    if config["diffacc"].get("enabled", False):
        out += [
            f"{RESULTS}/diffacc/differential_accessibility.tsv",
            f"{RESULTS}/diffacc/tested_peaks.bed",
        ]
        # Publication figure suite (volcano, MA, PCA, correlation/distance, DA
        # heatmap, scree) produced by the DESeq2 rule.
        out += expand(
            f"{RESULTS}/diffacc/{{f}}.pdf",
            f=["MA_plot", "volcano_plot", "PCA_plot", "scree_plot",
               "sample_correlation_heatmap", "sample_distance_heatmap",
               "differential_peaks_heatmap"],
        )

    if config["annotation"].get("enabled", False):
        out += [f"{RESULTS}/annotation/consensus_peaks.annotated.tsv"]

    if config.get("functional_enrichment", {}).get("enabled", False) and \
            config["annotation"].get("enabled", False) and \
            config["diffacc"].get("enabled", False):
        out += [
            f"{RESULTS}/enrichment/enrichment_dotplots.pdf",
            f"{RESULTS}/enrichment/enrichment_status.tsv",
            f"{RESULTS}/enrichment/tables",
        ]

    if config["motif"].get("enabled", False):
        out += [f"{RESULTS}/motif/up/homerResults.html"]
        out += [f"{RESULTS}/motif/down/homerResults.html"]

    if config["footprinting"].get("enabled", False):
        out += [f"{RESULTS}/footprint/bindetect/bindetect_results.txt"]
        out += [f"{RESULTS}/footprint/bindetect/bindetect_figures.pdf"]

    if config["chromvar"].get("enabled", False):
        out += [f"{RESULTS}/chromvar/chromvar_deviations.tsv"]
        out += [f"{RESULTS}/chromvar/chromvar_variability.pdf"]

    # ENCODE QC completeness + aggregate publication figures (report.smk).
    out += report_outputs()

    # Run record of what generated the results.
    out += [f"{RESULTS}/provenance/run_manifest.json"]
    out += [f"{RESULTS}/provenance/effective_config.yaml"]
    out += [f"{RESULTS}/provenance/samples.tsv"]
    out += [f"{RESULTS}/provenance/software_environments.sha256.tsv"]
    out += [f"{RESULTS}/provenance/raw_inputs.sha256.tsv"]

    return out
