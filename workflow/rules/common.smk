# =============================================================================
# common.smk — sample-sheet parsing, genome presets, helper functions
# =============================================================================
import os
import pandas as pd

# -----------------------------------------------------------------------------
# Output directory layout
# -----------------------------------------------------------------------------
RESULTS = "results"
LOGS = "logs"
REF = "resources/reference"

# -----------------------------------------------------------------------------
# Sample sheet
# -----------------------------------------------------------------------------
samples = (
    pd.read_csv(
        config["samples"],
        sep="\t",
        dtype=str,
        comment="#",
    )
    .set_index("sample", drop=False)
    .sort_index()
)
samples = samples.fillna("")

SAMPLES = list(samples["sample"])


def is_sra(sample):
    """A sample is SRA-mode if it has no local fq1 but does have an accession."""
    row = samples.loc[sample]
    return (row.get("fq1", "") == "") and (row.get("sra", "") != "")


def sra_accession(sample):
    return samples.loc[sample, "sra"]


def raw_fastqs(sample):
    """Return (R1, R2) FASTQ paths for a sample, whether local or downloaded."""
    if is_sra(sample):
        acc = sra_accession(sample)
        return (
            f"{RESULTS}/fastq/{sample}_R1.fastq.gz",
            f"{RESULTS}/fastq/{sample}_R2.fastq.gz",
        )
    row = samples.loc[sample]
    return (row["fq1"], row["fq2"])


def conditions():
    return sorted(set(samples["condition"]))


def samples_in_condition(cond):
    return list(samples.loc[samples["condition"] == cond, "sample"])


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
        "effective_genome_size": 2913022398,
        "macs_gsize": "hs",
        "blacklist_url": "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz",
        "txdb": "TxDb.Hsapiens.UCSC.hg38.knownGene",
        "orgdb": "org.Hs.eg.db",
        "homer_genome": "hg38",
        "bsgenome": "BSgenome.Hsapiens.UCSC.hg38",
    },
    "mouse": {  # GRCm39 / mm39
        "ensembl_species": "mus_musculus",
        "ensembl_assembly": "GRCm39",
        "ensembl_release": 111,
        "effective_genome_size": 2654621783,
        "macs_gsize": "mm",
        "blacklist_url": "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz",
        "txdb": "TxDb.Mmusculus.UCSC.mm39.refGene",
        "orgdb": "org.Mm.eg.db",
        "homer_genome": "mm39",
        "bsgenome": "BSgenome.Mmusculus.UCSC.mm39",
    },
    "mouse_mm10": {
        "ensembl_species": "mus_musculus",
        "ensembl_assembly": "GRCm38",
        "ensembl_release": 102,
        "effective_genome_size": 2652783500,
        "macs_gsize": "mm",
        "blacklist_url": "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz",
        "txdb": "TxDb.Mmusculus.UCSC.mm10.knownGene",
        "orgdb": "org.Mm.eg.db",
        "homer_genome": "mm10",
        "bsgenome": "BSgenome.Mmusculus.UCSC.mm10",
    },
    "rat": {  # mRatBN7.2 / rn7
        "ensembl_species": "rattus_norvegicus",
        "ensembl_assembly": "mRatBN7.2",
        "ensembl_release": 111,
        "effective_genome_size": 2626580772,
        "macs_gsize": "2.6e9",
        "blacklist_url": "",  # No official ENCODE blacklist for rn7; see README.
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


def blacklist_bed():
    if config["genome"]["blacklist"]:
        return config["genome"]["blacklist"]
    if GENOME.get("blacklist_url"):
        return f"{REF}/blacklist.bed"
    return ""   # no blacklist available -> filtering step is skipped


def bowtie2_index_prefix():
    return f"{REF}/bowtie2/genome"


# -----------------------------------------------------------------------------
# Final-output collector — respects the feature switches in config.yaml
# -----------------------------------------------------------------------------
def collect_final_outputs():
    out = []

    # Always: per-sample filtered BAMs, peaks, bigWigs, and the MultiQC report.
    out += expand(f"{RESULTS}/filtered/{{s}}.filtered.bam", s=SAMPLES)
    out += expand(f"{RESULTS}/peaks/macs3/{{s}}_peaks.narrowPeak", s=SAMPLES)
    out += expand(f"{RESULTS}/coverage/{{s}}.cpm.bw", s=SAMPLES)
    out += [f"{RESULTS}/qc/multiqc_report.html"]
    out += [f"{RESULTS}/consensus/consensus_peaks.bed"]
    out += [f"{RESULTS}/counts/consensus_counts.tsv"]

    if config["peaks"].get("run_genrich", False):
        out += expand(f"{RESULTS}/peaks/genrich/{{c}}.narrowPeak", c=conditions())

    if config["diffacc"].get("enabled", False):
        out += [f"{RESULTS}/diffacc/differential_accessibility.tsv"]
        out += [f"{RESULTS}/diffacc/MA_plot.pdf"]

    if config["annotation"].get("enabled", False):
        out += [f"{RESULTS}/annotation/consensus_peaks.annotated.tsv"]

    if config["motif"].get("enabled", False):
        out += [f"{RESULTS}/motif/up/homerResults.html"]
        out += [f"{RESULTS}/motif/down/homerResults.html"]

    if config["footprinting"].get("enabled", False):
        out += expand(f"{RESULTS}/footprint/{{c}}_footprints.bed", c=conditions())

    if config["chromvar"].get("enabled", False):
        out += [f"{RESULTS}/chromvar/chromvar_deviations.tsv"]

    return out
