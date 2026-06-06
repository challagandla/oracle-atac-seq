#!/usr/bin/env Rscript
# =============================================================================
# chromVAR motif accessibility deviations from the consensus peak count matrix.
# Computes bias-corrected per-sample deviations and ranks motifs by variability
# (Schep et al. 2017, Nat Methods). Requires a BSgenome matching the genome and
# JASPAR2020 motifs.
# =============================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(chromVAR)
  library(motifmatchr)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(JASPAR2020)
  library(TFBSTools)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts"),
  make_option("--bed"),
  make_option("--samples"),
  make_option("--bsgenome", default = ""),
  make_option("--out-dev", dest = "out_dev"),
  make_option("--out-var", dest = "out_var")
)))

if (!nzchar(opt$bsgenome) || !requireNamespace(opt$bsgenome, quietly = TRUE)) {
  stop("chromVAR requires the BSgenome package '", opt$bsgenome,
       "'. Install it (see envs/r.yaml) and re-run.")
}
library(opt$bsgenome, character.only = TRUE)
genome <- get(opt$bsgenome)

# ---- peaks & counts ---------------------------------------------------------
bed <- read.delim(opt$bed, header = FALSE)
peaks <- GRanges(bed[[1]], IRanges(bed[[2]] + 1, bed[[3]]), peak = bed[[4]])

fc <- read.delim(opt$counts, comment.char = "#", check.names = FALSE)
counts <- as.matrix(fc[, 7:ncol(fc), drop = FALSE])
colnames(counts) <- sub("\\.filtered\\.bam$", "", basename(colnames(counts)))
rownames(counts) <- fc$Geneid

se <- SummarizedExperiment(assays = list(counts = counts), rowRanges = peaks)
se <- addGCBias(se, genome = genome)
se <- filterPeaks(se, non_overlapping = TRUE)

# ---- motifs -----------------------------------------------------------------
motifs <- getMatrixSet(JASPAR2020,
                       list(species = 9606, collection = "CORE"))
match <- matchMotifs(motifs, se, genome = genome)

dev <- computeDeviations(object = se, annotations = match)
variability <- computeVariability(dev)

write.table(data.frame(motif = rownames(deviations(dev)),
                       deviations(dev), check.names = FALSE),
            opt$out_dev, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(variability, opt$out_var, sep = "\t", quote = FALSE,
            row.names = FALSE)
cat("chromVAR deviations written.\n")
