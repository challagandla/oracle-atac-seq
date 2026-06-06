#!/usr/bin/env Rscript
# =============================================================================
# Annotate consensus ATAC peaks with ChIPseeker.
# Uses a Bioconductor TxDb if available (by name), otherwise builds a TxDb from
# the supplied GTF on the fly. Emits an annotated TSV + feature-distribution PDF.
# =============================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(ChIPseeker)
  library(GenomicFeatures)
  library(GenomicRanges)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--bed"),
  make_option("--gtf"),
  make_option("--txdb", default = ""),
  make_option("--orgdb", default = ""),
  make_option("--tss-up", type = "integer", default = 3000, dest = "tss_up"),
  make_option("--tss-down", type = "integer", default = 3000, dest = "tss_down"),
  make_option("--out-tsv", dest = "out_tsv"),
  make_option("--out-plot", dest = "out_plot")
)))

peaks <- readPeakFile(opt$bed)

# Resolve TxDb: try the named Bioconductor package, else build from GTF.
txdb <- NULL
if (nzchar(opt$txdb) && requireNamespace(opt$txdb, quietly = TRUE)) {
  library(opt$txdb, character.only = TRUE)
  txdb <- get(opt$txdb)
  message("Using TxDb package: ", opt$txdb)
} else {
  message("TxDb package unavailable; building TxDb from GTF: ", opt$gtf)
  txdb <- suppressWarnings(makeTxDbFromGFF(opt$gtf, format = "gtf"))
}

anno <- annotatePeak(
  peaks, TxDb = txdb,
  tssRegion = c(-opt$tss_up, opt$tss_down),
  level = "gene", verbose = FALSE,
  annoDb = if (nzchar(opt$orgdb) &&
               requireNamespace(opt$orgdb, quietly = TRUE)) opt$orgdb else NULL
)

df <- as.data.frame(anno)
write.table(df, opt$out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

pdf(opt$out_plot, width = 8, height = 5)
print(plotAnnoBar(anno))
print(plotDistToTSS(anno))
dev.off()

cat(sprintf("Annotated %d peaks -> %s\n", length(peaks), opt$out_tsv))
