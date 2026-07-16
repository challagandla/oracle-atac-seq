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

# Harmonise chromosome naming without a species-specific "standard chromosome"
# filter. Try each TxDb style on a copy and retain the representation that maps
# the largest number of peaks. This preserves arbitrary custom contigs while
# still converting common Ensembl/UCSC differences such as 1 versus chr1.
n_before <- length(peaks)
tx_levels <- GenomeInfoDb::seqlevels(txdb)
shared_peak_count <- function(x) {
  sum(as.character(GenomicRanges::seqnames(x)) %in% tx_levels)
}
best <- peaks
best_mode <- "direct chromosome names"
best_n <- shared_peak_count(best)
target_styles <- tryCatch(
  GenomeInfoDb::seqlevelsStyle(txdb),
  error = function(e) character()
)
for (style in unique(target_styles)) {
  candidate <- tryCatch({
    x <- peaks
    suppressWarnings(GenomeInfoDb::seqlevelsStyle(x) <- style)
    x
  }, error = function(e) NULL)
  if (!is.null(candidate)) {
    n_shared <- shared_peak_count(candidate)
    if (n_shared > best_n) {
      best <- candidate
      best_n <- n_shared
      best_mode <- paste0("seqlevelsStyle=", style)
    }
  }
}
shared <- intersect(GenomeInfoDb::seqlevels(best), tx_levels)
if (best_n == 0L || length(shared) == 0L) {
  stop(sprintf(paste0(
    "peaks and the TxDb share no chromosome names after direct matching or ",
    "style conversion. peaks: %s ... TxDb: %s ..."),
    paste(utils::head(GenomeInfoDb::seqlevels(peaks), 3), collapse = ", "),
    paste(utils::head(tx_levels, 3), collapse = ", ")))
}
peaks <- GenomeInfoDb::keepSeqlevels(best, shared, pruning.mode = "coarse")
if (length(peaks) == 0L) {
  stop("no peaks remain on TxDb-annotated contigs")
}
retained_fraction <- length(peaks) / n_before
if (retained_fraction < 0.5) {
  stop(sprintf(
    paste0("annotation retained only %d/%d peaks (%.1f%%) on TxDb contigs; ",
           "refusing a severely biased subset"),
    length(peaks), n_before, 100 * retained_fraction
  ))
}
if (retained_fraction < 0.8) {
  warning(sprintf("annotation retained %d/%d peaks (%.1f%%) on TxDb contigs",
                  length(peaks), n_before, 100 * retained_fraction))
}
message(sprintf(
  "chromosome naming: %s; %d shared seqlevels; retained %d/%d peaks (%.1f%%)",
  best_mode, length(shared), length(peaks), n_before, 100 * retained_fraction
))

anno <- annotatePeak(
  peaks, TxDb = txdb,
  tssRegion = c(-opt$tss_up, opt$tss_down),
  level = "gene", verbose = FALSE,
  annoDb = if (nzchar(opt$orgdb) &&
               requireNamespace(opt$orgdb, quietly = TRUE)) opt$orgdb else NULL
)

df <- as.data.frame(anno)

# A run where nothing landed near a gene is a naming or annotation failure, not
# a biological finding: ATAC peaks are strongly promoter-enriched.
if (nrow(df) == 0L) {
  stop("ChIPseeker annotated 0 peaks; check the TxDb and the peak file.")
}
prom <- sum(grepl("^Promoter", df$annotation)) / nrow(df)
message(sprintf("annotated %d/%d peaks; %.1f%% in promoters", nrow(df), n_before, 100 * prom))
if (prom == 0) {
  stop(paste0("no peak was assigned to a promoter. ATAC peaks are promoter-enriched;\n",
              "  0% indicates the peaks and the annotation do not describe the same genome."))
}

write.table(df, opt$out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

pdf(opt$out_plot, width = 8, height = 5)
print(plotAnnoBar(anno))
print(plotDistToTSS(anno))
plotAnnoPie(anno)                     # base graphics: print() would be a no-op
try(print(upsetplot(anno)), silent = TRUE)   # cosmetic; absent in older ChIPseeker
dev.off()

cat(sprintf("Annotated %d peaks -> %s\n", nrow(df), opt$out_tsv))
