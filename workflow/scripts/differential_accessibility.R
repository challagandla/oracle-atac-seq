#!/usr/bin/env Rscript
# =============================================================================
# Differential accessibility with DESeq2 over consensus ATAC peaks.
# Reads a featureCounts matrix + sample sheet, fits the configured design,
# extracts the requested contrast, and writes results, MA/PCA plots, and
# up/down peak BEDs for downstream motif analysis.
# =============================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(ggplot2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts"),
  make_option("--samples"),
  make_option("--design", default = "~condition"),
  make_option("--factor", default = "condition"),
  make_option("--numerator"),
  make_option("--denominator"),
  make_option("--fdr", type = "double", default = 0.05),
  make_option("--lfc", type = "double", default = 0),
  make_option("--outdir")
)))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- featureCounts matrix ---------------------------------------------------
fc <- read.delim(opt$counts, comment.char = "#", check.names = FALSE)
# Columns: GeneID Chr Start End Strand Length <bam1> <bam2> ...
peak_info <- fc[, c("Geneid", "Chr", "Start", "End")]
colnames(peak_info)[1] <- "peak"
count_cols <- fc[, 7:ncol(fc), drop = FALSE]
# featureCounts names columns by BAM path; reduce to sample name.
clean <- basename(colnames(count_cols))
clean <- sub("\\.filtered\\.bam$", "", clean)
colnames(count_cols) <- clean
mat <- as.matrix(count_cols)
rownames(mat) <- fc$Geneid

# ---- sample sheet -----------------------------------------------------------
coldata <- read.delim(opt$samples, comment.char = "#", stringsAsFactors = FALSE)
rownames(coldata) <- coldata$sample
coldata <- coldata[colnames(mat), , drop = FALSE]
coldata[[opt$factor]] <- factor(coldata[[opt$factor]])
# Make the denominator the reference level.
coldata[[opt$factor]] <- relevel(coldata[[opt$factor]], ref = opt$denominator)

stopifnot(all(colnames(mat) == rownames(coldata)))

# ---- DESeq2 -----------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = mat,
  colData = coldata,
  design = as.formula(opt$design)
)
# Pre-filter very low-count peaks.
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c(opt$factor, opt$numerator, opt$denominator),
  alpha = opt$fdr
)
res <- res[order(res$padj), ]

# ---- write results ----------------------------------------------------------
res_df <- as.data.frame(res)
res_df$peak <- rownames(res_df)
res_df <- merge(peak_info, res_df, by = "peak")
res_df <- res_df[order(res_df$padj), ]
write.table(res_df, file.path(opt$outdir, "differential_accessibility.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# normalised counts
nc <- counts(dds, normalized = TRUE)
write.table(data.frame(peak = rownames(nc), nc, check.names = FALSE),
            file.path(opt$outdir, "normalized_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- up / down BEDs ---------------------------------------------------------
sig <- subset(res_df, !is.na(padj) & padj < opt$fdr &
                abs(log2FoldChange) >= opt$lfc)
up <- subset(sig, log2FoldChange > 0)
down <- subset(sig, log2FoldChange < 0)
write_bed <- function(d, path) {
  if (nrow(d) == 0) { file.create(path); return(invisible()) }
  write.table(d[, c("Chr", "Start", "End", "peak")], path,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
write_bed(up, file.path(opt$outdir, "up_peaks.bed"))
write_bed(down, file.path(opt$outdir, "down_peaks.bed"))

# ---- plots ------------------------------------------------------------------
pdf(file.path(opt$outdir, "MA_plot.pdf")); plotMA(res, ylim = c(-5, 5)); dev.off()

vsd <- vst(dds, blind = TRUE)
p <- plotPCA(vsd, intgroup = opt$factor) + theme_bw() +
  ggtitle("ATAC-seq accessibility PCA")
ggsave(file.path(opt$outdir, "PCA_plot.pdf"), p, width = 6, height = 5)

cat(sprintf("DESeq2 done: %d significant peaks (%d up, %d down) at FDR %.3g\n",
            nrow(sig), nrow(up), nrow(down), opt$fdr))
