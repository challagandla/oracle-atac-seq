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
  library(pheatmap)
})
source("workflow/scripts/atac_theme.R")

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts"),
  make_option("--bed"),
  make_option("--samples"),
  make_option("--bsgenome", default = ""),
  make_option("--taxid", type = "integer", default = 9606,
              help = "NCBI taxonomy id for JASPAR motif selection (9606 human, 10090 mouse, 10116 rat)"),
  make_option("--out-dev", dest = "out_dev"),
  make_option("--out-var", dest = "out_var"),
  make_option("--out-var-plot", dest = "out_var_plot", default = ""),
  make_option("--out-heatmap", dest = "out_heatmap", default = ""),
  make_option("--top-n", type = "integer", default = 30, dest = "top_n")
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
# Harmonise chromosome naming to the BSgenome (peaks use the genome's native
# names, e.g. Ensembl "1"; a UCSC BSgenome uses "chr1") and drop scaffolds the
# BSgenome doesn't carry, so GC bias + motif matching find the sequences.
suppressWarnings({
  GenomeInfoDb::seqlevelsStyle(se) <- GenomeInfoDb::seqlevelsStyle(genome)[1]
  se <- GenomeInfoDb::keepStandardChromosomes(se, pruning.mode = "coarse")
})
se <- addGCBias(se, genome = genome)
se <- filterPeaks(se, non_overlapping = TRUE)

# ---- motifs -----------------------------------------------------------------
motifs <- getMatrixSet(JASPAR2020,
                       list(species = opt$taxid, collection = "CORE"))
motif_hits <- matchMotifs(motifs, se, genome = genome)

dev <- computeDeviations(object = se, annotations = motif_hits)
variability <- computeVariability(dev)

write.table(data.frame(motif = rownames(deviations(dev)),
                       deviations(dev), check.names = FALSE),
            opt$out_dev, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(variability, opt$out_var, sep = "\t", quote = FALSE,
            row.names = FALSE)

# Variability figure: top motifs ranked by accessibility variability.
if (nzchar(opt$out_var_plot)) {
  vp <- plotVariability(variability, use_plotly = FALSE) + theme_atac(base_size = 10)
  suppressWarnings(save_fig(vp, opt$out_var_plot, width = 8, height = 5, mqc = TRUE))
}

# Heatmap of the most variable motifs (bias-corrected deviation z-scores) across
# samples, annotated by condition — the standard chromVAR summary figure.
if (nzchar(opt$out_heatmap)) {
  vv <- variability[order(-variability$variability), , drop = FALSE]
  top_ids <- head(rownames(vv), opt$top_n)
  dm <- deviations(dev)[top_ids, , drop = FALSE]
  # Label rows by TF name where available, keeping ids unique.
  labs <- vv[top_ids, "name"]
  rownames(dm) <- ifelse(is.na(labs) | !nzchar(labs), top_ids,
                         paste0(labs, " (", top_ids, ")"))
  ann <- NULL; ann_colours <- NA
  coldata <- tryCatch(read.delim(opt$samples, comment.char = "#",
                                 stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(coldata) && all(c("sample", "condition") %in% colnames(coldata))) {
    rownames(coldata) <- coldata$sample
    common <- intersect(colnames(dm), rownames(coldata))
    if (length(common)) {
      ann <- data.frame(condition = coldata[common, "condition"])
      rownames(ann) <- common
      ann_colours <- list(condition = atac_condition_colours(ann$condition))
    }
  }
  lim <- max(abs(quantile(dm, c(0.02, 0.98), na.rm = TRUE)))
  brk <- seq(-lim, lim, length.out = 101)
  save_base_fig(function() {
    grid::grid.draw(pheatmap(dm, annotation_col = ann,
             annotation_colors = ann_colours,
             color = ATAC_DIVERGING(100), breaks = brk, border_color = NA,
             main = sprintf("Top %d variable TF motifs (chromVAR)", nrow(dm)),
             fontsize_row = 7, silent = TRUE)$gtable)
  }, opt$out_heatmap, width = 8, height = 9, mqc = TRUE)
}
cat("chromVAR deviations written.\n")
