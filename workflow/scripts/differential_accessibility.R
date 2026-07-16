#!/usr/bin/env Rscript
# =============================================================================
# Differential accessibility with DESeq2 over consensus ATAC peaks.
# Reads a featureCounts matrix + sample sheet, fits the configured design,
# extracts the requested contrast, applies adaptive LFC shrinkage (ashr) for
# ranking/plots, and writes:
#   - differential_accessibility.tsv   (full results, unshrunken + shrunken LFC)
#   - normalized_counts.tsv            (DESeq2 size-factor-normalised counts)
#   - tested_peaks.bed                 (finite-p-value opportunity universe)
#   - up_peaks.bed / down_peaks.bed    (significant peaks for motif/enrichment)
#   - diffacc_summary.tsv              (counts of tested/sig/up/down)
#   Analysis figures (vector PDF + PNG; headline ones also *_mqc.png):
#   - PCA_plot.pdf + scree_plot.pdf        sample structure
#   - sample_correlation_heatmap.pdf       Spearman corr of VST
#   - sample_distance_heatmap.pdf          Euclidean VST distances
#   - MA_plot.pdf                          shrunken LFC vs mean accessibility
#   - volcano_plot.pdf                     signed significance
#   - differential_peaks_heatmap.pdf       z-scored VST of top DA peaks
#   - pvalue_histogram.pdf                 p-value calibration diagnostic
# =============================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})
source("workflow/scripts/atac_theme.R")

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts"),
  make_option("--samples"),
  make_option("--design", default = "~condition"),
  make_option("--factor", default = "condition"),
  make_option("--numerator"),
  make_option("--denominator"),
  make_option("--fdr", type = "double", default = 0.05),
  make_option("--lfc", type = "double", default = 0),
  make_option("--topn", type = "integer", default = 60,
              help = "top significant peaks shown in the DA heatmap"),
  make_option("--label-n", type = "integer", default = 12, dest = "label_n",
              help = "peaks labelled on the volcano plot, split evenly up/down"),
  make_option("--ntop-pca", type = "integer", default = 2000, dest = "ntop_pca",
              help = "most-variable peaks used for PCA / correlation"),
  make_option("--outdir")
)))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
P <- function(f) file.path(opt$outdir, f)

# A matrix with no peaks means the upstream peak universe failed. Finding zero
# significant peaks after testing a non-empty matrix is a valid result.
fail_degenerate <- function(msg) {
  stop(sprintf(paste0(
    "%s.\n",
    "  DESeq2 has nothing to test. This means peak calling or the consensus step\n",
    "  produced an empty peak set -- check the configured consensus_peaks.bed\n",
    "  and the per-sample MACS3 logs. It is not a biological result."), msg),
    call. = FALSE)
}

# ---- featureCounts matrix ---------------------------------------------------
fc <- read.delim(opt$counts, comment.char = "#", check.names = FALSE)
if (nrow(fc) == 0) fail_degenerate("No consensus peaks")
feature_columns <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
if (ncol(fc) <= length(feature_columns) ||
    !identical(colnames(fc)[seq_along(feature_columns)], feature_columns)) {
  stop("merged count matrix lacks the six feature columns or any sample columns")
}

peak_info <- fc[, c("Geneid", "Chr", "Start", "End")]
colnames(peak_info)[1] <- "peak"
count_cols <- fc[, 7:ncol(fc), drop = FALSE]
count_samples <- colnames(count_cols)
if (anyDuplicated(count_samples) || any(!nzchar(count_samples))) {
  stop("merged count matrix sample identifiers must be non-empty and unique")
}
mat <- as.matrix(count_cols)
if (!is.numeric(mat) || any(!is.finite(mat)) || any(mat < 0) ||
    any(mat != floor(mat)) || any(mat > .Machine$integer.max)) {
  stop("merged count matrix values must be finite, non-negative integer counts")
}
storage.mode(mat) <- "integer"
rownames(mat) <- fc$Geneid

# ---- sample sheet -----------------------------------------------------------
coldata <- read.delim(opt$samples, comment.char = "#", stringsAsFactors = FALSE)
if (anyDuplicated(coldata$sample)) stop("sample sheet contains duplicate sample identifiers")
selected <- rep(TRUE, nrow(coldata))
if ("include" %in% colnames(coldata)) {
  include <- tolower(trimws(as.character(coldata$include)))
  include[is.na(include)] <- ""
  selected <- include %in% c("", "1", "true", "yes", "y")
}
expected_samples <- as.character(coldata$sample[selected])
if (!setequal(count_samples, expected_samples)) {
  stop(sprintf(
    paste0("merged count matrix and included sample sheet disagree; ",
           "missing from counts: [%s]; unexpected in counts: [%s]"),
    paste(setdiff(expected_samples, count_samples), collapse = ", "),
    paste(setdiff(count_samples, expected_samples), collapse = ", ")
  ))
}
coldata <- coldata[match(count_samples, coldata$sample), , drop = FALSE]
rownames(coldata) <- coldata$sample
stopifnot(identical(count_samples, rownames(coldata)))
coldata[[opt$factor]] <- factor(coldata[[opt$factor]])
coldata[[opt$factor]] <- relevel(coldata[[opt$factor]], ref = opt$denominator)
if ("replicate" %in% colnames(coldata)) {
  coldata$replicate <- factor(coldata$replicate)
}
stopifnot(identical(colnames(mat), rownames(coldata)))

# ---- DESeq2 -----------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = mat, colData = coldata,
                              design = as.formula(opt$design))
dds <- dds[rowSums(counts(dds)) >= 10, ]
if (nrow(dds) < 2) fail_degenerate("Too few peaks after count filtering")
dds <- DESeq(dds)

if (opt$lfc > 0) {
  res <- results(dds, contrast = c(opt$factor, opt$numerator, opt$denominator),
                 alpha = opt$fdr, lfcThreshold = opt$lfc,
                 altHypothesis = "greaterAbs")
} else {
  res <- results(dds, contrast = c(opt$factor, opt$numerator, opt$denominator),
                 alpha = opt$fdr)
}

# Adaptive shrinkage of LFCs (ashr supports arbitrary contrasts). It is a
# required, pinned part of this analysis: silently substituting unshrunken
# effects would mislabel the output column and MA/volcano axes.
res_shrunk <- lfcShrink(
  dds, contrast = c(opt$factor, opt$numerator, opt$denominator),
  type = "ashr", res = res
)

res_df <- as.data.frame(res)
res_df$peak <- rownames(res_df)
res_df$lfcShrink <- res_shrunk$log2FoldChange[match(res_df$peak, rownames(res_shrunk))]
res_df <- merge(peak_info, res_df, by = "peak")
res_df <- res_df[order(res_df$padj), ]
res_df <- res_df[, c("peak", "Chr", "Start", "End", "baseMean", "log2FoldChange",
                     "lfcShrink", "lfcSE", "stat", "pvalue", "padj")]

# Downstream over-representation tests and motif enrichment must use the loci
# that actually had an opportunity to be significant. Keep finite raw-p-value
# hypotheses (including independently filtered adjusted p-values), and exclude
# low-count prefilter failures and Cook's-distance outliers with no test.
tested <- res_df[is.finite(res_df$pvalue), , drop = FALSE]
if (nrow(tested) < 2L) fail_degenerate("Fewer than two peaks received finite DESeq2 p-values")
tested_bed <- tested[, c("Chr", "Start", "End", "peak")]
tested_bed$Start <- tested_bed$Start - 1L
write.table(tested_bed, P("tested_peaks.bed"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

nc <- counts(dds, normalized = TRUE)
write.table(data.frame(peak = rownames(nc), nc, check.names = FALSE),
            P("normalized_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- significance + up/down BEDs (ranked by shrunken LFC) -------------------
# One significance definition is shared by BEDs, summaries, and figures.
eligible <- !is.na(res_df$padj) & res_df$padj < opt$fdr &
  is.finite(res_df$log2FoldChange) & abs(res_df$log2FoldChange) >= opt$lfc
res_df$direction <- "ns"
res_df$direction[eligible & res_df$log2FoldChange > 0] <- "up"
res_df$direction[eligible & res_df$log2FoldChange < 0] <- "down"

# Written only now that `direction` exists, so the table carries the same call
# the BEDs and the figures make. Re-deriving it downstream means re-implementing
# the FDR and LFC thresholds, and the two drift.
write.table(res_df, P("differential_accessibility.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sig <- subset(res_df, direction != "ns")
up   <- subset(res_df, direction == "up")
down <- subset(res_df, direction == "down")
write_bed <- function(d, path) {
  if (nrow(d) == 0) { file.create(path); return(invisible()) }
  d <- d[order(-abs(d$lfcShrink)), ]
  # featureCounts reports SAF coordinates as 1-based inclusive. BED is 0-based
  # half-open, so convert the start back before downstream tools consume it.
  bed <- d[, c("Chr", "Start", "End", "peak")]
  bed$Start <- bed$Start - 1L
  write.table(bed, path,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
write_bed(up, P("up_peaks.bed")); write_bed(down, P("down_peaks.bed"))
write.table(data.frame(
  metric = c("tested", "adjusted_p_available", "significant", "up", "down"),
  value = c(nrow(tested), sum(!is.na(res_df$padj)),
            nrow(sig), nrow(up), nrow(down))),
  P("diffacc_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# =============================================================================
# Figures
# =============================================================================
cond_col <- atac_condition_colours(levels(coldata[[opt$factor]]))
contrast_lab <- sprintf("%s vs %s", opt$numerator, opt$denominator)

vsd <- tryCatch(vst(dds, blind = TRUE),
                error = function(e) varianceStabilizingTransformation(dds, blind = TRUE))
vmat <- assay(vsd)

# ---- PCA + scree (manual prcomp so both share one basis) --------------------
rv <- matrixStats::rowVars(vmat)
sel <- order(rv, decreasing = TRUE)[seq_len(min(opt$ntop_pca, sum(rv > 0)))]
pca <- prcomp(t(vmat[sel, , drop = FALSE]))
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pdat <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                   sample = colnames(vmat),
                   condition = coldata[[opt$factor]])
p_pca <- ggplot(pdat, aes(PC1, PC2, colour = condition)) +
  geom_point(size = 3.2, alpha = 0.9) +
  scale_colour_manual(values = cond_col, name = opt$factor) +
  labs(title = "ATAC accessibility PCA", subtitle = "variance-stabilised counts",
       x = sprintf("PC1 (%.1f%%)", pct[1]), y = sprintf("PC2 (%.1f%%)", pct[2])) +
  theme_atac()
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_pca <- p_pca + ggrepel::geom_text_repel(aes(label = sample), size = 2.6,
                                            colour = "#52514e", max.overlaps = 20,
                                            show.legend = FALSE)
} else {
  p_pca <- p_pca + geom_text(aes(label = sample), size = 2.4, vjust = -0.9,
                             colour = "#52514e", show.legend = FALSE)
}
save_fig(p_pca, P("PCA_plot.pdf"), width = 7, height = 5.5, mqc = TRUE)

scree <- data.frame(PC = factor(paste0("PC", seq_along(pct)),
                                levels = paste0("PC", seq_along(pct))),
                    pct = pct)[seq_len(min(10, length(pct))), ]
p_scree <- ggplot(scree, aes(PC, pct)) +
  geom_col(fill = ATAC_CATEGORICAL[1], width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", pct)), vjust = -0.4, size = 2.8,
            colour = "#52514e") +
  labs(title = "Scree plot", x = NULL, y = "Variance explained (%)") +
  theme_atac()
save_fig(p_scree, P("scree_plot.pdf"), width = 6.5, height = 4.5)

# ---- sample-sample Spearman correlation heatmap -----------------------------
ann <- data.frame(coldata[[opt$factor]], check.names = FALSE)
colnames(ann) <- opt$factor
rownames(ann) <- colnames(vmat)
ann_colours <- setNames(list(cond_col), opt$factor)
cormat <- cor(vmat, method = "spearman")
save_base_fig(function() {
  grid::grid.draw(pheatmap(cormat, annotation_col = ann, annotation_row = ann,
           annotation_colors = ann_colours,
           color = ATAC_SEQ(100), border_color = NA,
           main = "Sample-sample correlation (Spearman, VST)",
           display_numbers = ncol(vmat) <= 16,
           number_color = "#0b0b0b", fontsize_number = 6, silent = TRUE)$gtable)
}, P("sample_correlation_heatmap.pdf"), width = 7.5, height = 6.5, mqc = TRUE)

# ---- sample distance heatmap ------------------------------------------------
sampleDist <- as.matrix(dist(t(vmat)))
save_base_fig(function() {
  grid::grid.draw(pheatmap(sampleDist, annotation_col = ann,
           annotation_colors = ann_colours,
           color = rev(ATAC_SEQ(100)), border_color = NA,
           main = "Sample distances (Euclidean, VST)", silent = TRUE)$gtable)
}, P("sample_distance_heatmap.pdf"), width = 7.5, height = 6.5)

# ---- MA plot (shrunken LFC) -------------------------------------------------
ma <- res_df
ma$sig <- ma$direction
ma <- ma[is.finite(ma$baseMean) & ma$baseMean > 0, ]
p_ma <- ggplot(ma[order(ma$sig != "ns"), ],
               aes(baseMean, lfcShrink, colour = sig)) +
  geom_hline(yintercept = 0, colour = "#9a988f", linewidth = 0.4) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_x_log10() +
  scale_colour_manual(values = ATAC_STATUS,
                      breaks = c("up", "down", "ns"),
                      labels = c(sprintf("up (%d)", nrow(up)),
                                 sprintf("down (%d)", nrow(down)), "n.s."),
                      name = NULL) +
  labs(title = "MA plot", subtitle = contrast_lab,
       x = "mean normalised accessibility", y = "log2 fold-change (shrunken)") +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  theme_atac()
save_fig(p_ma, P("MA_plot.pdf"), width = 6.5, height = 5, mqc = TRUE)

# ---- volcano plot -----------------------------------------------------------
vol <- res_df[!is.na(res_df$padj), ]
if (nrow(vol) > 0) {
  positive_padj <- vol$padj[is.finite(vol$padj) & vol$padj > 0]
  padj_floor <- if (length(positive_padj)) {
    max(min(positive_padj) / 10, .Machine$double.xmin)
  } else {
    .Machine$double.xmin
  }
  vol$logp <- -log10(pmax(vol$padj, padj_floor))
  vol$sig <- vol$direction
  vol$locus <- sprintf("%s:%d", vol$Chr, vol$Start)
  # Label each direction separately so one well-powered direction cannot occupy
  # every annotation slot.
  per_side <- max(1L, opt$label_n %/% 2L)
  lab <- do.call(rbind, lapply(c("up", "down"), function(dir) {
    side <- vol[vol$sig == dir, ]
    head(side[order(side$padj), ], per_side)
  }))
  call_lab <- sprintf("Calls use padj < %.3g and |raw LFC| >= %.3g",
                      opt$fdr, opt$lfc)
  p_vol <- ggplot(vol, aes(lfcShrink, logp, colour = sig)) +
    geom_hline(yintercept = -log10(opt$fdr), linetype = "dashed",
               colour = "#9a988f", linewidth = 0.4) +
    geom_point(size = 0.7, alpha = 0.6) +
    scale_colour_manual(values = ATAC_STATUS, breaks = c("up", "down", "ns"),
                        labels = c(sprintf("up (%d)", nrow(up)),
                                   sprintf("down (%d)", nrow(down)), "n.s."),
                        name = NULL) +
    labs(title = "Differential accessibility",
         subtitle = paste(contrast_lab, call_lab, sep = "\n"),
         x = "log2 fold-change (shrunken)", y = "-log10 adjusted p") +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    theme_atac()
  if (nrow(lab) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    p_vol <- p_vol + ggrepel::geom_text_repel(
      data = lab, aes(label = locus), size = 2.3, colour = "#2b2b2b",
      max.overlaps = 30, min.segment.length = 0, show.legend = FALSE)
  }
  save_fig(p_vol, P("volcano_plot.pdf"), width = 6.5, height = 6, mqc = TRUE)

} else {
  placeholder_fig(P("volcano_plot.pdf"), "No peaks with adjusted p-values", mqc = TRUE)
}

# Use every finite raw p-value. This diagnostic is deliberately independent of
# adjusted-p-value availability because DESeq2 independent filtering may leave
# valid tested hypotheses with NA padj values.
pvals <- res_df[is.finite(res_df$pvalue), , drop = FALSE]
if (nrow(pvals) > 0) {
  p_ph <- ggplot(pvals, aes(pvalue)) +
    geom_histogram(breaks = seq(0, 1, 0.025), fill = ATAC_CATEGORICAL[1],
                   colour = "white", linewidth = 0.2) +
    labs(title = "P-value distribution", x = "p-value", y = "peaks") +
    theme_atac()
  save_fig(p_ph, P("pvalue_histogram.pdf"), width = 6, height = 4.5)
} else {
  placeholder_fig(P("pvalue_histogram.pdf"), "No finite raw p-values")
}

# ---- heatmap of top differential peaks (z-scored VST) -----------------------
# Balance the heatmap across directions so both sides remain visible when their
# statistical power differs.
top_sig <- do.call(rbind, lapply(c("up", "down"), function(dir) {
  side <- sig[sig$direction == dir, ]
  head(side[order(side$padj), ], max(1L, opt$topn %/% 2L))
}))
if (nrow(top_sig) >= 2) {
  hm <- vmat[top_sig$peak, , drop = FALSE]
  hm <- t(scale(t(hm)))                       # z-score per peak
  hm[!is.finite(hm)] <- 0
  lim <- max(abs(quantile(hm, c(0.01, 0.99), na.rm = TRUE)))
  brk <- seq(-lim, lim, length.out = 101)
  save_base_fig(function() {
    grid::grid.draw(pheatmap(hm, annotation_col = ann,
             annotation_colors = ann_colours,
             color = ATAC_DIVERGING(100), breaks = brk, border_color = NA,
             show_rownames = FALSE, cluster_cols = TRUE,
             clustering_method = "ward.D2",
             main = sprintf("Top %d peaks per direction (z-scored VST)\n%s",
                            max(sum(top_sig$direction == "up"),
                                sum(top_sig$direction == "down")),
                            contrast_lab), silent = TRUE)$gtable)
  }, P("differential_peaks_heatmap.pdf"), width = 7.5, height = 8, mqc = TRUE)
} else {
  placeholder_fig(P("differential_peaks_heatmap.pdf"),
                  "Fewer than 2 significant peaks", mqc = TRUE)
}

cat(sprintf("DESeq2 done: %d significant peaks (%d up, %d down) at FDR %.3g\n",
            nrow(sig), nrow(up), nrow(down), opt$fdr))
