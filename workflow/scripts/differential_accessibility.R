#!/usr/bin/env Rscript
# =============================================================================
# Differential accessibility with DESeq2 over consensus ATAC peaks.
# Reads a featureCounts matrix + sample sheet, fits the configured design,
# extracts the requested contrast, applies adaptive LFC shrinkage (ashr) for
# ranking/plots, and writes:
#   - differential_accessibility.tsv   (full results, unshrunken + shrunken LFC)
#   - normalized_counts.tsv            (DESeq2 size-factor-normalised counts)
#   - up_peaks.bed / down_peaks.bed    (significant peaks for motif/enrichment)
#   - diffacc_summary.tsv              (counts of tested/sig/up/down)
#   Publication-grade figures (vector PDF + PNG; headline ones also *_mqc.png):
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
  make_option("--label-n", type = "integer", default = 20, dest = "label_n",
              help = "top peaks labelled on the volcano plot"),
  make_option("--ntop-pca", type = "integer", default = 2000, dest = "ntop_pca",
              help = "most-variable peaks used for PCA / correlation"),
  make_option("--outdir")
)))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
P <- function(f) file.path(opt$outdir, f)

# All figure/table targets, so we can always create them (even when degenerate).
fig_pdfs <- c("PCA_plot.pdf", "scree_plot.pdf", "sample_correlation_heatmap.pdf",
              "sample_distance_heatmap.pdf", "MA_plot.pdf", "volcano_plot.pdf",
              "differential_peaks_heatmap.pdf", "pvalue_histogram.pdf")

emit_degenerate <- function(msg) {
  empty_res <- data.frame(
    peak = character(), Chr = character(), Start = integer(), End = integer(),
    baseMean = numeric(), log2FoldChange = numeric(), lfcShrink = numeric(),
    lfcSE = numeric(), stat = numeric(), pvalue = numeric(), padj = numeric()
  )
  write.table(empty_res, P("differential_accessibility.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(peak = character()), P("normalized_counts.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(metric = c("tested", "significant", "up", "down"),
                         value = c(0L, 0L, 0L, 0L)),
              P("diffacc_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  file.create(P("up_peaks.bed")); file.create(P("down_peaks.bed"))
  # These five must always emit a *_mqc.png (declared MultiQC inputs).
  mqc_figs <- c("PCA_plot.pdf", "MA_plot.pdf", "volcano_plot.pdf",
                "sample_correlation_heatmap.pdf", "differential_peaks_heatmap.pdf")
  for (f in fig_pdfs) placeholder_fig(P(f), msg, mqc = f %in% mqc_figs)
  cat(sprintf("DESeq2 skipped: %s\n", msg))
  quit(save = "no", status = 0)
}

# ---- featureCounts matrix ---------------------------------------------------
fc <- read.delim(opt$counts, comment.char = "#", check.names = FALSE)
if (nrow(fc) == 0) emit_degenerate("No consensus peaks")

peak_info <- fc[, c("Geneid", "Chr", "Start", "End")]
colnames(peak_info)[1] <- "peak"
count_cols <- fc[, 7:ncol(fc), drop = FALSE]
clean <- sub("\\.filtered\\.bam$", "", basename(colnames(count_cols)))
colnames(count_cols) <- clean
mat <- as.matrix(count_cols)
rownames(mat) <- fc$Geneid

# ---- sample sheet -----------------------------------------------------------
coldata <- read.delim(opt$samples, comment.char = "#", stringsAsFactors = FALSE)
rownames(coldata) <- coldata$sample
coldata <- coldata[colnames(mat), , drop = FALSE]
coldata[[opt$factor]] <- factor(coldata[[opt$factor]])
coldata[[opt$factor]] <- relevel(coldata[[opt$factor]], ref = opt$denominator)
if ("replicate" %in% colnames(coldata)) {
  coldata$replicate <- factor(coldata$replicate)
}
stopifnot(all(colnames(mat) == rownames(coldata)))

# ---- DESeq2 -----------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = mat, colData = coldata,
                              design = as.formula(opt$design))
dds <- dds[rowSums(counts(dds)) >= 10, ]
if (nrow(dds) < 2) emit_degenerate("Too few peaks after count filtering")
dds <- DESeq(dds)

res <- results(dds, contrast = c(opt$factor, opt$numerator, opt$denominator),
               alpha = opt$fdr)

# Adaptive shrinkage of LFCs (ashr supports arbitrary contrasts). Used for
# ranking + MA/volcano so weakly-supported peaks are pulled toward zero.
res_shrunk <- tryCatch(
  lfcShrink(dds, contrast = c(opt$factor, opt$numerator, opt$denominator),
            type = "ashr", res = res),
  error = function(e) { message("lfcShrink failed (", conditionMessage(e),
                                "); using unshrunken LFC."); res })

res_df <- as.data.frame(res)
res_df$peak <- rownames(res_df)
res_df$lfcShrink <- res_shrunk$log2FoldChange[match(res_df$peak, rownames(res_shrunk))]
res_df <- merge(peak_info, res_df, by = "peak")
res_df <- res_df[order(res_df$padj), ]
res_df <- res_df[, c("peak", "Chr", "Start", "End", "baseMean", "log2FoldChange",
                     "lfcShrink", "lfcSE", "stat", "pvalue", "padj")]
write.table(res_df, P("differential_accessibility.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

nc <- counts(dds, normalized = TRUE)
write.table(data.frame(peak = rownames(nc), nc, check.names = FALSE),
            P("normalized_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- significance + up/down BEDs (ranked by shrunken LFC) -------------------
sig <- subset(res_df, !is.na(padj) & padj < opt$fdr & abs(log2FoldChange) >= opt$lfc)
up   <- subset(sig, log2FoldChange > 0)
down <- subset(sig, log2FoldChange < 0)
write_bed <- function(d, path) {
  if (nrow(d) == 0) { file.create(path); return(invisible()) }
  d <- d[order(-abs(d$lfcShrink)), ]
  write.table(d[, c("Chr", "Start", "End", "peak")], path,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
write_bed(up, P("up_peaks.bed")); write_bed(down, P("down_peaks.bed"))
write.table(data.frame(
  metric = c("tested", "significant", "up", "down"),
  value = c(sum(!is.na(res_df$padj)), nrow(sig), nrow(up), nrow(down))),
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
  geom_col(fill = "#2a78d6", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", pct)), vjust = -0.4, size = 2.8,
            colour = "#52514e") +
  labs(title = "Scree plot", x = NULL, y = "Variance explained (%)") +
  theme_atac()
save_fig(p_scree, P("scree_plot.pdf"), width = 6.5, height = 4.5)

# ---- sample-sample Spearman correlation heatmap -----------------------------
ann <- data.frame(condition = coldata[[opt$factor]])
rownames(ann) <- colnames(vmat)
ann_colours <- list(condition = cond_col)
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
ma$sig <- ifelse(is.na(ma$padj) | ma$padj >= opt$fdr, "ns",
                 ifelse(ma$lfcShrink > 0, "up", "down"))
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
  minp <- min(vol$padj[vol$padj > 0], na.rm = TRUE)
  vol$logp <- -log10(pmax(vol$padj, minp / 10))
  vol$sig <- ifelse(vol$padj < opt$fdr & vol$lfcShrink >= opt$lfc, "up",
             ifelse(vol$padj < opt$fdr & vol$lfcShrink <= -opt$lfc, "down", "ns"))
  vol$locus <- sprintf("%s:%d", vol$Chr, vol$Start)
  lab <- head(vol[vol$sig != "ns", ][order(vol[vol$sig != "ns", ]$padj), ], opt$label_n)
  p_vol <- ggplot(vol, aes(lfcShrink, logp, colour = sig)) +
    geom_hline(yintercept = -log10(opt$fdr), linetype = "dashed",
               colour = "#9a988f", linewidth = 0.4) +
    { if (opt$lfc > 0) geom_vline(xintercept = c(-opt$lfc, opt$lfc),
                                  linetype = "dashed", colour = "#9a988f",
                                  linewidth = 0.4) else NULL } +
    geom_point(size = 0.7, alpha = 0.6) +
    scale_colour_manual(values = ATAC_STATUS, breaks = c("up", "down", "ns"),
                        labels = c(sprintf("up (%d)", nrow(up)),
                                   sprintf("down (%d)", nrow(down)), "n.s."),
                        name = NULL) +
    labs(title = "Differential accessibility", subtitle = contrast_lab,
         x = "log2 fold-change (shrunken)", y = "-log10 adjusted p") +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    theme_atac()
  if (nrow(lab) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    p_vol <- p_vol + ggrepel::geom_text_repel(
      data = lab, aes(label = locus), size = 2.3, colour = "#2b2b2b",
      max.overlaps = 30, min.segment.length = 0, show.legend = FALSE)
  }
  save_fig(p_vol, P("volcano_plot.pdf"), width = 6.5, height = 6, mqc = TRUE)

  p_ph <- ggplot(vol, aes(pvalue)) +
    geom_histogram(breaks = seq(0, 1, 0.025), fill = "#2a78d6",
                   colour = "white", linewidth = 0.2) +
    labs(title = "P-value distribution", x = "p-value", y = "peaks") +
    theme_atac()
  save_fig(p_ph, P("pvalue_histogram.pdf"), width = 6, height = 4.5)
} else {
  placeholder_fig(P("volcano_plot.pdf"), "No testable peaks", mqc = TRUE)
  placeholder_fig(P("pvalue_histogram.pdf"), "No testable peaks")
}

# ---- heatmap of top differential peaks (z-scored VST) -----------------------
top_sig <- head(sig[order(sig$padj), ], opt$topn)
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
             main = sprintf("Top %d differential peaks (z-scored VST)",
                            nrow(top_sig)), silent = TRUE)$gtable)
  }, P("differential_peaks_heatmap.pdf"), width = 7.5, height = 8, mqc = TRUE)
} else {
  placeholder_fig(P("differential_peaks_heatmap.pdf"),
                  "Fewer than 2 significant peaks", mqc = TRUE)
}

cat(sprintf("DESeq2 done: %d significant peaks (%d up, %d down) at FDR %.3g\n",
            nrow(sig), nrow(up), nrow(down), opt$fdr))
