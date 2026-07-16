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
  make_option("--taxid", type = "integer", default = 0,
              help = "NCBI taxonomy id for JASPAR motif selection (9606 human, 10090 mouse, 10116 rat)"),
  make_option("--out-dev", dest = "out_dev"),
  make_option("--out-var", dest = "out_var"),
  make_option("--out-var-plot", dest = "out_var_plot", default = ""),
  make_option("--out-heatmap", dest = "out_heatmap", default = ""),
  make_option("--top-n", type = "integer", default = 30, dest = "top_n"),
  make_option("--seed", type = "integer", default = 1,
              help = "RNG seed; chromVAR samples background peaks and bootstraps")
)))

# chromVAR draws GC/depth-matched backgrounds and bootstraps variability
# intervals. A fixed seed makes repeated runs reproducible.
set.seed(opt$seed)
# The environment's automatic BiocParallel backend can select a socket cluster,
# which is fragile on restricted compute nodes and does not make the bootstrap
# RNG contract explicit. A seeded serial backend is deterministic and portable;
# users can scale this only after validating a seeded parallel backend locally.
BiocParallel::register(
  BiocParallel::SerialParam(RNGseed = opt$seed, progressbar = FALSE),
  default = TRUE
)

if (is.na(opt$taxid) || opt$taxid < 1L) {
  stop("chromVAR requires a positive NCBI taxonomy id; got ", opt$taxid)
}
if (!nzchar(opt$bsgenome) || !requireNamespace(opt$bsgenome, quietly = TRUE)) {
  stop("chromVAR requires the BSgenome package '", opt$bsgenome,
       "'. Install it (see envs/r.yaml) and re-run.")
}
library(opt$bsgenome, character.only = TRUE)
genome <- get(opt$bsgenome)

# ---- peaks & counts ---------------------------------------------------------
bed <- read.delim(opt$bed, header = FALSE)
if (ncol(bed) < 4L) stop("consensus BED must contain peak identifiers in column 4")

fc <- read.delim(opt$counts, comment.char = "#", check.names = FALSE)
feature_columns <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
if (ncol(fc) <= length(feature_columns) ||
    !identical(colnames(fc)[seq_along(feature_columns)], feature_columns)) {
  stop("merged count matrix lacks the six feature columns or any sample columns")
}
bed_ids <- as.character(bed[[4]])
count_ids <- as.character(fc$Geneid)
if (anyDuplicated(bed_ids) || anyDuplicated(count_ids)) {
  stop("consensus BED and featureCounts Geneid values must each be unique")
}
if (!setequal(bed_ids, count_ids)) {
  stop(sprintf(
    paste0("consensus BED and featureCounts peak IDs disagree: ",
           "%d count IDs missing from BED; %d BED IDs missing from counts"),
    length(setdiff(count_ids, bed_ids)), length(setdiff(bed_ids, count_ids))
  ))
}
# featureCounts commonly preserves SAF order, but chromVAR must never rely on
# that implementation detail: align genomic ranges explicitly to count rows.
bed <- bed[match(count_ids, bed_ids), , drop = FALSE]
stopifnot(identical(as.character(bed[[4]]), count_ids))
peaks <- GRanges(bed[[1]], IRanges(bed[[2]] + 1, bed[[3]]), peak = bed[[4]])
counts <- as.matrix(fc[, 7:ncol(fc), drop = FALSE])
if (!is.numeric(counts) || any(!is.finite(counts)) || any(counts < 0) ||
    any(counts != floor(counts))) {
  stop("merged count matrix values must be finite, non-negative integer counts")
}
count_samples <- colnames(counts)
if (anyDuplicated(count_samples) || any(!nzchar(count_samples))) {
  stop("merged count matrix sample identifiers must be non-empty and unique")
}
sample_table <- read.delim(opt$samples, comment.char = "#", stringsAsFactors = FALSE)
if (anyDuplicated(sample_table$sample)) {
  stop("sample sheet contains duplicate sample identifiers")
}
selected <- rep(TRUE, nrow(sample_table))
if ("include" %in% colnames(sample_table)) {
  include <- tolower(trimws(as.character(sample_table$include)))
  include[is.na(include)] <- ""
  selected <- include %in% c("", "1", "true", "yes", "y")
}
expected_samples <- as.character(sample_table$sample[selected])
if (!setequal(count_samples, expected_samples)) {
  stop(sprintf(
    paste0("merged count matrix and included sample sheet disagree; ",
           "missing from counts: [%s]; unexpected in counts: [%s]"),
    paste(setdiff(expected_samples, count_samples), collapse = ", "),
    paste(setdiff(count_samples, expected_samples), collapse = ", ")
  ))
}
sample_table <- sample_table[match(count_samples, sample_table$sample), , drop = FALSE]
stopifnot(identical(count_samples, as.character(sample_table$sample)))
rownames(counts) <- count_ids

se <- SummarizedExperiment(assays = list(counts = counts), rowRanges = peaks)
# Harmonise names without assuming a human/mouse standard-chromosome set. Keep
# direct custom contig names when they work; otherwise try BSgenome styles on
# copies and retain the representation that maps the largest number of peaks.
n_before <- nrow(se)
genome_levels <- GenomeInfoDb::seqlevels(genome)
shared_peak_count <- function(x) {
  sum(as.character(GenomicRanges::seqnames(SummarizedExperiment::rowRanges(x))) %in%
        genome_levels)
}
best <- se
best_mode <- "direct chromosome names"
best_n <- shared_peak_count(best)
target_styles <- tryCatch(
  GenomeInfoDb::seqlevelsStyle(genome),
  error = function(e) character()
)
for (style in unique(target_styles)) {
  candidate <- tryCatch({
    x <- se
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
shared <- intersect(GenomeInfoDb::seqlevels(best), genome_levels)
if (best_n == 0L || length(shared) == 0L) {
  stop(sprintf(
    "chromVAR peaks and BSgenome share no contigs. peaks: %s ... BSgenome: %s ...",
    paste(utils::head(GenomeInfoDb::seqlevels(se), 3), collapse = ", "),
    paste(utils::head(genome_levels, 3), collapse = ", ")
  ))
}
se <- GenomeInfoDb::keepSeqlevels(best, shared, pruning.mode = "coarse")
if (nrow(se) == 0L) stop("chromVAR retained zero peaks on BSgenome contigs")
retained_fraction <- nrow(se) / n_before
if (retained_fraction < 0.5) {
  stop(sprintf(
    paste0("chromVAR retained only %d/%d peaks (%.1f%%) on BSgenome contigs; ",
           "refusing a severely biased subset"),
    nrow(se), n_before, 100 * retained_fraction
  ))
}
if (retained_fraction < 0.8) {
  warning(sprintf("chromVAR retained %d/%d peaks (%.1f%%) on BSgenome contigs",
                  nrow(se), n_before, 100 * retained_fraction))
}
message(sprintf(
  "chromVAR chromosome naming: %s; %d shared contigs; retained %d/%d peaks (%.1f%%)",
  best_mode, length(shared), nrow(se), n_before, 100 * retained_fraction
))
se <- addGCBias(se, genome = genome)
se <- filterPeaks(se, non_overlapping = TRUE)

# ---- motifs -----------------------------------------------------------------
motifs <- getMatrixSet(JASPAR2020,
                       list(species = opt$taxid, collection = "CORE"))
if (length(motifs) == 0L) {
  stop("JASPAR2020 returned no CORE motifs for taxid ", opt$taxid)
}
motif_hits <- matchMotifs(motifs, se, genome = genome)

dev <- computeDeviations(object = se, annotations = motif_hits)
variability <- computeVariability(dev)
variability$motif <- rownames(variability)
variability <- variability[, c("motif", setdiff(colnames(variability), "motif")),
                           drop = FALSE]

write.table(data.frame(motif = rownames(deviations(dev)),
                       deviations(dev), check.names = FALSE),
            opt$out_dev, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(variability, opt$out_var, sep = "\t", quote = FALSE,
            row.names = FALSE)

# Variability figure: the most variable motifs, with their bootstrap intervals.
#
# chromVAR's own plotVariability() draws an error bar for every motif. With the
# full JASPAR CORE set that is 600+ vertical bars packed into the left edge of the
# panel: they merge into a solid block, the tail is unreadable, and the handful of
# labels it prints land on top of one another.
#
# Drawn horizontally, one row per motif, the intervals are legible -- and they
# need to be legible, because they are wide. Variability is estimated across
# samples, so with a handful of libraries the bootstrap lower bound sits near zero
# even for the top motif. The ranking is informative; the spacing between adjacent
# ranks is not. Showing the intervals is what stops a reader over-reading the
# order. The full table is written to chromvar_variability.tsv.
if (nzchar(opt$out_var_plot)) {
  vd <- variability[order(-variability$variability), , drop = FALSE]
  n_top <- min(15L, nrow(vd))
  top <- vd[seq_len(n_top), , drop = FALSE]
  top$name <- factor(top$name, levels = rev(top$name))

  vp <- ggplot(top, aes(variability, name)) +
    geom_errorbarh(aes(xmin = bootstrap_lower_bound, xmax = bootstrap_upper_bound),
                   height = 0, colour = "#9ec5f4", linewidth = 0.9) +
    geom_point(colour = "#184f95", size = 2.2) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Most variable motifs (chromVAR)",
         subtitle = sprintf("top %d of %d JASPAR2020 CORE motifs; bars are bootstrap 95%% CI across %d samples",
                            n_top, nrow(vd), ncol(deviations(dev))),
         x = "Variability", y = NULL) +
    theme_atac(base_size = 10)
  suppressWarnings(save_fig(vp, opt$out_var_plot, width = 7, height = 5.5, mqc = TRUE))
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
  coldata <- sample_table
  if (all(c("sample", "condition") %in% colnames(coldata))) {
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
