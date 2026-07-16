#!/usr/bin/env Rscript
# =============================================================================
# Functional enrichment of genes near differentially-accessible peaks.
# Uses the ChIPseeker annotation (peak -> nearest gene) to build up/down gene
# sets, then runs GO (enrichGO, offline via OrgDb) and optionally KEGG
# (enrichKEGG, online) with the full annotated gene set as the universe.
# Emits per-set result tables and a multi-panel dot-plot PDF.
# =============================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
})
source("workflow/scripts/atac_theme.R")

opt <- parse_args(OptionParser(option_list = list(
  make_option("--annotation"),                       # ChIPseeker annotated TSV
  make_option("--tested"),                           # DESeq2 opportunity BED
  make_option("--up"), make_option("--down"),        # up/down peak BEDs
  make_option("--orgdb", default = ""),
  make_option("--taxid", type = "integer", default = 0),
  make_option("--ontologies", default = "BP,MF"),
  make_option("--kegg", default = "false"),
  make_option("--qvalue", type = "double", default = 0.05),
  make_option("--top-n", type = "integer", default = 20, dest = "top_n"),
  make_option("--out-pdf", dest = "out_pdf"),
  make_option("--out-status", dest = "out_status"),
  make_option("--out-dir", dest = "out_dir")
)))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
kegg_org <- c("9606" = "hsa", "10090" = "mmu", "10116" = "rno")[as.character(opt$taxid)]
onts <- trimws(strsplit(opt$ontologies, ",")[[1]])
if (tolower(opt$kegg) == "true" && is.na(kegg_org)) {
  stop("KEGG enrichment supports taxid 9606, 10090, or 10116; got ",
       opt$taxid, call. = FALSE)
}

have_org <- nzchar(opt$orgdb) && requireNamespace(opt$orgdb, quietly = TRUE)
if (have_org) suppressPackageStartupMessages(library(opt$orgdb, character.only = TRUE))
if (!have_org) {
  stop("functional enrichment requires a matching OrgDb to validate Entrez IDs: ",
       opt$orgdb,
       call. = FALSE)
}
orgdb <- if (have_org) get(opt$orgdb) else NULL
if (have_org && !"ENTREZID" %in% AnnotationDbi::keytypes(orgdb)) {
  stop("OrgDb does not expose ENTREZID keys: ", opt$orgdb, call. = FALSE)
}
valid_entrez_ids <- as.character(AnnotationDbi::keys(orgdb, keytype = "ENTREZID"))
if (length(valid_entrez_ids) == 0L) {
  stop("OrgDb exposes no ENTREZID values: ", opt$orgdb, call. = FALSE)
}

# ---- map peaks -> Entrez gene ids via the annotation ------------------------
# The annotation is a required input, not an optional one. Swallowing a read
# error here yields an empty gene set, an empty dot-plot, and a zero exit status:
# indistinguishable from "nothing was enriched".
ann <- read.delim(opt$annotation, stringsAsFactors = FALSE)
if (nrow(ann) == 0L) {
  stop("peak annotation table is empty: ", opt$annotation, call. = FALSE)
}

peak_col <- intersect(c("V4", "name", "peak"), colnames(ann))
gene_col <- intersect(c("ENTREZID", "geneId", "geneID"), colnames(ann))
if (!length(peak_col) || !length(gene_col)) {
  stop(sprintf(paste0(
    "annotation lacks a peak-name and/or gene-id column.\n",
    "  found: %s\n",
    "  a gene id column requires ChIPseeker to have run with an OrgDb (annotation.orgdb)."),
    paste(colnames(ann), collapse = ", ")), call. = FALSE)
}

tested_bed <- read.delim(opt$tested, header = FALSE, stringsAsFactors = FALSE)
if (nrow(tested_bed) == 0L || ncol(tested_bed) < 4L) {
  stop("tested peak BED is empty or lacks column-4 peak identifiers: ",
       opt$tested, call. = FALSE)
}
tested_peaks <- unique(as.character(tested_bed[[4]]))
tested_peaks <- tested_peaks[!is.na(tested_peaks) & nzchar(tested_peaks)]
if (length(tested_peaks) == 0L) {
  stop("tested peak BED contains no usable peak identifiers: ", opt$tested,
       call. = FALSE)
}

ann$.peak <- as.character(ann[[peak_col[1]]])
source_gene <- trimws(as.character(ann[[gene_col[1]]]))
present_gene <- !is.na(source_gene) & nzchar(source_gene)
entrez_like <- present_gene & grepl("^[0-9]+$", source_gene) &
  source_gene %in% valid_entrez_ids
# Missing annotation rows and missing/non-Entrez IDs are mapping failures. The
# denominator is the exact set of finite-p-value DESeq2 hypotheses, not every
# consensus peak and not just the small subset with a non-missing gene field.
mapped_tested_peaks <- unique(ann$.peak[
  ann$.peak %in% tested_peaks & entrez_like
])
mapping_coverage <- length(mapped_tested_peaks) / length(tested_peaks)
if (mapping_coverage == 0) {
  stop(sprintf(
    paste0("annotation gene column '%s' has no Entrez IDs. Supply a matching ",
           "OrgDb to annotation so it emits ENTREZID before enrichment."),
    gene_col[1]
  ), call. = FALSE)
}
if (mapping_coverage < 0.5) {
  stop(sprintf(
    "only %.1f%% of finite-p-value tested peaks map to Entrez IDs; refusing biased enrichment",
    100 * mapping_coverage
  ), call. = FALSE)
}
if (mapping_coverage < 0.8) {
  warning(sprintf("only %.1f%% of finite-p-value tested peaks map to Entrez IDs",
                  100 * mapping_coverage))
}
ann$.gene <- ifelse(entrez_like, source_gene, NA_character_)
ann <- ann[ann$.peak %in% tested_peaks, , drop = FALSE]
universe <- unique(ann$.gene[!is.na(ann$.gene) & nzchar(ann$.gene)])
if (length(universe) == 0L) {
  stop("no peak mapped to a gene id; the peaks and the annotation disagree.", call. = FALSE)
}
message(sprintf(paste0(
  "universe: %d Entrez genes from %d finite-p-value tested peaks ",
  "(%.1f%% peak-to-ID coverage via %s)"),
  length(universe), length(tested_peaks), 100 * mapping_coverage, gene_col[1]
))

# An empty up/down BED is a real DESeq2 result (nothing passed the FDR), so it is
# allowed -- but a peak set that maps to no gene at all is not.
gene_of <- function(bed, direction) {
  if (!file.exists(bed) || file.info(bed)$size == 0) return(character(0))
  p <- unique(as.character(read.delim(bed, header = FALSE)[[4]]))
  valid <- !is.na(ann$.gene) & nzchar(ann$.gene)
  mapped_peaks <- unique(as.character(ann$.peak[ann$.peak %in% p & valid]))
  peak_mapping_coverage <- length(intersect(p, mapped_peaks)) / length(p)
  if (length(mapped_peaks) == 0L) {
    stop(sprintf("none of the %d peaks in %s matched a peak name in the annotation.",
                 length(p), bed), call. = FALSE)
  }
  if (peak_mapping_coverage < 0.5) {
    stop(sprintf(
      paste0("only %.1f%% of %s differential peaks have a valid Entrez mapping; ",
             "refusing biased enrichment"),
      100 * peak_mapping_coverage, direction
    ), call. = FALSE)
  }
  if (peak_mapping_coverage < 0.8) {
    warning(sprintf("only %.1f%% of %s differential peaks have a valid Entrez mapping",
                    100 * peak_mapping_coverage, direction))
  }
  message(sprintf(
    "%s target mapping: %d/%d unique peaks (%.1f%%) have valid Entrez IDs",
    direction, length(intersect(p, mapped_peaks)), length(p),
    100 * peak_mapping_coverage
  ))
  hit <- unique(ann$.gene[ann$.peak %in% p & valid])
  hit
}

sets <- list(up = gene_of(opt$up, "up"), down = gene_of(opt$down, "down"))

run_go <- function(genes, ont) {
  if (!have_org || length(genes) < 5) return(NULL)
  enrichGO(gene = genes, OrgDb = orgdb, keyType = "ENTREZID",
           ont = ont, universe = universe, qvalueCutoff = opt$qvalue,
           readable = TRUE)
}
run_kegg <- function(genes) {
  if (tolower(opt$kegg) != "true" || is.na(kegg_org) || length(genes) < 5) return(NULL)
  # Use the same measured-gene universe as GO. Testing against every gene in
  # KEGG would overstate pathways whose genes are simply more likely to be
  # assigned an accessible peak.
  enrichKEGG(gene = genes, organism = kegg_org, universe = universe,
             qvalueCutoff = opt$qvalue)
}

# ---- collect, write tables, plot --------------------------------------------
panels <- list()
statuses <- data.frame(analysis = character(), status = character(),
                       terms = integer(), stringsAsFactors = FALSE)
add_panel <- function(res, title, skipped = "skipped_fewer_than_5_genes") {
  if (is.null(res)) {
    statuses <<- rbind(statuses, data.frame(analysis = title, status = skipped, terms = 0L))
    return(invisible())
  }
  df <- as.data.frame(res)
  state <- if (nrow(df) == 0L) "no_hits" else "success"
  statuses <<- rbind(statuses, data.frame(analysis = title, status = state,
                                         terms = nrow(df)))
  if (nrow(df) == 0L) return(invisible())
  write.table(df, file.path(opt$out_dir, paste0(gsub("[^A-Za-z0-9]+", "_", title),
              ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  # p.adjust is a magnitude, so it gets the sequential ramp: dark = more
  # significant. enrichplot's default is a red-to-blue gradient, which reads as a
  # diverging scale with no meaningful midpoint -- and red and blue already mean
  # "up" and "down" on the volcano, so a red dot in a panel titled "down peaks"
  # said the opposite of what it meant.
  #
  # dotplot.enrichResult maps adjusted p-values to the fill aesthetic.
  p <- dotplot(res, showCategory = opt$top_n) +
    ggtitle(title) + theme_atac(base_size = 10) +
    scale_fill_gradientn(colours = rev(ATAC_SEQ(256)), name = "adjusted p")
  panels[[length(panels) + 1]] <<- p
}

for (dir in names(sets)) {
  for (ont in onts) add_panel(run_go(sets[[dir]], ont),
                              sprintf("%s peaks — GO:%s", dir, ont))
  kegg_skip <- if (tolower(opt$kegg) == "true")
    "skipped_fewer_than_5_genes" else "disabled"
  add_panel(run_kegg(sets[[dir]]), sprintf("%s peaks — KEGG", dir), kegg_skip)
}

write.table(statuses, opt$out_status, sep = "\t", quote = FALSE, row.names = FALSE)

.atac_pdf(opt$out_pdf, width = 8.5, height = 7, onefile = TRUE)
if (length(panels) == 0) {
  print(ggplot() + annotate("text", 0, 0,
        label = "No enriched terms passed the configured threshold",
        colour = "#52514e", size = 5) + theme_void())
} else {
  for (p in panels) print(p)
}
grDevices::dev.off()
cat(sprintf("Functional enrichment: %d dot-plot panel(s) written.\n", length(panels)))
