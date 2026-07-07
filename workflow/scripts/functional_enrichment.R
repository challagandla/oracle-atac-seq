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
  make_option("--up"), make_option("--down"),        # up/down peak BEDs
  make_option("--orgdb", default = ""),
  make_option("--taxid", type = "integer", default = 9606),
  make_option("--ontologies", default = "BP,MF"),
  make_option("--kegg", default = "true"),
  make_option("--qvalue", type = "double", default = 0.05),
  make_option("--top-n", type = "integer", default = 20, dest = "top_n"),
  make_option("--out-pdf", dest = "out_pdf"),
  make_option("--out-dir", dest = "out_dir")
)))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
kegg_org <- c("9606" = "hsa", "10090" = "mmu", "10116" = "rno")[as.character(opt$taxid)]
onts <- trimws(strsplit(opt$ontologies, ",")[[1]])

have_org <- nzchar(opt$orgdb) && requireNamespace(opt$orgdb, quietly = TRUE)
if (have_org) suppressPackageStartupMessages(library(opt$orgdb, character.only = TRUE))

# ---- map peaks -> Entrez gene ids via the annotation ------------------------
ann <- tryCatch(read.delim(opt$annotation, stringsAsFactors = FALSE),
                error = function(e) data.frame())
peak_col <- intersect(c("V4", "name", "peak"), colnames(ann))
gene_col <- intersect(c("geneId", "ENTREZID", "geneID"), colnames(ann))
gene_of <- function(bed) character(0)
universe <- character(0)
if (length(peak_col) && length(gene_col) && nrow(ann) > 0) {
  ann$.peak <- ann[[peak_col[1]]]; ann$.gene <- as.character(ann[[gene_col[1]]])
  universe <- unique(ann$.gene[!is.na(ann$.gene) & nzchar(ann$.gene)])
  gene_of <- function(bed) {
    if (!file.exists(bed) || file.info(bed)$size == 0) return(character(0))
    p <- read.delim(bed, header = FALSE)[[4]]
    unique(ann$.gene[ann$.peak %in% p & !is.na(ann$.gene) & nzchar(ann$.gene)])
  }
}

sets <- list(up = gene_of(opt$up), down = gene_of(opt$down))

run_go <- function(genes, ont) {
  if (!have_org || length(genes) < 5) return(NULL)
  tryCatch(enrichGO(gene = genes, OrgDb = get(opt$orgdb), keyType = "ENTREZID",
                    ont = ont, universe = universe, qvalueCutoff = opt$qvalue,
                    readable = TRUE),
           error = function(e) { message("enrichGO(", ont, ") failed: ",
                                         conditionMessage(e)); NULL })
}
run_kegg <- function(genes) {
  if (tolower(opt$kegg) != "true" || is.na(kegg_org) || length(genes) < 5) return(NULL)
  tryCatch(enrichKEGG(gene = genes, organism = kegg_org, qvalueCutoff = opt$qvalue),
           error = function(e) { message("enrichKEGG failed (needs network): ",
                                         conditionMessage(e)); NULL })
}

# ---- collect, write tables, plot --------------------------------------------
panels <- list()
add_panel <- function(res, title) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(invisible())
  df <- as.data.frame(res)
  write.table(df, file.path(opt$out_dir, paste0(gsub("[^A-Za-z0-9]+", "_", title),
              ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  p <- dotplot(res, showCategory = opt$top_n) +
    ggtitle(title) + theme_atac(base_size = 10) +
    scale_colour_gradient(low = "#e34948", high = "#2a78d6")
  panels[[length(panels) + 1]] <<- p
}

for (dir in names(sets)) {
  for (ont in onts) add_panel(run_go(sets[[dir]], ont),
                              sprintf("%s peaks — GO:%s", dir, ont))
  add_panel(run_kegg(sets[[dir]]), sprintf("%s peaks — KEGG", dir))
}

.atac_pdf(opt$out_pdf, width = 8.5, height = 7, onefile = TRUE)
if (length(panels) == 0) {
  print(ggplot() + annotate("text", 0, 0,
        label = "No enriched terms (or OrgDb/gene sets unavailable)",
        colour = "#52514e", size = 5) + theme_void())
} else {
  for (p in panels) print(p)
}
grDevices::dev.off()
cat(sprintf("Functional enrichment: %d dot-plot panel(s) written.\n", length(panels)))
