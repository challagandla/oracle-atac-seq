#!/usr/bin/env Rscript
# Count one paired-end ATAC-seq library with Rsubread::featureCounts and write
# the standard featureCounts table/summary layout consumed by the checked merge.

suppressPackageStartupMessages(library(Rsubread))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(paste(
    "usage: run_featurecounts.R SAMPLE SAF BAM THREADS TMPDIR",
    "OUT_COUNTS OUT_SUMMARY"
  ), call. = FALSE)
}
sample <- args[[1]]
saf_path <- args[[2]]
bam_path <- args[[3]]
threads <- suppressWarnings(as.integer(args[[4]]))
tmp_dir <- args[[5]]
out_counts <- args[[6]]
out_summary <- args[[7]]

if (!nzchar(sample)) stop("sample identifier is empty", call. = FALSE)
if (is.na(threads) || threads < 1L) stop("threads must be a positive integer", call. = FALSE)
if (!file.exists(saf_path) || !file.exists(bam_path)) {
  stop("SAF annotation and BAM input must both exist", call. = FALSE)
}
if (!dir.exists(tmp_dir)) stop("temporary directory does not exist", call. = FALSE)

saf <- read.delim(
  saf_path,
  colClasses = "character",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
expected_saf <- c("GeneID", "Chr", "Start", "End", "Strand")
if (!identical(colnames(saf), expected_saf) || nrow(saf) == 0L) {
  stop("SAF must contain GeneID, Chr, Start, End, Strand and at least one row",
       call. = FALSE)
}
if (anyDuplicated(saf$GeneID) || any(!nzchar(saf$GeneID)) || any(!nzchar(saf$Chr))) {
  stop("SAF feature identifiers must be non-empty and unique", call. = FALSE)
}
if (any(!grepl("^[0-9]+$", saf$Start)) || any(!grepl("^[0-9]+$", saf$End))) {
  stop("SAF coordinates must be positive integers", call. = FALSE)
}
start <- as.integer(saf$Start)
end <- as.integer(saf$End)
if (anyNA(start) || anyNA(end) || any(start < 1L) || any(end < start)) {
  stop("SAF coordinates are outside the valid 1-based inclusive range", call. = FALSE)
}
if (any(!saf$Strand %in% c("+", "-", "."))) {
  stop("SAF strand must be +, -, or .", call. = FALSE)
}
annotation <- data.frame(
  GeneID = saf$GeneID,
  Chr = saf$Chr,
  Start = start,
  End = end,
  Strand = saf$Strand,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# These arguments are the Rsubread equivalents of CLI featureCounts
# -p --countReadPairs -B -C -F SAF. Multi-mapping and multi-overlap assignment
# retain the featureCounts defaults. The filtered BAM contract has already
# removed secondary, supplementary, duplicate, QC-failed, and improper records.
fc <- featureCounts(
  files = bam_path,
  annot.ext = annotation,
  isGTFAnnotationFile = FALSE,
  useMetaFeatures = TRUE,
  allowMultiOverlap = FALSE,
  countMultiMappingReads = FALSE,
  isPairedEnd = TRUE,
  countReadPairs = TRUE,
  requireBothEndsMapped = TRUE,
  checkFragLength = FALSE,
  countChimericFragments = FALSE,
  autosort = TRUE,
  nthreads = threads,
  tmpDir = tmp_dir,
  verbose = TRUE
)

expected_annotation <- data.frame(
  GeneID = annotation$GeneID,
  Chr = annotation$Chr,
  Start = annotation$Start,
  End = annotation$End,
  Strand = annotation$Strand,
  Length = annotation$End - annotation$Start + 1L,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!identical(fc$annotation, expected_annotation)) {
  stop("Rsubread returned feature metadata that differs from the SAF input",
       call. = FALSE)
}
if (!is.matrix(fc$counts) || nrow(fc$counts) != nrow(annotation) ||
    ncol(fc$counts) != 1L || any(!is.finite(fc$counts)) ||
    any(fc$counts < 0) || any(fc$counts != floor(fc$counts))) {
  stop("Rsubread returned an invalid single-library count matrix", call. = FALSE)
}
if (!is.data.frame(fc$stat) || ncol(fc$stat) != 2L ||
    colnames(fc$stat)[[1]] != "Status" || anyDuplicated(fc$stat$Status) ||
    any(!nzchar(fc$stat$Status))) {
  stop("Rsubread returned an invalid assignment summary", call. = FALSE)
}
stat_values <- suppressWarnings(as.numeric(fc$stat[[2]]))
if (anyNA(stat_values) || any(!is.finite(stat_values)) || any(stat_values < 0) ||
    any(stat_values != floor(stat_values))) {
  stop("Rsubread summary values must be non-negative integers", call. = FALSE)
}
assigned_index <- match("Assigned", fc$stat$Status)
if (is.na(assigned_index) || stat_values[[assigned_index]] <= 0) {
  stop("Rsubread assigned zero fragments to the consensus peak universe",
       call. = FALSE)
}
if (sum(fc$counts[, 1]) != stat_values[[assigned_index]]) {
  stop("count-column sum differs from the Rsubread Assigned summary", call. = FALSE)
}

count_table <- data.frame(
  Geneid = fc$annotation$GeneID,
  Chr = fc$annotation$Chr,
  Start = fc$annotation$Start,
  End = fc$annotation$End,
  Strand = fc$annotation$Strand,
  Length = fc$annotation$Length,
  count = sprintf("%.0f", fc$counts[, 1]),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
colnames(count_table)[[7]] <- sample
summary_table <- data.frame(
  Status = fc$stat$Status,
  value = sprintf("%.0f", stat_values),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
colnames(summary_table)[[2]] <- sample

dir.create(dirname(out_counts), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_summary), recursive = TRUE, showWarnings = FALSE)
counts_tmp <- tempfile(pattern = paste0(".", basename(out_counts), "."),
                       tmpdir = dirname(out_counts))
summary_tmp <- tempfile(pattern = paste0(".", basename(out_summary), "."),
                        tmpdir = dirname(out_summary))
on.exit(unlink(c(counts_tmp, summary_tmp)), add = TRUE)
counts_connection <- file(counts_tmp, open = "wt")
tryCatch({
  writeLines(
    sprintf("# Program:Rsubread featureCounts v%s; sample=%s",
            as.character(packageVersion("Rsubread")), sample),
    counts_connection
  )
  write.table(count_table, counts_connection, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
}, finally = close(counts_connection))
write.table(summary_table, summary_tmp, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
if (!file.rename(counts_tmp, out_counts) || !file.rename(summary_tmp, out_summary)) {
  stop("could not atomically publish featureCounts outputs", call. = FALSE)
}
