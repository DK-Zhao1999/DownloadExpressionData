#!/usr/bin/env Rscript
# Remove normal / adjacent / user-specified samples from a prepared RData.
#
# This is intentionally a separate step.  The agent must ask the user first;
# there is no implicit "drop normal" default anywhere in the pipeline.
#
# Example:
#   Rscript subset_expression.R \
#     --rdata GSE50760/GSE50760_expression.RData \
#     --outdir GSE50760/tumor_only \
#     --column tissue_type \
#     --drop-regex "normal|adjacent|non-tumor|healthy|benign" \
#     --label tumor_only
#
# Use --dry-run first to show exactly which samples would be removed.

options(stringsAsFactors = FALSE, warn = 1)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) stop(sprintf("unexpected argument: %s", token), call. = FALSE)
    key <- sub("^--", "", token)
    if (grepl("=", key, fixed = TRUE)) {
      parts <- strsplit(key, "=", fixed = TRUE)[[1]]
      out[[parts[[1]]]] <- paste(parts[-1], collapse = "=")
      i <- i + 1L
    } else if (i < length(args) && !startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    } else {
      out[[key]] <- TRUE
      i <- i + 1L
    }
  }
  out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
arg <- function(name, default = NULL) if (!is.null(opt[[name]])) opt[[name]] else default
flag <- function(name) isTRUE(opt[[name]])
require_arg <- function(name) {
  value <- arg(name)
  if (is.null(value) || identical(value, TRUE) || !nzchar(as.character(value))) {
    stop(sprintf("missing required --%s", name), call. = FALSE)
  }
  as.character(value)
}
stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

rdata_path <- normalizePath(require_arg("rdata"), mustWork = TRUE)
outdir <- require_arg("outdir")
label <- gsub("[^A-Za-z0-9_.-]+", "_", arg("label", "tumor_only"))
if (!nzchar(label)) label <- "tumor_only"
force <- flag("force")
dry_run <- flag("dry-run")
min_samples <- suppressWarnings(as.integer(arg("min-samples", 2)))
if (is.na(min_samples) || min_samples < 1L) min_samples <- 2L

env <- new.env(parent = baseenv())
load(rdata_path, envir = env)
if (is.null(env$rt) || !(is.data.frame(env$rt) || is.matrix(env$rt))) {
  stopf("RData does not contain a data.frame/matrix 'rt'")
}
if (is.null(env$metadata) || !is.data.frame(env$metadata)) stopf("RData does not contain data.frame 'metadata'")
rt <- as.data.frame(env$rt, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- env$metadata
dataset_info <- env$dataset_info
if (is.null(dataset_info)) dataset_info <- list()

if (!"sample_id" %in% colnames(metadata)) stopf("metadata has no sample_id column")
if (!identical(colnames(rt), as.character(metadata$sample_id))) {
  stopf("metadata sample_id does not match expression columns; refuse to subset")
}
if (!all(vapply(rt, is.numeric, logical(1)))) stopf("expression data.frame has non-numeric columns")
if (anyNA(rt)) stopf("expression data.frame contains NA; run prepare_expression.R first")
if (anyDuplicated(rownames(rt))) stopf("expression data.frame has duplicate gene symbols")
if (any(grepl("[[:space:]]", rownames(rt)))) stopf("gene symbols contain whitespace")

n_before <- ncol(rt)
removed <- character()
selected <- rep(TRUE, n_before)

keep_file <- arg("keep-file", NULL)
column <- arg("column", NULL)
drop_regex <- arg("drop-regex", NULL)
keep_regex <- arg("keep-regex", NULL)

if (!is.null(keep_file) && nzchar(as.character(keep_file))) {
  ids <- trimws(readLines(normalizePath(as.character(keep_file), mustWork = TRUE), warn = FALSE))
  ids <- ids[nzchar(ids) & ids != "sample_id"]
  if (!length(ids)) stopf("--keep-file contains no sample IDs")
  selected <- colnames(rt) %in% ids
  if (!any(selected)) stopf("no samples from --keep-file are present in the data.frame")
  removed <- colnames(rt)[!selected]
  message(sprintf("[subset] keep-file selected %d / %d samples", sum(selected), n_before))
} else {
  if (is.null(column) || !nzchar(as.character(column))) {
    stopf("provide --column with --drop-regex/--keep-regex, or provide --keep-file")
  }
  if (!column %in% colnames(metadata)) stopf("metadata has no column: %s", column)
  values <- as.character(metadata[[column]])
  message(sprintf("[subset] values in metadata$%s:", column))
  print(table(values, useNA = "ifany"))

  if (!is.null(keep_regex) && nzchar(as.character(keep_regex))) {
    selected <- grepl(keep_regex, values, ignore.case = TRUE, perl = TRUE)
    message(sprintf("[subset] keep-regex matched %d sample(s)", sum(selected)))
  }
  if (!is.null(drop_regex) && nzchar(as.character(drop_regex))) {
    drop <- grepl(drop_regex, values, ignore.case = TRUE, perl = TRUE)
    message(sprintf("[subset] drop-regex matched %d sample(s)", sum(drop)))
    selected <- selected & !drop
  }
  if (is.null(keep_regex) && is.null(drop_regex)) {
    stopf("provide at least one of --keep-regex or --drop-regex")
  }
  removed <- colnames(rt)[!selected]
}

if (!any(selected)) stopf("no samples remain after filtering; refusing to write an empty data.frame")
if (sum(selected) < min_samples) {
  stopf("%d sample(s) remain, below --min-samples=%d", sum(selected), min_samples)
}
if (sum(selected) == n_before) {
  stopf("the filter removed no samples; check --column and the regex")
}

message(sprintf("[subset] removing %d sample(s): %s", length(removed), paste(removed, collapse = ", ")))
message(sprintf("[subset] keeping %d sample(s)", sum(selected)))

if (dry_run) {
  message("[subset] dry run; no files written")
  quit(save = "no", status = 0)
}

rt <- rt[, selected, drop = FALSE]
metadata <- metadata[selected, , drop = FALSE]
rownames(metadata) <- NULL
if (!identical(colnames(rt), as.character(metadata$sample_id))) {
  stopf("internal error: subset metadata does not match data.frame columns")
}

dataset_info$n_samples <- ncol(rt)
dataset_info$transformations <- c(dataset_info$transformations, paste0("subset_", label))
dataset_info$subset <- list(
  label = label,
  column = if (!is.null(column)) column else NA_character_,
  keep_regex = if (!is.null(keep_regex)) keep_regex else NA_character_,
  drop_regex = if (!is.null(drop_regex)) drop_regex else NA_character_,
  keep_file = if (!is.null(keep_file)) normalizePath(as.character(keep_file), mustWork = FALSE) else NA_character_,
  n_before = n_before,
  n_after = ncol(rt),
  removed_samples = removed,
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
prefix <- dataset_info$accession
if (is.null(prefix) || is.na(prefix) || !nzchar(prefix)) {
  prefix <- tools::file_path_sans_ext(basename(rdata_path))
}
prefix <- gsub("[^A-Za-z0-9_.-]+", "_", prefix)
out_rdata <- file.path(outdir, paste0(prefix, "_", label, ".RData"))
if (file.exists(out_rdata) && !force) {
  stopf("output exists: %s (use --force to overwrite)", out_rdata)
}
save(rt, metadata, dataset_info, file = out_rdata, compress = "gzip")

removed_table <- env$metadata[!selected, , drop = FALSE]
write.table(
  removed_table,
  file = file.path(outdir, paste0(prefix, "_", label, "_removed_samples.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
summary_table <- as.data.frame(table(metadata[[if (!is.null(column)) column else "sample_id"]]))
colnames(summary_table) <- c("value", "n_after")
write.table(
  summary_table,
  file = file.path(outdir, paste0(prefix, "_", label, "_summary.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

message(sprintf("[subset] saved: %s", out_rdata))
message(sprintf("[subset] removed-samples table: %s", file.path(outdir, paste0(prefix, "_", label, "_removed_samples.tsv"))))
message("[subset] next: ask whether clinical/prognosis data are needed.")

