#!/usr/bin/env Rscript
# Prepare a bulk gene-expression matrix for downstream analysis.
#
# Contract:
#   input  : genes x samples (or samples x genes with --transpose yes)
#   output : <prefix>_expression.RData containing
#              rt           data.frame, rows = clean gene symbols,
#                           columns = sample names, no NA/Inf
#              metadata     data.frame, one row per sample, sample_id first
#              dataset_info list describing every transformation
#
# The processing order is deliberate and follows the analysis SOP:
#   1. harmonize IDs; remove whitespace; keep first symbol before "///";
#      merge duplicate gene symbols by mean
#   2. KNN-impute missing values with impute::impute.knn
#   3. if min(rt) < 0, shift rt <- rt - min(rt)
#   4. if range is not inside [0, 30], apply log2(rt + 1)
#   5. save RData + QC + processing log
#
# This script does NOT remove normal samples. That is a separate, explicit
# decision performed by subset_expression.R after the user is asked.

options(stringsAsFactors = FALSE, warn = 1)

# ---------------------------------------------------------------------------
# CLI helpers
# ---------------------------------------------------------------------------

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) {
      stop(sprintf("unexpected argument: %s", token), call. = FALSE)
    }
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

arg <- function(name, default = NULL) {
  if (!is.null(opt[[name]])) opt[[name]] else default
}
flag <- function(name) isTRUE(opt[[name]])
require_arg <- function(name) {
  value <- arg(name)
  if (is.null(value) || identical(value, TRUE) || !nzchar(as.character(value))) {
    stop(sprintf("missing required --%s", name), call. = FALSE)
  }
  as.character(value)
}
stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

as_int <- function(x, default) {
  if (is.null(x) || !nzchar(as.character(x))) return(as.integer(default))
  value <- suppressWarnings(as.integer(x))
  if (is.na(value)) stopf("expected an integer, got: %s", x)
  value
}
as_num <- function(x, default) {
  if (is.null(x) || !nzchar(as.character(x))) return(as.numeric(default))
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value)) stopf("expected a number, got: %s", x)
  value
}

input_path <- normalizePath(require_arg("input"), mustWork = TRUE)
outdir <- require_arg("outdir")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
accession <- arg("accession", "")
prefix <- arg("prefix", if (nzchar(accession)) accession else tools::file_path_sans_ext(basename(input_path)))
prefix <- gsub("[^A-Za-z0-9_.-]+", "_", prefix)
if (!nzchar(prefix)) prefix <- "expression"

value_type <- tolower(arg("value-type", "unknown"))
valid_types <- c(
  "tpm", "fpkm", "rpkm", "counts", "array_linear", "array_log2",
  "normalized", "unknown"
)
if (!value_type %in% valid_types) {
  stopf("--value-type must be one of: %s", paste(valid_types, collapse = ", "))
}
if (value_type == "counts" && !flag("allow-counts")) {
  stop(
    paste0(
      "--value-type counts requires --allow-counts.\n",
      "Raw counts are not TPM.  Prefer an author-provided TPM matrix, ",
      "recount3/ARCHS4, or quantify raw reads first.  If a count-based ",
      "analysis is truly intended, pass --allow-counts and report that ",
      "log2(counts + 1) is not TPM."
    ),
    call. = FALSE
  )
}

log_lines <- character()
log_step <- function(fmt, ...) {
  msg <- sprintf(fmt, ...)
  log_lines <<- c(log_lines, sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg))
  message(msg)
}

# ---------------------------------------------------------------------------
# Input readers
# ---------------------------------------------------------------------------

read_table_auto <- function(path, sep = NULL, nrows = -1L) {
  lower <- tolower(path)
  if (grepl("\\.xlsx?$", lower)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stopf("reading %s needs the openxlsx package", basename(path))
    }
    return(openxlsx::read.xlsx(path, sheet = 1, check.names = FALSE))
  }
  skiprows <- 0L
  if (is.null(sep)) {
    con <- if (grepl("\\.gz$", lower)) gzfile(path, "rt") else file(path, "rt")
    raw_lines <- tryCatch(readLines(con, n = 500L, warn = FALSE), finally = close(con))
    if (length(raw_lines)) raw_lines[[1]] <- sub("^\ufeff", "", raw_lines[[1]])
    header_index <- NA_integer_
    for (index in seq_along(raw_lines)) {
      line <- trimws(raw_lines[[index]])
      if (!nzchar(line) || grepl("^[!#^]", line)) next
      header_index <- index
      break
    }
    if (is.na(header_index)) stopf("no tabular content found in %s", basename(path))
    header_line <- raw_lines[[header_index]]
    counts <- c(
      tab = lengths(regmatches(header_line, gregexpr("\t", header_line, fixed = TRUE))),
      comma = lengths(regmatches(header_line, gregexpr(",", header_line, fixed = TRUE))),
      semi = lengths(regmatches(header_line, gregexpr(";", header_line, fixed = TRUE)))
    )
    sep <- c("\t", ",", ";")[[which.max(counts)]]
    if (max(counts) == 0L) sep <- ""
    skiprows <- header_index - 1L
  }
  out <- read.delim(
    path,
    sep = sep,
    header = TRUE,
    check.names = FALSE,
    quote = "\"",
    comment.char = "",
    na.strings = c("NA", "N/A", "null", "NULL", "NaN", "nan", ""),
    nrows = nrows,
    skip = skiprows,
    stringsAsFactors = FALSE
  )
  if (ncol(out) > 0L) colnames(out)[[1]] <- sub("^\ufeff", "", colnames(out)[[1]])
  out
}

extract_matrix_object <- function(obj) {
  if (is.matrix(obj)) return(obj)
  if (is.data.frame(obj)) return(obj)
  if (is.list(obj)) {
    for (name in c("rt", "expr", "expression", "matrix", "data")) {
      if (!is.null(obj[[name]]) && (is.matrix(obj[[name]]) || is.data.frame(obj[[name]]))) {
        return(obj[[name]])
      }
    }
  }
  if (requireNamespace("SummarizedExperiment", quietly = TRUE) &&
      methods::is(obj, "SummarizedExperiment")) {
    return(SummarizedExperiment::assay(obj))
  }
  NULL
}

read_input_object <- function(path) {
  lower <- tolower(path)
  if (grepl("\\.rds$", lower)) {
    obj <- readRDS(path)
    mat <- extract_matrix_object(obj)
    if (is.null(mat)) stopf("no matrix/data.frame found in RDS: %s", basename(path))
    return(list(data = mat, metadata = NULL))
  }
  if (grepl("\\.(rdata|rda)$", lower)) {
    env <- new.env(parent = baseenv())
    load(path, envir = env)
    mat <- NULL
    for (name in c("rt", "expr", "expression", "matrix", "data")) {
      if (!is.null(env[[name]])) {
        mat <- extract_matrix_object(env[[name]])
        if (!is.null(mat)) break
      }
    }
    if (is.null(mat)) {
      for (name in ls(env)) {
        mat <- extract_matrix_object(env[[name]])
        if (!is.null(mat)) break
      }
    }
    if (is.null(mat)) stopf("no matrix/data.frame found in RData: %s", basename(path))
    metadata <- NULL
    if (!is.null(env$metadata) && is.data.frame(env$metadata)) metadata <- env$metadata
    dataset_info <- if (!is.null(env$dataset_info)) env$dataset_info else NULL
    return(list(data = mat, metadata = metadata, dataset_info = dataset_info))
  }
  list(data = read_table_auto(path), metadata = NULL)
}

resolve_column <- function(frame, spec, default = NULL) {
  if (is.null(spec) || !nzchar(as.character(spec))) {
    if (is.null(default)) stopf("no column specified")
    return(default)
  }
  spec <- as.character(spec)
  if (spec %in% colnames(frame)) return(spec)
  index <- suppressWarnings(as.integer(spec))
  if (!is.na(index) && index >= 1L && index <= ncol(frame)) return(colnames(frame)[[index]])
  stopf("column not found: %s", spec)
}

read_metadata <- function(path, sample_column = NULL) {
  frame <- read_table_auto(path)
  if (nrow(frame) == 0L) stopf("metadata file is empty: %s", basename(path))
  candidates <- c("sample_id", "sample", "SampleID", "geo_accession", "GSM", "title", "sample_name")
  chosen <- NULL
  if (!is.null(sample_column) && nzchar(as.character(sample_column))) {
    chosen <- resolve_column(frame, sample_column)
  } else {
    hit <- candidates[candidates %in% colnames(frame)]
    if (length(hit)) chosen <- hit[[1]]
  }
  if (is.null(chosen)) {
    warning("no sample-ID column recognized; using the first metadata column")
    chosen <- colnames(frame)[[1]]
  }
  ids <- trimws(as.character(frame[[chosen]]))
  if (any(!nzchar(ids))) stopf("metadata contains empty sample IDs in %s", chosen)
  if (anyDuplicated(ids)) stopf("metadata contains duplicate sample IDs in %s", chosen)
  frame[[chosen]] <- ids
  if (chosen != "sample_id") colnames(frame)[colnames(frame) == chosen] <- "sample_id"
  frame <- frame[, c("sample_id", setdiff(colnames(frame), "sample_id")), drop = FALSE]
  rownames(frame) <- NULL
  frame
}

# ---------------------------------------------------------------------------
# Orientation and gene IDs
# ---------------------------------------------------------------------------

looks_like_sample_id <- function(x) {
  values <- as.character(x)
  mean(grepl("^(GSM\\d+|TCGA-|E-MTAB-|PAT\\d+|S\\d+$|Sample)", values)) > 0.8
}

looks_like_id_column <- function(frame, column = 1L) {
  values <- frame[[column]]
  if (is.character(values) || is.factor(values)) return(TRUE)
  name <- colnames(frame)[[column]]
  grepl("(^|_)(id|gene|probe|symbol|entrez)(_|$)", tolower(name)) ||
    grepl("ID_REF|Probe", name, ignore.case = TRUE)
}

clean_gene_symbols <- function(x) {
  # Final object contract: one row per gene symbol, no whitespace, and for
  # composite GEO symbols such as "A /// B" keep only the first symbol.
  x <- as.character(x)
  x <- sub("///.*$", "", x)
  x <- gsub("[[:space:]\u00A0\u3000]+", "", x, perl = TRUE)
  trimws(x)
}

fix_excel_date_genes <- function(x) {
  # Excel turns some valid gene symbols into dates.  This is a fixed,
  # deliberately narrow map: only the MARCH, SEPT and DEC families.
  #   MARCH1 -> 1-Mar, ..., MARCH11 -> 11-Mar
  #   SEPT1  -> 1-Sep, ..., SEPT15  -> 15-Sep
  #   DEC1   -> 1-Dec, DEC2   -> 2-Dec
  # Other date-like strings are NOT guessed; they are reported as unresolved.
  map <- c(
    "1-Mar" = "MARCH1", "2-Mar" = "MARCH2", "3-Mar" = "MARCH3",
    "4-Mar" = "MARCH4", "5-Mar" = "MARCH5", "6-Mar" = "MARCH6",
    "7-Mar" = "MARCH7", "8-Mar" = "MARCH8", "9-Mar" = "MARCH9",
    "10-Mar" = "MARCH10", "11-Mar" = "MARCH11",
    "1-Sep" = "SEPT1", "2-Sep" = "SEPT2", "3-Sep" = "SEPT3",
    "4-Sep" = "SEPT4", "5-Sep" = "SEPT5", "6-Sep" = "SEPT6",
    "7-Sep" = "SEPT7", "8-Sep" = "SEPT8", "9-Sep" = "SEPT9",
    "10-Sep" = "SEPT10", "11-Sep" = "SEPT11", "12-Sep" = "SEPT12",
    "13-Sep" = "SEPT13", "14-Sep" = "SEPT14", "15-Sep" = "SEPT15",
    "1-Dec" = "DEC1", "2-Dec" = "DEC2"
  )
  values <- trimws(as.character(x))
  keys <- sub("^0([1-9])-", "\\1-", values)
  index <- match(tolower(keys), tolower(names(map)))
  hit <- !is.na(index)
  values[hit] <- unname(map[index[hit]])
  list(values = values, n_fixed = sum(hit))
}

orient_input <- function(frame, metadata, transpose = "auto", gene_column = NULL) {
  frame <- as.data.frame(frame, check.names = FALSE, stringsAsFactors = FALSE)
  if (inherits(frame[[1]], c("Date", "POSIXct", "POSIXlt"))) {
    stop(
      paste0(
        "the first column was read as a Date, so Excel date-corrupted gene ",
        "symbols cannot be recovered safely. Re-export the gene column as text."
      ),
      call. = FALSE
    )
  }
  first_is_id <- looks_like_id_column(frame, 1L)
  data_columns <- if (first_is_id) colnames(frame)[-1] else colnames(frame)
  candidate_rows <- if (first_is_id) as.character(frame[[1]]) else rownames(frame)

  do_transpose <- FALSE
  if (identical(transpose, "yes")) {
    do_transpose <- TRUE
  } else if (identical(transpose, "auto")) {
    if (!is.null(metadata) && nrow(metadata) > 0L) {
      column_match <- mean(data_columns %in% metadata$sample_id)
      row_match <- mean(candidate_rows %in% metadata$sample_id)
      if (row_match >= 0.8 && column_match < 0.5) do_transpose <- TRUE
    } else {
      if (looks_like_sample_id(candidate_rows) && !looks_like_sample_id(data_columns)) {
        do_transpose <- TRUE
      }
    }
  } else if (!identical(transpose, "no")) {
    stopf("--transpose must be auto, yes, or no")
  }

  if (do_transpose) {
    sample_ids <- if (first_is_id) as.character(frame[[1]]) else rownames(frame)
    if (is.null(sample_ids) || any(!nzchar(sample_ids))) {
      stopf("cannot find sample IDs in the first column/row names for transposition")
    }
    if (anyDuplicated(sample_ids)) stopf("duplicated sample IDs after transposition")
    gene_ids <- if (first_is_id) colnames(frame)[-1] else colnames(frame)
    mat <- as.matrix(frame[, if (first_is_id) -1L else TRUE, drop = FALSE])
    rt <- t(mat)
    rownames(rt) <- gene_ids
    colnames(rt) <- sample_ids
    log_step("transposed input: %d genes x %d samples", nrow(rt), ncol(rt))
    return(list(rt = rt, gene_ids = gene_ids))
  }

  if (!is.null(gene_column) && nzchar(as.character(gene_column))) {
    gene_col <- resolve_column(frame, gene_column)
    if (inherits(frame[[gene_col]], c("Date", "POSIXct", "POSIXlt"))) {
      stopf("gene column %s was read as a Date; re-export it as text", gene_col)
    }
    gene_ids <- as.character(frame[[gene_col]])
    keep <- setdiff(colnames(frame), gene_col)
    mat <- as.matrix(frame[, keep, drop = FALSE])
    log_step("gene IDs taken from column: %s", gene_col)
  } else if (first_is_id) {
    gene_ids <- as.character(frame[[1]])
    mat <- as.matrix(frame[, -1L, drop = FALSE])
  } else {
    gene_ids <- rownames(frame)
    mat <- as.matrix(frame)
    if (is.null(gene_ids) || identical(gene_ids, as.character(seq_len(nrow(frame))))) {
      stop(
        paste0(
          "no gene-ID column or row names found.  Use --gene-column or ",
          "provide a matrix with gene symbols/IDs in the first column."
        ),
        call. = FALSE
      )
    }
  }
  rownames(mat) <- gene_ids
  list(rt = mat, gene_ids = gene_ids)
}

read_annotation <- function(path, probe_column = NULL, gene_column = NULL) {
  frame <- read_table_auto(path)
  if (nrow(frame) == 0L) stopf("annotation file is empty: %s", basename(path))
  if (is.null(probe_column) || !nzchar(as.character(probe_column))) {
    probe_candidates <- c("ID", "ID_REF", "Probe", "probe_id", "Probe ID", "Affy Probe Set Name")
    hit <- probe_candidates[probe_candidates %in% colnames(frame)]
    probe_column <- if (length(hit)) hit[[1]] else colnames(frame)[[1]]
  }
  if (is.null(gene_column) || !nzchar(as.character(gene_column))) {
    gene_candidates <- c(
      "Gene symbol", "GENE_SYMBOL", "gene_symbol", "Gene Symbol",
      "Symbol", "symbol", "Gene", "gene_name", "gene"
    )
    hit <- gene_candidates[gene_candidates %in% colnames(frame)]
    if (!length(hit)) {
      stop(
        paste0(
          "cannot find a gene-symbol column in the annotation.  Pass ",
          "--annotation-gene-col with the exact column name."
        ),
        call. = FALSE
      )
    }
    gene_column <- hit[[1]]
  }
  probe_col <- resolve_column(frame, probe_column)
  gene_col <- resolve_column(frame, gene_column)
  mapping <- data.frame(
    probe = trimws(as.character(frame[[probe_col]])),
    gene = trimws(as.character(frame[[gene_col]])),
    stringsAsFactors = FALSE
  )
  mapping <- mapping[
    nzchar(mapping$probe) & nzchar(mapping$gene) &
      !grepl("^[!#^]", mapping$probe),
    ,
    drop = FALSE
  ]
  mapping <- mapping[!duplicated(mapping$probe), , drop = FALSE]
  log_step(
    "annotation: %d probe-to-gene mappings from columns %s -> %s",
    nrow(mapping), probe_col, gene_col
  )
  mapping
}

# ---------------------------------------------------------------------------
# Matrix preparation
# ---------------------------------------------------------------------------

metadata_arg <- arg("metadata", NULL)
metadata <- NULL
if (!is.null(metadata_arg) && nzchar(as.character(metadata_arg))) {
  metadata <- read_metadata(
    normalizePath(as.character(metadata_arg), mustWork = TRUE),
    sample_column = arg("sample-column", NULL)
  )
}
input <- read_input_object(input_path)
if (!is.null(input$dataset_info) &&
    !is.null(input$dataset_info$output_value_type) &&
    nzchar(as.character(input$dataset_info$output_value_type))) {
  stop(
    paste0(
      "input RData is already a processed expression object (",
      input$dataset_info$output_value_type,
      ").  Do not re-run prepare_expression.R on it: start from the raw ",
      "matrix, or use subset_expression.R to remove samples."
    ),
    call. = FALSE
  )
}
if (is.null(metadata) && !is.null(input$metadata)) {
  metadata <- input$metadata
  if (!"sample_id" %in% colnames(metadata)) {
    metadata <- NULL
  } else {
    metadata$sample_id <- trimws(as.character(metadata$sample_id))
    metadata <- metadata[, c("sample_id", setdiff(colnames(metadata), "sample_id")), drop = FALSE]
  }
}
if (is.null(metadata)) {
  warning("no metadata supplied; metadata will contain sample_id only")
}

oriented <- orient_input(
  input$data,
  metadata,
  transpose = tolower(arg("transpose", "auto")),
  gene_column = arg("gene-column", NULL)
)
rt <- oriented$rt
gene_ids <- oriented$gene_ids

if (nrow(rt) < 2L) stopf("matrix has fewer than 2 genes")
if (ncol(rt) < 2L) stopf("matrix has fewer than 2 samples")
if (is.null(rownames(rt)) || any(!nzchar(trimws(rownames(rt))))) {
  stopf("gene IDs are missing from row names/first column")
}
if (is.null(colnames(rt)) || any(!nzchar(trimws(colnames(rt))))) {
  stopf("sample IDs are missing from column names")
}
if (anyDuplicated(colnames(rt))) stopf("duplicate sample IDs in the matrix")

# Convert to numeric while preserving the original non-numeric count.
numeric_rt <- suppressWarnings(matrix(
  as.numeric(rt),
  nrow = nrow(rt),
  ncol = ncol(rt),
  dimnames = dimnames(rt)
))
nonnumeric_cells <- sum(is.na(numeric_rt) & !is.na(rt))
if (nonnumeric_cells > 0L) {
  log_step("converted %d non-numeric cell(s) to NA", nonnumeric_cells)
  if (nonnumeric_cells > 0.1 * length(rt)) {
    stop(
      paste0(
        "more than 10% of cells are non-numeric; the input may contain ",
        "annotation columns or the wrong separator.  Inspect the file and ",
        "use --gene-column / a clean matrix."
      ),
      call. = FALSE
    )
  }
}
rt <- numeric_rt
mapping_rate <- NA_real_
n_unmapped <- 0L
n_multi_symbol <- 0L
n_duplicates_collapsed <- 0L
collapse_method <- "none"
n_date_genes_fixed <- 0L
n_unresolved_date_like <- 0L
n_all_na_samples_removed <- 0L

bad_finite <- !is.na(rt) & !is.finite(rt)
if (any(bad_finite)) {
  stopf("matrix contains %d Inf/-Inf cell(s); fix the input before imputation", sum(bad_finite))
}

# Annotation-based probe -> gene mapping.
annotation_arg <- arg("annotation", NULL)
if (!is.null(annotation_arg) && nzchar(as.character(annotation_arg))) {
  mapping <- read_annotation(
    normalizePath(as.character(annotation_arg), mustWork = TRUE),
    probe_column = arg("annotation-probe-col", NULL),
    gene_column = arg("annotation-gene-col", NULL)
  )
  mapped <- mapping$gene[match(gene_ids, mapping$probe)]
  mapping_rate <- mean(!is.na(mapped))
  min_rate <- as_num(arg("min-mapping-rate", 0.5), 0.5)
  log_step("probe mapping rate: %.1f%% (minimum %.1f%%)", 100 * mapping_rate, 100 * min_rate)
  if (mapping_rate < min_rate) {
    stopf(
      paste0(
        "probe mapping rate %.1f%% is below --min-mapping-rate %.1f%%. ",
        "Check the annotation platform and probe IDs."
      ),
      100 * mapping_rate, 100 * min_rate
    )
  }
  keep <- !is.na(mapped)
  n_unmapped <- sum(!keep)
  rt <- rt[keep, , drop = FALSE]
  gene_ids <- mapped[keep]
  log_step("dropped %d unmapped probe(s)", n_unmapped)
}

# Clean gene symbols before duplicate handling.  Final contract:
#   - "A /// B" -> "A"
#   - all whitespace removed
#   - no NA/empty symbols
n_multi_symbol <- sum(grepl("///", gene_ids, fixed = TRUE), na.rm = TRUE)
gene_ids <- clean_gene_symbols(gene_ids)
if (n_multi_symbol > 0L) {
  log_step("kept the first symbol for %d composite gene ID(s)", n_multi_symbol)
}

# Restore known Excel date-corrupted gene symbols before duplicate merging.
date_fixed <- fix_excel_date_genes(gene_ids)
gene_ids <- date_fixed$values
n_date_genes_fixed <- date_fixed$n_fixed
if (n_date_genes_fixed > 0L) {
  log_step("restored %d Excel date-like gene symbol(s) to MARCH/SEPT/DEC names", n_date_genes_fixed)
}

keep <- !is.na(gene_ids) & nzchar(gene_ids)
if (!all(keep)) {
  log_step("dropped %d row(s) with empty gene IDs after symbol cleaning", sum(!keep))
  rt <- rt[keep, , drop = FALSE]
  gene_ids <- gene_ids[keep]
}
if (any(grepl("///", gene_ids, fixed = TRUE))) {
  stopf("composite '///' gene symbols remain after cleaning")
}
if (any(grepl("[[:space:]]", gene_ids))) {
  stopf("whitespace remains in gene symbols after cleaning")
}

date_like <- grepl(
  "^[0-9]{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$",
  gene_ids,
  ignore.case = TRUE
)
n_unresolved_date_like <- sum(date_like, na.rm = TRUE)
if (n_unresolved_date_like > 0L) {
  example_dates <- paste(head(unique(gene_ids[date_like]), 10), collapse = ", ")
  warning(
    sprintf(
      "%d date-like gene symbol(s) are outside the fixed MARCH/SEPT/DEC map: %s",
      n_unresolved_date_like, example_dates
    )
  )
}

# Final requirement: duplicate gene symbols are always merged by the mean value.
if (anyDuplicated(gene_ids)) {
  n_before <- length(gene_ids)
  groups <- split(seq_len(nrow(rt)), gene_ids)
  rows <- lapply(groups, function(index) {
    if (length(index) == 1L) return(rt[index, , drop = FALSE])
    apply(rt[index, , drop = FALSE], 2, mean, na.rm = TRUE)
  })
  rt <- do.call(rbind, rows)
  gene_ids <- names(groups)
  n_duplicates_collapsed <- n_before - length(gene_ids)
  collapse_method <- "mean"
  log_step(
    "collapsed %d duplicate gene ID row(s) using mean",
    n_duplicates_collapsed
  )
}
rownames(rt) <- gene_ids
if (anyDuplicated(rownames(rt))) stopf("duplicate gene IDs remain after collapsing")
if (any(grepl("[[:space:]]", rownames(rt)))) stopf("gene symbols contain whitespace")

# ---------------------------------------------------------------------------
# Metadata alignment
# ---------------------------------------------------------------------------

if (is.null(metadata)) {
  metadata <- data.frame(sample_id = colnames(rt), stringsAsFactors = FALSE)
} else {
  metadata$sample_id <- trimws(as.character(metadata$sample_id))
  if (anyDuplicated(metadata$sample_id)) stopf("metadata has duplicate sample_id values")
  missing_meta <- setdiff(colnames(rt), metadata$sample_id)
  extra_meta <- setdiff(metadata$sample_id, colnames(rt))
  if (length(missing_meta)) {
    warning(sprintf("%d matrix sample(s) lack metadata; adding NA rows", length(missing_meta)))
    add <- data.frame(sample_id = missing_meta, stringsAsFactors = FALSE)
    for (column in setdiff(colnames(metadata), "sample_id")) add[[column]] <- NA
    metadata <- rbind(metadata, add[, colnames(metadata), drop = FALSE])
  }
  if (length(extra_meta)) {
    log_step("metadata has %d sample(s) not in the matrix; they are ignored", length(extra_meta))
  }
  metadata <- metadata[match(colnames(rt), metadata$sample_id), , drop = FALSE]
  rownames(metadata) <- NULL
}
if (!identical(colnames(rt), as.character(metadata$sample_id))) {
  stopf("internal error: metadata order does not match matrix columns")
}

sample_qc <- function(mat, suffix) {
  data.frame(
    sample_id = colnames(mat),
    n_na = colSums(is.na(mat)),
    pct_na = round(100 * colMeans(is.na(mat)), 2),
    min = apply(mat, 2, min, na.rm = TRUE),
    median = apply(mat, 2, stats::median, na.rm = TRUE),
    max = apply(mat, 2, max, na.rm = TRUE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
qc_before <- sample_qc(rt, "before")

# KNN needs at least one observation per row/column after all-NA filtering.
all_na_rows <- rowSums(!is.na(rt)) == 0L
if (any(all_na_rows)) {
  log_step("dropped %d all-NA gene row(s); KNN cannot impute them", sum(all_na_rows))
  rt <- rt[!all_na_rows, , drop = FALSE]
}
all_na_samples <- colSums(!is.na(rt)) == 0L
if (any(all_na_samples)) {
  if (!flag("drop-all-na-samples")) {
    stopf(
      "these samples are entirely NA: %s. Fix the matrix or pass --drop-all-na-samples.",
      paste(colnames(rt)[all_na_samples], collapse = ", ")
    )
  }
  n_all_na_samples_removed <- sum(all_na_samples)
  log_step("dropped %d all-NA sample(s)", n_all_na_samples_removed)
  rt <- rt[, !all_na_samples, drop = FALSE]
  metadata <- metadata[!all_na_samples, , drop = FALSE]
  rownames(metadata) <- NULL
}
if (nrow(rt) < 2L || ncol(rt) < 2L) stopf("fewer than 2 genes or 2 samples remain")

n_na_input <- sum(is.na(rt))
if (n_na_input > 0L) {
  if (!requireNamespace("impute", quietly = TRUE)) {
    stop(
      paste0(
        "missing values detected but the Bioconductor package 'impute' is not installed.\n",
        "Install it with: BiocManager::install('impute')"
      ),
      call. = FALSE
    )
  }
  knn_k <- as_int(arg("knn-k", 10), 10)
  knn_rowmax <- as_num(arg("knn-rowmax", 0.5), 0.5)
  knn_colmax <- as_num(arg("knn-colmax", 0.8), 0.8)
  knn_maxp <- as_int(arg("knn-maxp", 1500), 1500)
  col_missing <- colSums(is.na(rt))
  bad_columns <- colnames(rt)[col_missing > knn_colmax * nrow(rt)]
  if (length(bad_columns)) {
    stopf(
      paste0(
        "these samples exceed --knn-colmax (%.0f%% missing): %s. ",
        "Remove them explicitly or raise --knn-colmax."
      ),
      100 * knn_colmax, paste(bad_columns, collapse = ", ")
    )
  }
  row_missing_fraction <- rowMeans(is.na(rt))
  mean_imputed_rows <- sum(row_missing_fraction > knn_rowmax)
  set.seed(1)
  imputed <- impute::impute.knn(
    rt,
    k = knn_k,
    rowmax = knn_rowmax,
    colmax = knn_colmax,
    maxp = knn_maxp,
    rng.seed = 1
  )
  rt <- imputed$data
  if (anyNA(rt)) stopf("KNN imputation left %d NA cell(s)", sum(is.na(rt)))
  log_step(
    "KNN-imputed %d cell(s): k=%d, rowmax=%.2f, colmax=%.2f, maxp=%d; %d row(s) used the package's mean fallback",
    n_na_input, knn_k, knn_rowmax, knn_colmax, knn_maxp, mean_imputed_rows
  )
} else {
  knn_k <- NA_integer_
  knn_rowmax <- NA_real_
  knn_colmax <- NA_real_
  knn_maxp <- NA_integer_
  mean_imputed_rows <- 0L
  log_step("no missing values; KNN imputation skipped")
}

# Negative values: shift by the global minimum, exactly as specified.
negative_shift <- 0
min_before_shift <- min(rt)
if (min_before_shift < 0) {
  negative_shift <- -min_before_shift
  rt <- rt + negative_shift
  log_step("negative minimum %.6g; shifted all values by +%.6g", min_before_shift, negative_shift)
} else {
  log_step("no negative values; no shift applied")
}
rt[rt < 0 & rt > -1e-8] <- 0

# Range rule: if any value is outside [0, 30], apply log2(x + 1).
min_before_log <- min(rt)
max_before_log <- max(rt)
log2_applied <- min_before_log < 0 || max_before_log > 30
if (log2_applied) {
  rt <- log2(rt + 1)
  log_step(
    "range [%.6g, %.6g] is outside [0, 30]; applied log2(x + 1)",
    min_before_log, max_before_log
  )
} else {
  log_step("range [%.6g, %.6g] is inside [0, 30]; log2 skipped", min_before_log, max_before_log)
}
if (any(!is.finite(rt))) stopf("non-finite values produced by transformation")
if (anyNA(rt)) stopf("NA values remain after processing")

# ---------------------------------------------------------------------------
# Output metadata / QC
# ---------------------------------------------------------------------------

base_type <- switch(
  value_type,
  tpm = "TPM",
  fpkm = "FPKM",
  rpkm = "RPKM",
  counts = "counts",
  array_linear = "intensity",
  array_log2 = "log2(intensity)",
  normalized = "normalized expression",
  "expression"
)
output_type <- if (log2_applied) sprintf("log2(%s + 1)", base_type) else base_type
if (identical(value_type, "counts")) {
  warning(
    paste0(
      "input is raw counts: log2(counts + 1) is NOT TPM and is not ",
      "library-size normalized.  Use this matrix only for an explicitly ",
      "count-based workflow."
    )
  )
}
if (identical(value_type, "array_log2") && log2_applied) {
  warning(
    paste0(
      "--value-type array_log2 but the range triggered another log2. ",
      "Check whether the input was already log2; a double log is usually wrong."
    )
  )
}

transformations <- c(
  if (n_date_genes_fixed > 0L) "fix_excel_date_gene_symbols",
  if (n_na_input > 0L) "knn_impute",
  if (negative_shift > 0) "negative_shift",
  if (log2_applied) "log2",
  if (n_duplicates_collapsed > 0L) paste0("collapse_duplicates_by_", collapse_method)
)
if (!length(transformations)) transformations <- "none"

qc_after <- sample_qc(rt, "after")
sample_qc_out <- merge(
  qc_before,
  qc_after,
  by = "sample_id",
  suffixes = c("_before", "_after"),
  all = TRUE,
  sort = FALSE
)
sample_qc_out <- sample_qc_out[match(colnames(rt), sample_qc_out$sample_id), , drop = FALSE]

mapping_rate_pct <- if (is.na(mapping_rate)) NA_real_ else round(100 * mapping_rate, 2)
qc_summary <- data.frame(
  metric = c(
    "n_genes_output", "n_samples_output", "n_na_input", "pct_na_input",
    "n_nonnumeric_cells", "n_all_na_genes_removed", "n_all_na_samples_removed",
    "n_date_genes_fixed", "n_unresolved_date_like",
    "n_duplicate_genes_collapsed", "n_multi_symbol_rows", "n_unmapped_probes",
    "probe_mapping_rate_pct", "n_negative_shift", "min_before_shift",
    "min_before_log", "max_before_log", "log2_applied",
    "min_output", "max_output", "median_output", "n_na_output",
    "knn_k", "knn_rowmax", "knn_colmax", "knn_maxp",
    "knn_mean_fallback_rows", "value_type_input", "value_type_output"
  ),
  value = c(
    nrow(rt), ncol(rt), n_na_input, round(100 * n_na_input / (nrow(rt) * ncol(rt)), 4),
    nonnumeric_cells, sum(all_na_rows), n_all_na_samples_removed,
    n_date_genes_fixed, n_unresolved_date_like,
    n_duplicates_collapsed, n_multi_symbol, n_unmapped,
    mapping_rate_pct, negative_shift, min_before_shift,
    min_before_log, max_before_log, log2_applied,
    min(rt), max(rt), stats::median(rt), sum(is.na(rt)),
    knn_k, knn_rowmax, knn_colmax, knn_maxp,
    mean_imputed_rows, value_type, output_type
  ),
  stringsAsFactors = FALSE
)

# Final delivery object is a data.frame, not a matrix:
#   rows    = unique gene symbols (no whitespace, no "///")
#   columns = sample IDs only
rt <- as.data.frame(rt, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(vapply(rt, is.numeric, logical(1)))) {
  stopf("internal error: not all expression columns are numeric")
}
if (anyNA(rt)) stopf("internal error: NA values remain in the final data.frame")
if (anyDuplicated(rownames(rt))) stopf("internal error: duplicate gene symbols remain")
if (any(grepl("[[:space:]]", rownames(rt)))) {
  stopf("internal error: gene symbols still contain whitespace")
}

dataset_info <- list(
  accession = if (nzchar(accession)) accession else NA_character_,
  source = arg("source", "GEO"),
  platform = arg("platform", NA_character_),
  technology = arg(
    "technology",
    if (value_type %in% c("tpm", "fpkm", "rpkm", "counts")) "bulk RNA-seq" else if (
      value_type %in% c("array_linear", "array_log2")
    ) "microarray" else "unknown"
  ),
  species = arg("species", NA_character_),
  input_value_type = value_type,
  output_value_type = output_type,
  output_class = "data.frame",
  gene_id_rule = paste(
    "row names are gene symbols; Excel MARCH/SEPT/DEC dates restored;",
    "first symbol before '///'; all whitespace removed;",
    if (n_duplicates_collapsed > 0L) {
      paste("duplicates merged by", collapse_method)
    } else {
      "no duplicate gene symbols"
    }
  ),
  date_gene_fix = list(
    families = c("MARCH1-11", "SEPT1-15", "DEC1-2"),
    n_fixed = n_date_genes_fixed,
    n_unresolved_date_like = n_unresolved_date_like
  ),
  collapse = list(method = collapse_method, n_collapsed = n_duplicates_collapsed),
  n_genes = nrow(rt),
  n_samples = ncol(rt),
  transformations = transformations,
  negative_shift = negative_shift,
  log2_applied = log2_applied,
  knn = list(
    method = if (n_na_input > 0L) "impute::impute.knn" else "skipped",
    k = knn_k,
    rowmax = knn_rowmax,
    colmax = knn_colmax,
    maxp = knn_maxp,
    n_imputed_cells = n_na_input,
    mean_fallback_rows = mean_imputed_rows
  ),
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version = R.version.string,
  input_files = normalizePath(input_path, mustWork = FALSE),
  input_md5 = unname(tools::md5sum(input_path))
)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

out_rdata <- file.path(outdir, paste0(prefix, "_expression.RData"))
if (file.exists(out_rdata) && !flag("force")) {
  stopf("output exists: %s (use --force to overwrite)", out_rdata)
}
compress_arg <- tolower(arg("compress", "gzip"))
if (!compress_arg %in% c("gzip", "xz", "none")) stopf("--compress must be gzip, xz, or none")
save(
  rt, metadata, dataset_info,
  file = out_rdata,
  compress = if (identical(compress_arg, "none")) FALSE else compress_arg
)
log_step("saved RData: %s", out_rdata)

write.table(
  qc_summary,
  file = file.path(outdir, paste0(prefix, "_qc_summary.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  sample_qc_out,
  file = file.path(outdir, paste0(prefix, "_sample_qc.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  metadata,
  file = file.path(outdir, paste0(prefix, "_metadata.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
writeLines(capture.output(sessionInfo()), file.path(outdir, paste0(prefix, "_sessionInfo.txt")))

if (flag("export-tsv")) {
  export_path <- file.path(outdir, paste0(prefix, "_expression.tsv.gz"))
  con <- gzfile(export_path, "wt")
  write.table(
    data.frame(gene_id = rownames(rt), rt, check.names = FALSE),
    con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
  close(con)
  log_step("exported TSV: %s", export_path)
}

if (flag("rds")) {
  saveRDS(
    list(rt = rt, metadata = metadata, dataset_info = dataset_info),
    file = file.path(outdir, paste0(prefix, "_expression.rds"))
  )
  log_step("saved RDS: %s", file.path(outdir, paste0(prefix, "_expression.rds")))
}

writeLines(log_lines, file.path(outdir, paste0(prefix, "_processing_log.txt")))

message("")
message("============================================================")
message(sprintf("Done: %d genes x %d samples", nrow(rt), ncol(rt)))
message(sprintf("Value type: %s", output_type))
message(sprintf("Transformations: %s", paste(transformations, collapse = ", ")))
message(sprintf("Range: [%.4f, %.4f]", min(rt), max(rt)))
message(sprintf("RData: %s", out_rdata))
message("Next: ask the user about normal-tissue removal and clinical data.")
message("============================================================")
