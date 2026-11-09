#!/usr/bin/env Rscript
# Validate and standardize a clinical / prognosis table, optionally matching it
# to a prepared expression RData.
#
# Standard survival columns produced:
#   OS_time_months,  OS_event
#   PFS_time_months, PFS_event
#   DFS_time_months, DFS_event
#   DSS_time_months, DSS_event
#
# Event coding: 1 = event (death / progression / recurrence), 0 = censored.
# Time unit: suffix _months or _days is mandatory in the output.  If the input
# has bare OS_time, pass --time-unit months (or days); the script will rename it.
#
# Example:
#   Rscript validate_clinical.R \
#     --clinical clinical.csv \
#     --rdata GSE50760/GSE50760_expression.RData \
#     --outdir GSE50760/clinical \
#     --time-unit months

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

read_table_auto <- function(path) {
  lower <- tolower(path)
  if (grepl("\\.xlsx?$", lower)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) stopf("reading xlsx needs openxlsx")
    return(openxlsx::read.xlsx(path, sheet = 1, check.names = FALSE))
  }
  con <- if (grepl("\\.gz$", lower)) gzfile(path, "rt") else file(path, "rt")
  lines <- tryCatch(readLines(con, n = 20L, warn = FALSE), finally = close(con))
  lines <- lines[nzchar(trimws(lines))]
  counts <- c(
    tab = lengths(regmatches(lines[[1]], gregexpr("\t", lines[[1]], fixed = TRUE))),
    comma = lengths(regmatches(lines[[1]], gregexpr(",", lines[[1]], fixed = TRUE))),
    semi = lengths(regmatches(lines[[1]], gregexpr(";", lines[[1]], fixed = TRUE)))
  )
  sep <- c("\t", ",", ";")[[which.max(counts)]]
  out <- read.delim(
    path, sep = sep, header = TRUE, check.names = FALSE, quote = "\"",
    comment.char = "", na.strings = c("NA", "N/A", "null", "NULL", "NaN", ""),
    stringsAsFactors = FALSE
  )
  if (ncol(out) > 0L) colnames(out)[[1]] <- sub("^\ufeff", "", colnames(out)[[1]])
  out
}

clinical_path <- normalizePath(require_arg("clinical"), mustWork = TRUE)
outdir <- require_arg("outdir")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

frame <- read_table_auto(clinical_path)
if (nrow(frame) == 0L) stopf("clinical table is empty")
if (!"patient_id" %in% colnames(frame)) {
  stopf("clinical table must contain a patient_id column")
}
frame$patient_id <- trimws(as.character(frame$patient_id))
if (any(!nzchar(frame$patient_id))) stopf("empty patient_id values")
if (!"sample_id" %in% colnames(frame)) frame$sample_id <- NA_character_
frame$sample_id <- trimws(as.character(frame$sample_id))

time_unit <- tolower(arg("time-unit", ""))
if (nzchar(time_unit) && !time_unit %in% c("months", "days")) {
  stopf("--time-unit must be months or days")
}

checks <- data.frame(check = character(), status = character(), detail = character())
add_check <- function(check, status, detail) {
  checks <<- rbind(
    checks,
    data.frame(check = check, status = status, detail = detail, stringsAsFactors = FALSE)
  )
}

endpoints <- c("OS", "PFS", "DFS", "DSS")
standard_cols <- c("patient_id", "sample_id")
for (endpoint in endpoints) {
  time_candidates <- c(
    paste0(endpoint, "_time_", c("months", "days")),
    paste0(endpoint, ".time.", c("months", "days")),
    paste0(endpoint, "_time"),
    paste0(endpoint, ".time")
  )
  event_candidates <- c(paste0(endpoint, "_event"), paste0(endpoint, ".event"))
  time_col <- time_candidates[time_candidates %in% colnames(frame)][1]
  event_col <- event_candidates[event_candidates %in% colnames(frame)][1]
  if (is.na(time_col) && is.na(event_col)) next
  if (is.na(time_col) || is.na(event_col)) {
    add_check(endpoint, "ERROR", "time and event must both be present")
    next
  }

  time_suffix <- if (grepl("_months$|\\.months$", time_col)) {
    "months"
  } else if (grepl("_days$|\\.days$", time_col)) {
    "days"
  } else {
    time_unit
  }
  if (!nzchar(time_suffix)) {
    add_check(endpoint, "ERROR", sprintf("%s has no time unit; pass --time-unit", time_col))
    next
  }
  if (nzchar(time_unit) && !identical(time_suffix, time_unit)) {
    add_check(
      endpoint, "ERROR",
      sprintf("%s is in %s but --time-unit %s was supplied", time_col, time_suffix, time_unit)
    )
    next
  }

  event_values <- frame[[event_col]]
  numeric_event <- suppressWarnings(as.numeric(event_values))
  valid_event <- !is.na(numeric_event) & numeric_event %in% c(0, 1)
  if (!all(valid_event | is.na(event_values))) {
    bad <- unique(event_values[!valid_event & !is.na(event_values)])
    add_check(
      endpoint, "ERROR",
      sprintf("%s must be coded 0/1; invalid: %s", event_col, paste(bad, collapse = ", "))
    )
  } else if (any(is.na(event_values))) {
    add_check(endpoint, "WARN", sprintf("%s has missing event values", event_col))
  } else {
    add_check(endpoint, "OK", sprintf("%s: 0/1 coding", event_col))
  }

  time_values <- suppressWarnings(as.numeric(frame[[time_col]]))
  if (all(is.na(time_values))) {
    add_check(endpoint, "ERROR", sprintf("%s has no numeric values", time_col))
  } else if (any(time_values < 0, na.rm = TRUE)) {
    add_check(endpoint, "ERROR", sprintf("%s has negative times", time_col))
  } else if (any(time_values == 0, na.rm = TRUE)) {
    add_check(endpoint, "WARN", sprintf("%s has zero follow-up", time_col))
  } else {
    add_check(endpoint, "OK", sprintf("%s: numeric %s", time_col, time_suffix))
  }
  if (any(numeric_event == 1 & is.na(time_values), na.rm = TRUE)) {
    add_check(endpoint, "ERROR", "event=1 but time is missing")
  }
  if (any(numeric_event == 0 & is.na(time_values), na.rm = TRUE)) {
    add_check(endpoint, "WARN", "event=0 but time is missing")
  }

  new_time <- paste0(endpoint, "_time_", time_suffix)
  new_event <- paste0(endpoint, "_event")
  frame[[new_time]] <- time_values
  frame[[new_event]] <- numeric_event
  if (!identical(time_col, new_time)) frame[[time_col]] <- NULL
  standard_cols <- c(standard_cols, new_time, new_event)
}

if (length(standard_cols) == 2L) {
  add_check("survival", "WARN", "no OS/PFS/DFS/DSS time+event pair found")
}

# Patient-level uniqueness.
if (anyDuplicated(frame$patient_id)) {
  if (!flag("allow-multiple-samples")) {
    add_check("patient_id", "ERROR", "duplicate patient_id; use --allow-multiple-samples")
  } else if (all(is.na(frame$sample_id))) {
    add_check("patient_id", "ERROR", "multiple rows per patient require sample_id")
  } else {
    add_check("patient_id", "OK", "multiple samples per patient allowed")
  }
} else {
  add_check("patient_id", "OK", "unique")
}

# Age and sex are optional; check them only if present.
if ("age_years" %in% colnames(frame) || "age" %in% colnames(frame)) {
  age_col <- if ("age_years" %in% colnames(frame)) "age_years" else "age"
  age <- suppressWarnings(as.numeric(frame[[age_col]]))
  if (any(age < 0 | age > 120, na.rm = TRUE)) {
    add_check("age", "WARN", sprintf("%s outside 0-120", age_col))
  } else {
    add_check("age", "OK", age_col)
  }
}
if ("sex" %in% colnames(frame)) {
  add_check("sex", "OK", paste(unique(frame$sex[!is.na(frame$sex)]), collapse = ", "))
}

# Match sample IDs to the expression matrix, when available.
rdata_path <- arg("rdata", NULL)
matched <- NULL
if (!is.null(rdata_path) && nzchar(as.character(rdata_path))) {
  env <- new.env(parent = baseenv())
  load(normalizePath(as.character(rdata_path), mustWork = TRUE), envir = env)
  if (is.null(env$rt) || !(is.data.frame(env$rt) || is.matrix(env$rt))) {
    stopf("RData has no data.frame/matrix 'rt'")
  }
  matrix_samples <- colnames(env$rt)
  has_sample <- !is.na(frame$sample_id) & nzchar(frame$sample_id)
  in_matrix <- has_sample & frame$sample_id %in% matrix_samples
  frame$in_expression <- in_matrix
  matched <- data.frame(
    patient_id = frame$patient_id,
    sample_id = frame$sample_id,
    in_expression = in_matrix,
    stringsAsFactors = FALSE
  )
  n_unmatched <- sum(has_sample & !in_matrix)
  add_check(
    "expression_match", if (n_unmatched == 0L) "OK" else "WARN",
    sprintf("%d clinical sample(s) not found in expression matrix", n_unmatched)
  )
  if (flag("require-match") && n_unmatched > 0L) {
    stopf("%d clinical sample(s) do not match the expression matrix", n_unmatched)
  }
  missing_clinical <- setdiff(matrix_samples, frame$sample_id[has_sample])
  if (length(missing_clinical)) {
    add_check(
      "expression_match", "WARN",
      sprintf("%d expression sample(s) have no clinical row", length(missing_clinical))
    )
  }
}

errors <- checks$status == "ERROR"
if (any(errors)) {
  print(checks)
  stopf("clinical validation failed; fix the ERROR rows above")
}

# Put standardized columns first, preserve all other columns.
extras <- setdiff(colnames(frame), standard_cols)
clinical <- frame[, c(standard_cols, extras), drop = FALSE]
rownames(clinical) <- NULL

clinical_info <- list(
  n_rows = nrow(clinical),
  n_patients = length(unique(clinical$patient_id)),
  time_unit = if (nzchar(time_unit)) time_unit else "per-column suffix",
  endpoints = endpoints[
    paste0(endpoints, "_event") %in% colnames(clinical)
  ],
  source_file = normalizePath(clinical_path, mustWork = TRUE),
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  validation = checks
)

write.csv(
  clinical,
  file = file.path(outdir, "clinical_clean.csv"),
  row.names = FALSE, na = ""
)
write.table(
  checks,
  file = file.path(outdir, "clinical_validation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
if (!is.null(matched)) {
  write.table(
    matched,
    file = file.path(outdir, "sample_map.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}
if (!is.null(matched)) {
  sample_map <- matched
  save(clinical, clinical_info, sample_map, file = file.path(outdir, "clinical.RData"), compress = "gzip")
} else {
  save(clinical, clinical_info, file = file.path(outdir, "clinical.RData"), compress = "gzip")
}

# A compact analysis-ready survival table.
survival_cols <- c("patient_id", "sample_id")
for (endpoint in endpoints) {
  time_cols <- c(paste0(endpoint, "_time_months"), paste0(endpoint, "_time_days"))
  time_col <- time_cols[time_cols %in% colnames(clinical)][1]
  event_col <- paste0(endpoint, "_event")
  if (!is.na(time_col) && event_col %in% colnames(clinical)) {
    survival_cols <- c(survival_cols, time_col, event_col)
  }
}
if (length(survival_cols) > 2L) {
  write.csv(
    clinical[, survival_cols, drop = FALSE],
    file = file.path(outdir, "survival_analysis.csv"),
    row.names = FALSE, na = ""
  )
}

print(checks)
message("")
message(sprintf("[clinical] %d rows, %d patients", nrow(clinical), clinical_info$n_patients))
message(sprintf("[clinical] saved: %s", file.path(outdir, "clinical.RData")))
message("[clinical] times are explicitly suffixed; events are 0/1.")
