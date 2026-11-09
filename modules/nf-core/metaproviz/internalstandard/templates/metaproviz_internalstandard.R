#!/usr/bin/env Rscript
# =============================================================================
# metaproviz_internalstandard.R — technical-QC check on the internal-standard spike-ins
#
# Approach: re-use MetaProViz::pool_estimation on the SE subsetted to just the
# IS features, with every sample flagged as "Pool" (because for IS we want the
# CV across all samples — IS should be constant, so a high CV reveals
# technical issues, not biology).
#
# Picking the IS features (top precedence first):
#   1. --internal_standards  valine-d8,hippuric acid-d5,...     (explicit list)
#   2. --pattern             "regex"                            (user regex)
#   3. Default auto-detect on the common IS naming conventions:
#         -dN suffix         (deuterated, e.g. valine-d8)
#         -13CN suffix       (carbon-13 labelled)
#         ISTD anywhere
#         IS_ prefix / _IS suffix / bare 'IS'
#
# The CV is computed across all samples of the input. Pool/QC samples must be removed before this module (for example
# with the metaproviz/poolestimation module), they are not detected here and would be used like all other samples.
#
# The IS are removed from the SE, which is written as <prefix>.is_filtered.rds. The removed IS are the rows of
# <prefix>.cv.tsv (one row per IS, with the CV and the HighVar flag), their values are in the input data.
#
# No hard stop when the IS have no usable values (all not detected, one value only, identical values): the report and
# the log say which of these cases applies. The module stops only if no IS is found or if every feature is an IS.
# =============================================================================

set.seed(42)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(magrittr)
  library(ggplot2); library(SummarizedExperiment); library(S4Vectors)
})

suppressPackageStartupMessages(library(MetaProViz))

# ── PARSE PARAMETERS FROM NEXTFLOW ──────────────────────────────────────────
nullify <- function(x) if (is.null(x) || x == "" || x == "null" || x == "[]") NULL else x

# Options in list_opts hold a comma-separated list. Their value is everything after the option name, split at the
# commas, so names with spaces need no quoting. A name with a comma goes in double quotes.
parse_args <- function(x, list_opts = character(0)) {
  if (is.null(x) || !nzchar(trimws(x))) return(list())
  args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
  parsed <- list()
  for (a in args_list) {
    nm <- scan(text = a, what = 'character', quiet = TRUE, nmax = 1)
    if (length(nm) == 0) next
    if (nm %in% list_opts) {
      rest <- trimws(sub("^[^[:space:]]+", "", a))
      parsed[[nm]] <- scan(text = rest, what = 'character', sep = ",", quote = '"', strip.white = TRUE, quiet = TRUE)
    } else {
      z <- scan(text = a, what = 'character', quiet = TRUE)
      length(z) <- 2
      if (!is.na(z[2])) parsed[[z[1]]] <- z[2]
    }
  }
  parsed
}

opt <- list(
  internal_standards  = character(0),
  pattern             = "",
  cutoff_cv           = "30",
  conditions_col      = "Conditions",
  id_columns          = character(0)
)
# R raw string: quotes and apostrophes inside the arguments (feature names like 5'-AMP) cannot break the R code.
args_opt <- parse_args(r"---($task.ext.args)---", list_opts = c("internal_standards", "id_columns"))
for (ao in names(args_opt)) {
  if (!ao %in% names(opt)) stop(paste("Invalid option:", ao))
  opt[[ao]] <- args_opt[[ao]]
}

se_path             <- nullify('$se_rds')
data_matrix_path    <- nullify('$data_matrix')
feature_matrix_path <- nullify('$feature_matrix')
sample_matrix_path  <- nullify('$sample_matrix')
internal_standards_list <- opt\$internal_standards
is_list_str          <- paste(internal_standards_list, collapse = ",")
user_pattern          <- opt\$pattern
cutoff_cv             <- as.numeric(opt\$cutoff_cv)

# colData column holding experimental condition, used to break the IS CV out
# per condition (condition_cv output). Column not found => condition_cv is
# written header-only.
conditions_col        <- opt\$conditions_col
id_columns_list <- opt\$id_columns
id_columns_str        <- paste(id_columns_list, collapse = ",")

prefix               <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
# output_format is a val input of the module.
# se (or an empty value): the SE without the IS as one RDS. tsv/csv: three files (data, feature and sample matrix), the same shape as the
# input matrices.
output_format_in    <- trimws('$output_format')
output_format       <- if (output_format_in %in% c("", "null", "[]")) "se" else tolower(output_format_in)
out_files           <- function(base) if (output_format == "se") paste0(base, ".rds") else
  paste0(base, c(".data_matrix.", ".feature_matrix.", ".sample_matrix."), output_format)
se_out              <- out_files(paste0(prefix, ".is_filtered"))
cv_out              <- paste0(prefix, ".cv.tsv")
condition_cv_out    <- paste0(prefix, ".condition_cv.tsv")
plots_out           <- paste0(prefix, ".plots.rds")
report_out          <- paste0(prefix, ".report.html")
log_out             <- paste0(prefix, ".log")

# Default regex — case-insensitive PCRE. Covers common IS labels.
DEFAULT_PATTERN <- "-d\\\\d+|-13C\\\\d*|ISTD|^IS_|_IS\$|^IS\$"
# Feature metadata columns that hold database IDs or annotations such as pathways (found by name, or set with
# --id_columns).
ID_PATTERN <- "hmdb|kegg|chebi|pubchem|inchi|smiles|lipidmaps|swisslipids|metacyc|biocyc|chembl|drugbank|pathway|(^|[._ -])cas([._ -]|\$)"

# ── LOGGING ─────────────────────────────────────────────────────────────────
ll <- character(0); wc <- 0; ec <- 0
log_msg <- function(lev = "INFO", ...) {
  txt <- paste0("[", lev, "] ", paste(..., sep = ""))
  ll <<- c(ll, txt); message(txt)
  if (lev == "WARN")  wc <<- wc + 1
  if (lev == "ERROR") ec <<- ec + 1
}
log_section <- function(t) {
  sep <- paste(rep("─", 60), collapse = "")
  ll <<- c(ll, "", sep, paste0("  ", t), sep)
  message(sep, "\\n  ", t, "\\n", sep)
}
flush_log <- function(p = log_out) {
  writeLines(c(ll, "",
    sprintf("SUMMARY: %d warning(s), %d error(s)", wc, ec)), p)
}
abort <- function(...) { log_msg("ERROR", ...); flush_log(); stop(paste0(...), call. = FALSE) }
if (!output_format %in% c("se", "tsv", "csv"))
  abort("output_format must be se, tsv or csv (empty means se), not '", output_format_in, "'.")
esc <- function(x) gsub(">", "&gt;", gsub("<", "&lt;",
                    gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE), fixed = TRUE)
# Names for the report: at most 25, and a pointer to the log file only if not all of them are shown.
names_html <- function(x, max_n = 25) paste0(
  paste0("<code>", esc(head(x, max_n)), "</code>", collapse = " "),
  if (length(x) > max_n) paste0(" and ", length(x) - max_n, " more (see <code>", esc(log_out), "</code>)") else "")

# IS names with a detail in brackets (the value and the sample), at most 25 shown.
info_html <- function(ids, info, max_n = 25) paste0(
  paste0("<code>", esc(head(ids, max_n)), "</code> (", esc(head(info, max_n)), ")", collapse = ", "),
  if (length(ids) > max_n) paste0(" and ", length(ids) - max_n, " more (see <code>", esc(log_out), "</code>)") else "")

# ── LOAD SE ─────────────────────────────────────────────────────────────────
log_section("Loading input data")

# Two acceptable input shapes, mutually exclusive: --se (a SummarizedExperiment
# .rds), or all three of --data_matrix/--feature_matrix/--sample_matrix (plain
# TSV/CSV, language-independent). Test data is provided in both shapes; this
# module accepts either so users aren't forced into one or the other.
have_se  <- !is.null(se_path) && nzchar(se_path)
have_csv <- !is.null(data_matrix_path) && nzchar(data_matrix_path) &&
            !is.null(feature_matrix_path) && nzchar(feature_matrix_path) &&
            !is.null(sample_matrix_path) && nzchar(sample_matrix_path)
have_partial_csv <- !have_csv && (
  (!is.null(data_matrix_path) && nzchar(data_matrix_path)) ||
  (!is.null(feature_matrix_path) && nzchar(feature_matrix_path)) ||
  (!is.null(sample_matrix_path) && nzchar(sample_matrix_path))
)

if (have_se && have_csv)
  abort("Provide either --se, or all three of --data_matrix/--feature_matrix/",
        "--sample_matrix, not both.")
if (have_se && have_partial_csv)
  abort("--se was provided together with some of --data_matrix/--feature_matrix/",
        "--sample_matrix; provide either --se alone, or all three of those, not a mix.")
if (have_partial_csv)
  abort("--data_matrix, --feature_matrix, and --sample_matrix must all be ",
        "provided together.")
if (!have_se && !have_csv)
  abort("Provide either --se <path>, or all three of --data_matrix/",
        "--feature_matrix/--sample_matrix.")

if (have_se) {
  if (!file.exists(se_path)) abort("SE file not found: ", se_path)
  se <- readRDS(se_path)
  if (!is(se, "SummarizedExperiment")) abort("Object is not a SummarizedExperiment.")
  log_msg("INFO", "Loaded SE from ", basename(se_path), ": ", nrow(se),
          " features × ", ncol(se), " samples.")
} else {
  for (p in c(data_matrix_path, feature_matrix_path, sample_matrix_path))
    if (!file.exists(p)) abort("Input file not found: ", p)

  read_flat <- function(path) {
    sep <- if (grepl("\\\\.csv\$", path, ignore.case = TRUE)) "," else "\\t"
    tryCatch(read.delim(path, sep = sep, check.names = FALSE,
                        stringsAsFactors = FALSE),
             error = function(e) abort("Failed to read ", path, ": ",
                                       conditionMessage(e)))
  }

  dm <- read_flat(data_matrix_path)
  fm <- read_flat(feature_matrix_path)
  sm <- read_flat(sample_matrix_path)

  if (ncol(dm) < 2)
    abort("--data_matrix must have a feature-ID column plus at least one ",
          "sample column.")

  feature_ids <- as.character(dm[[1]])
  assay_matrix <- as.matrix(dm[, -1, drop = FALSE])
  rownames(assay_matrix) <- feature_ids
  mode(assay_matrix) <- "numeric"

  sample_ids <- as.character(sm[[1]])
  col_data <- sm[, -1, drop = FALSE]
  rownames(col_data) <- sample_ids

  feature_ids2 <- as.character(fm[[1]])
  row_data <- fm[, -1, drop = FALSE]
  rownames(row_data) <- feature_ids2

  missing_samples <- setdiff(colnames(assay_matrix), sample_ids)
  if (length(missing_samples) > 0)
    abort("--data_matrix has sample column(s) not present in ",
          "--sample_matrix: ", paste(missing_samples, collapse = ", "))
  missing_features <- setdiff(rownames(assay_matrix), feature_ids2)
  if (length(missing_features) > 0)
    abort("--data_matrix has feature(s) not present in --feature_matrix: ",
          paste(head(missing_features, 5), collapse = ", "),
          if (length(missing_features) > 5) ", ..." else "")

  col_data <- col_data[colnames(assay_matrix), , drop = FALSE]
  row_data <- row_data[rownames(assay_matrix), , drop = FALSE]

  se <- SummarizedExperiment(assays = list(counts = assay_matrix),
                             colData = col_data, rowData = row_data)
  log_msg("INFO", "Built SE from data_matrix/feature_matrix/sample_matrix: ",
          nrow(se), " features × ", ncol(se), " samples.")
}

source_label <- if (have_se) {
  basename(se_path)
} else {
  paste(basename(data_matrix_path), basename(feature_matrix_path),
        basename(sample_matrix_path), sep = " + ")
}

n_samples_input <- ncol(se)

# ── IDENTIFY IS FEATURES ────────────────────────────────────────────────────
log_section("Identifying internal-standard features")

feat_names <- rownames(se)
is_feats <- character(0)
strategy <- "none"
is_attempted <- character(0)   # ways of finding the IS that were tried and found no feature

# 1) Explicit list
if (nzchar(is_list_str)) {
  asked <- trimws(internal_standards_list)
  asked <- asked[nzchar(asked)]
  is_feats <- intersect(asked, feat_names)
  missing  <- setdiff(asked, feat_names)
  if (length(missing) > 0)
    log_msg("WARN", "Requested IS feature(s) not found in SE: ",
            paste(missing, collapse = ", "))
  strategy <- "explicit_list"
  if (length(is_feats) == 0) is_attempted <- c(is_attempted, list = "explicit list (--internal_standards)")
}

# 2) User pattern
if (length(is_feats) == 0 && nzchar(user_pattern)) {
  hit <- grepl(user_pattern, feat_names, ignore.case = TRUE, perl = TRUE)
  is_feats <- feat_names[hit]
  strategy <- "user_pattern"
  if (length(is_feats) == 0) is_attempted <- c(is_attempted, pattern = sprintf("pattern '%s' (--pattern)", user_pattern))
}

# 3) Default pattern
id_status  <- "not run"
id_cols    <- character(0)
id_flagged <- character(0)
if (length(is_feats) == 0) {
  hit <- grepl(DEFAULT_PATTERN, feat_names, ignore.case = TRUE, perl = TRUE)
  is_feats <- feat_names[hit]
  strategy <- "default"

  # A name can look like a standard but belong to a real metabolite (Vitamin-D3 matches -d3). A real metabolite
  # usually has an entry in a database, a spike-in standard usually has none. So a feature found by the default pattern
  # that has a database ID in the feature metadata is flagged in the log and the report, but still used as IS. Only the
  # presence of an ID counts, a missing ID proves nothing. Names given with --internal_standards or --pattern are never
  # checked.
  if (length(is_feats) > 0) {
    feat_info <- as.data.frame(rowData(se), check.names = FALSE)
    id_cols <- if (nzchar(id_columns_str)) trimws(id_columns_list) else
      grep(ID_PATTERN, colnames(feat_info), ignore.case = TRUE, value = TRUE)
    id_cols <- intersect(id_cols, colnames(feat_info))
    has_id  <- function(x) {
      v <- trimws(iconv(as.character(x), "UTF-8", "UTF-8", sub = ""))   # invalid characters would stop trimws()
      !is.na(v) & !(toupper(v) %in% c("", "NA", "NAN", "N/A", "-"))
    }
    if (length(id_cols) == 0) {
      id_status <- "no_columns"
    } else if (!any(vapply(id_cols, function(cn) any(has_id(feat_info[[cn]])), logical(1)))) {
      id_status <- "empty"
    } else {
      id_status  <- "checked"
      id_matrix  <- do.call(cbind, lapply(id_cols, function(cn) has_id(feat_info[is_feats, cn])))
      dimnames(id_matrix) <- list(is_feats, id_cols)
      id_flagged <- is_feats[rowSums(id_matrix) > 0]
      # for each flagged name, the columns in which an ID was found
      id_found_cols <- lapply(setNames(id_flagged, id_flagged), function(f) id_cols[id_matrix[f, ]])
    }
  }
  if (id_status == "checked" && length(id_flagged) > 0) {
    log_msg("WARN", "Names with an entry in the database ID columns of the feature metadata: ",
            paste0(id_flagged, " (", vapply(id_found_cols, paste, character(1), collapse = ", "), ")", collapse = ", "),
            ". They are used as IS. Please check that they are real IS.")
  } else if (id_status == "checked") {
    log_msg("INFO", "Checked the names found by the default pattern against the feature metadata columns with database ",
            "IDs: ", paste(id_cols, collapse = ", "), ". None of the names has an entry in these columns. ",
            if (nzchar(id_columns_str)) "The columns were given with --id_columns." else
              "The columns were selected by default, --id_columns sets them directly.")
  } else if (id_status == "no_columns") {
    log_msg("INFO", "The feature metadata has no columns with database IDs, and the names found by the default pattern ",
            "could not be checked. Columns are recognised by names such as HMDB, KEGG or ChEBI, --id_columns sets them ",
            "directly.")
  } else if (id_status == "empty") {
    log_msg("INFO", "The feature metadata columns with database IDs (", paste(id_cols, collapse = ", "),
            ") are empty, and the names found by the default pattern could not be checked.")
  }
}
# Wording of the detection for the log and the report
is_label <- switch(strategy,
  explicit_list = "explicit list (--internal_standards)",
  user_pattern  = sprintf("pattern '%s' (--pattern)", user_pattern),
  "default naming conventions")
is_label_html <- switch(strategy,
  explicit_list = "explicit list <code>--internal_standards</code>",
  user_pattern  = sprintf("pattern <code>--pattern %s</code>", esc(user_pattern)),
  "default naming conventions: deuterated -dN, 13C, ISTD, IS_ prefix, _IS suffix")
is_attempted_txt  <- paste(is_attempted, collapse = " and the ")
is_attempted_html <- paste(vapply(names(is_attempted), function(k) if (k == "list") {
  "explicit list <code>--internal_standards</code>"
} else sprintf("pattern <code>--pattern %s</code>", esc(user_pattern)), character(1)), collapse = " and the ")
if (length(is_attempted) > 0 && strategy == "default" && length(is_feats) > 0)
  log_msg("WARN", "No IS features were found by the ", is_attempted_txt,
          ". The IS detection falls back to the default naming conventions.")

if (length(is_feats) == 0) {
  ways <- c(is_attempted, "default naming conventions")
  ways_txt <- if (length(ways) == 1) ways else
    paste0(paste(head(ways, -1), collapse = ", the "), " and the ", tail(ways, 1))
  abort("No IS features were found by the ", ways_txt, ". Please check the IS detection parameters of this module.")
}

log_msg("INFO", "IS detection strategy: ", is_label)
log_msg("INFO", "IS features (", length(is_feats), "): ",
        paste(head(is_feats, 12), collapse = ", "),
        if (length(is_feats) > 12) " ..." else "")

# ── SPLIT THE SE INTO IS AND THE REST ───────────────────────────────────────
log_section("Removing the IS from the SE")

if (length(is_feats) == nrow(se))
  abort("All ", nrow(se), " features were identified as IS (found by ", is_label, "). No feature would be left after ",
        "removing them. Please check the IS detection parameters of this module.")
se_is    <- se[is_feats, , drop = FALSE]
se_no_is <- se[setdiff(rownames(se), is_feats), , drop = FALSE]
log_msg("INFO", "Removed ", nrow(se_is), " IS from the SE, ", nrow(se_no_is), " features are kept.")

# ── PREPARE THE IS FOR THE CV ───────────────────────────────────────────────
log_section("Preparing data for MetaProViz::pool_estimation")

assay_df       <- t(assay(se_is, 1)) %>% as.data.frame(check.names = FALSE)
sample_info_df <- as.data.frame(colData(se_is), check.names = FALSE)

# For pool_estimation: all samples are "pools" (we want CV across everything).
# Preserve the original condition column (name set via --conditions_col) under
# a fixed internal name, since the "Conditions" column itself gets overwritten
# below for the pool_estimation() call.
no_condition_ids <- character(0)
if (conditions_col %in% colnames(sample_info_df)) {
  sample_info_df\$Conditions_original <- as.character(sample_info_df[[conditions_col]])
  # Samples without a value (NA or empty) are grouped as one condition of their own, so they stay in the
  # per-condition table and the box plots instead of breaking them.
  no_condition <- is.na(sample_info_df\$Conditions_original) | trimws(sample_info_df\$Conditions_original) == ""
  if (any(no_condition)) {
    no_condition_ids <- rownames(sample_info_df)[no_condition]
    sample_info_df\$Conditions_original[no_condition] <- "NA"
    log_msg("WARN", length(no_condition_ids), " sample(s) have no value in the ", conditions_col, " column: ",
            paste(no_condition_ids, collapse = ", "), ". They are grouped as the condition 'NA'.")
  }
}
sample_info_df\$Conditions <- "Pool"

# ── RUN pool_estimation ─────────────────────────────────────────────────────
log_section("Working out the CV of each IS")

# MetaProViz::pool_estimation() runs a PCA over the IS. It fails if only one IS is left and for an IS without values
# or without variation, and that would stop the whole run. These cases are handled here with the same formulas.
# Like MetaProViz, a value of 0 counts as not detected (same as NA), checked: 3 values set to 0 or to NA give the
# same CV and missing percentage in pool_estimation().
#   not detected   all values NA or 0: no CV exists, the IS is reported as not detected
#   one value      a single value only: no CV can be calculated
#   identical      the same value in every sample: CV is 0
#   usable         everything else goes to MetaProViz::pool_estimation(), as before
assay_nz        <- assay_df
assay_nz[is.na(assay_nz) | assay_nz == 0] <- NA
is_values       <- lapply(assay_nz, function(x) x[!is.na(x)])
is_not_detected <- names(Filter(function(v) length(v) == 0, is_values))
is_one_value    <- names(Filter(function(v) length(v) == 1, is_values))
is_identical    <- names(Filter(function(v) length(v) > 1 && sd(v) == 0, is_values))
is_usable       <- setdiff(colnames(assay_df), c(is_not_detected, is_one_value, is_identical))

# The sample that has the single value of an IS, named in the log and the report.
one_value_info <- vapply(is_one_value, function(f)
  paste0("sample ", rownames(assay_nz)[!is.na(assay_nz[[f]])]), character(1))

own_cv <- function(ids) {
  data.frame(
    Metabolite = ids,
    CV = vapply(ids, function(f) {
      v <- assay_nz[[f]]; v <- v[!is.na(v)]
      if (length(v) >= 2 && mean(v) != 0) 100 * sd(v) / mean(v) else NA_real_
    }, numeric(1)),
    MissingValuepercentage = vapply(ids, function(f) 100 * sum(is.na(assay_nz[[f]])) / nrow(assay_nz), numeric(1)),
    stringsAsFactors = FALSE, row.names = NULL
  ) %>%
    dplyr::mutate(HighVar = CV > cutoff_cv) %>%
    dplyr::select(Metabolite, CV, HighVar, MissingValuepercentage)
}

pe <- NULL
plots <- list()
if (length(is_usable) >= 2) {
  pe <- tryCatch(
    MetaProViz::pool_estimation(
      data            = assay_df[, is_usable, drop = FALSE],
      metadata_sample = sample_info_df,
      metadata_info   = c(PoolSamples = "Pool", Conditions = "Conditions"),
      cutoff_cv       = cutoff_cv,
      print_plot      = FALSE,
      save_plot       = NULL
    ),
    error = function(e) {
      log_msg("WARN", "pool_estimation() failed on the IS (", conditionMessage(e),
              "). Therefore the CV is calculated in the module instead.")
      NULL
    }
  )
} else {
  log_msg("INFO", "MetaProViz::pool_estimation() needs at least two IS with varying values. Therefore the CV is calculated ",
          "in the module.")
}

if (!is.null(pe)) {
  cv_metaproviz <- pe[["DF"]][["CV"]]
  if (is.null(cv_metaproviz)) abort("pool_estimation did not return DF\$CV.")
  plots   <- pe[["Plot"]]; if (is.null(plots)) plots <- list()
  own_ids <- setdiff(colnames(assay_df), is_usable)
  cv_df   <- if (length(own_ids) > 0) dplyr::bind_rows(cv_metaproviz, own_cv(own_ids)) else cv_metaproviz
} else {
  cv_df <- own_cv(colnames(assay_df))
}
cv_df <- as.data.frame(cv_df)[match(colnames(assay_df), cv_df\$Metabolite), , drop = FALSE]
rownames(cv_df) <- NULL

# Each case has its own message. If it applies to every IS, the message says so for all of them. Every message is a
# warning in the log (all names) and a bold sentence under the run summary in the report (at most 25 names, then a
# pointer to the log).
case_msgs <- function(cap_n) {
  more <- function(x) if (length(x) > cap_n) paste0(" and ", length(x) - cap_n, " more (see ", log_out, ")") else ""
  nm   <- function(x) paste0(paste(head(x, cap_n), collapse = ", "), more(x))
  them <- function(x) if (length(x) == 1) "it" else "them"
  nmi  <- function(x, info) paste0(paste0(head(x, cap_n), " (", head(info, cap_n), ")", collapse = ", "), more(x))
  n_is   <- length(is_feats)
  all_is <- if (n_is == 1) "The IS" else paste0("All ", n_is, " IS")
  has_v  <- if (n_is == 1) "has" else "have"
  m <- character(0)
  if (ncol(se) == 1) {
    m <- c(m, "The input has only one sample. No CV can be calculated.")
  } else if (length(is_one_value) == n_is) {
    m <- c(m, paste0(all_is, " ", has_v, " only one value across all samples: ", nmi(is_one_value, one_value_info),
      ". A CV needs at least two values. Therefore no CV can be calculated. Please check how many samples have a value for the IS."))
  }
  if (length(is_not_detected) == n_is) {
    m <- c(m, paste0(all_is, " ", has_v, " no value in any sample (all values missing or 0). No CV can be calculated. ",
      "Please check the spike-in, the extraction and the measurement, and whether the IS detection selected the right features."))
  }
  if (length(is_identical) == n_is) {
    m <- c(m, paste0(all_is, " ", has_v, " the identical value in every sample: ", nm(is_identical),
      ". The CV is 0", if (n_is == 1) "" else " for all of them", ". Exactly identical measured values are very unlikely. ",
      "Please check whether they were filled in or changed earlier."))
  }
  if (all(is.na(cv_df\$CV)) && length(m) == 0) {
    m <- c(m, paste0("No CV can be calculated for any of the IS (", length(is_not_detected), " not detected, ",
      length(is_one_value), " with one value only). Therefore no IS can be compared with the CV cutoff."))
  }
  if (length(is_not_detected) > 0 && length(is_not_detected) < n_is)
    m <- c(m, paste0("IS not detected in any sample (all values missing or 0): ", nm(is_not_detected),
                     ". No CV can be calculated for ", them(is_not_detected), "."))
  if (length(is_one_value) > 0 && length(is_one_value) < n_is && ncol(se) > 1)
    m <- c(m, paste0("IS with a single value only: ", nmi(is_one_value, one_value_info), ". No CV can be calculated for ", them(is_one_value), "."))
  if (length(is_identical) > 0 && length(is_identical) < n_is)
    m <- c(m, paste0("IS with identical values in every sample: ", nm(is_identical),
                     ". Exactly identical measured values are very unlikely. Please check whether they were filled in or changed earlier."))
  m
}
all_case_txt    <- case_msgs(.Machine\$integer.max)
for (m in all_case_txt) log_msg("WARN", m)
if (length(is_feats) == 1)
  log_msg("INFO", "With only one IS, a high CV cannot show whether the problem is specific to this compound or affects the whole run.")

# ── SUPERPLOT: ONE BOX PER INTERNAL STANDARD ────────────────────────────────
# viz_superplot()'s `Conditions` role becomes the x-axis grouping (one box per
# distinct value there) — we want exactly one box per standard, so every
# sample gets the same dummy value here. The real experimental condition (if
# present) goes into the `Superplot` role instead, which only colors the dots
log_section("Building internal-standard superplots")

superplot_meta <- sample_info_df
# Drop the "Conditions" column set earlier to "Pool" for the pool_estimation()
# call above — it's irrelevant here and its name collides with what
# viz_superplot() renames AllSamples to internally.
superplot_meta\$Conditions <- NULL
superplot_meta\$AllSamples <- "IS"
has_real_conditions <- "Conditions_original" %in% colnames(superplot_meta)
superplot_metadata_info <- if (has_real_conditions) {
  c(Conditions = "AllSamples", Superplot = "Conditions_original")
} else {
  c(Conditions = "AllSamples")
}

# zeros count as not detected, like NA, so they are not drawn as points at 0
superplot_data <- assay_nz[, setdiff(colnames(assay_nz), is_not_detected), drop = FALSE]
superplot_res <- if (ncol(superplot_data) == 0) NULL else tryCatch(
  MetaProViz::viz_superplot(
    data            = superplot_data,
    metadata_sample = superplot_meta,
    metadata_info   = superplot_metadata_info,
    plot_type       = "Box",
    print_plot      = FALSE,
    save_plot       = NULL
  ),
  error = function(e) {
    log_msg("WARN", "viz_superplot() failed: ", conditionMessage(e))
    NULL
  }
)
superplot_plots <- if (!is.null(superplot_res)) superplot_res[["Plot"]] else list()

# CV per standard, per condition — only meaningful with real condition
# metadata (with the dummy single-value column every condition would be
# identical to the overall CV already shown below). Computed once here so
# both the TSV output and the report table below reuse the same data.
condition_cv_df <- if (has_real_conditions) {
  long_df <- assay_nz %>%
    tibble::rownames_to_column("Sample") %>%
    tidyr::pivot_longer(-Sample, names_to = "Metabolite", values_to = "Value") %>%
    dplyr::mutate(Condition = sample_info_df[Sample, "Conditions_original"]) %>%
    dplyr::filter(!is.na(Value))

  long_df %>%
    dplyr::group_by(Metabolite, Condition) %>%
    dplyr::summarise(
      CV = if (dplyr::n() >= 2) 100 * sd(Value) / mean(Value) else NA_real_,
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(id_cols = Metabolite, names_from = Condition, values_from = CV)
} else {
  data.frame(Metabolite = character(0))
}
if (has_real_conditions) {
  absent <- setdiff(colnames(assay_df), condition_cv_df\$Metabolite)
  if (length(absent) > 0) condition_cv_df <- dplyr::bind_rows(condition_cv_df, tibble::tibble(Metabolite = absent))
  for (cn in setdiff(sort(unique(na.omit(as.character(sample_info_df\$Conditions_original)))), colnames(condition_cv_df)))
    condition_cv_df[[cn]] <- NA_real_
}

# ── WRITE OUTPUTS ──────────────────────────────────────────────────────────
log_section("Writing outputs")

high_var <- cv_df %>% filter(HighVar == TRUE) %>% pull(Metabolite)
log_msg("INFO", "IS with CV > ", cutoff_cv, "%: ",
        length(high_var), " / ", nrow(cv_df))

# SE as RDS, or the three matrices as TSV/CSV (only the first assay)
write_out <- function(x, files) {
  if (output_format == "se") return(saveRDS(x, files))
  flat <- function(d) {
    d <- as.data.frame(d, check.names = FALSE, stringsAsFactors = FALSE)
    d[] <- lapply(d, function(col) if (is.list(col)) vapply(col, paste, character(1), collapse = "; ") else col)
    d
  }
  sep <- if (output_format == "csv") "," else "\\t"
  for (i in 1:3) {
    tab <- switch(i,
      data.frame(Feature = rownames(x), flat(assay(x, 1)), check.names = FALSE),
      data.frame(Feature = rownames(x), flat(rowData(x)), check.names = FALSE),
      data.frame(Sample = colnames(x), flat(colData(x)), check.names = FALSE))
    write.table(tab, files[i], sep = sep, quote = output_format == "csv", row.names = FALSE, na = "NA")
  }
}
if (output_format != "se" && length(assays(se)) > 1)
  log_msg("WARN", "The input has ", length(assays(se)), " assays. Only the first one is written to the data matrix.")
write_out(se_no_is, se_out)
write.table(cv_df, cv_out,
            sep = "\\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(condition_cv_df, condition_cv_out,
            sep = "\\t", quote = FALSE, row.names = FALSE, na = "NA")
saveRDS(list(PoolEstimationStyle = plots, Superplot = superplot_plots),
        plots_out)
log_msg("INFO", "Written: ", paste(se_out, collapse = ", "), ", ", cv_out, ", ", condition_cv_out, ", ", plots_out)

# ── HTML REPORT ────────────────────────────────────────────────────────────
log_section("Building report")

png_to_b64 <- function(p, w = 8, h = 5) {
  tmp <- tempfile(fileext = ".png")
  ggplot2::ggsave(tmp, plot = p, width = w, height = h, dpi = 100,
                  units = "in", bg = "white")
  raw <- readBin(tmp, "raw", n = file.info(tmp)\$size); unlink(tmp)
  paste0("data:image/png;base64,", base64enc::base64encode(raw))
}

# One box plot per standard. Name + CV% go on the plot's own x-axis label.
superplot_html <- ""
if (length(superplot_plots) > 0) {
  cv_lookup <- setNames(cv_df\$CV, cv_df\$Metabolite)
  for (nm in names(superplot_plots)) {
    p <- superplot_plots[[nm]]
    if (!inherits(p, "ggplot")) next
    cv_val <- cv_lookup[[nm]]
    cv_str <- if (!is.null(cv_val) && !is.na(cv_val)) format(round(cv_val, 1)) else NA
    axis_label <- if (!is.na(cv_str)) sprintf("%s (%s%%)", nm, cv_str) else nm
    plot_title <- if (!is.na(cv_str)) {
      sprintf("IS: %s (%s%%)", nm, cv_str)
    } else {
      sprintf("IS: %s", nm)
    }
    p <- p +
      ggplot2::labs(title = plot_title, subtitle = NULL, color = "Condition") +
      ggplot2::scale_x_discrete(labels = axis_label) +
      ggplot2::xlab(NULL) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))
    src <- tryCatch(png_to_b64(p),
                    error = function(e) {
                      log_msg("WARN", "Could not render superplot for '", nm, "': ",
                              conditionMessage(e)); NA_character_
                    })
    if (!is.na(src)) {
      superplot_html <- paste0(
        superplot_html,
        sprintf('<img src="%s" alt="%s" />\\n', src, esc(axis_label))
      )
    }
  }
}

# CV per standard, per condition — only meaningful with real condition
# metadata (with the dummy single-value column every condition would be
# identical to the overall CV already shown above). condition_cv_df was
# already computed above, before "WRITE OUTPUTS", and saved as
# <prefix>.condition_cv.tsv; reused here for the report table.
condition_cv_section_html <- ""
if (has_real_conditions) {
  cond_cols <- setdiff(colnames(condition_cv_df), "Metabolite")
  cond_header <- paste0("<th>", esc(c("IS", cond_cols)), "</th>", collapse = "")
  cond_rows <- apply(condition_cv_df, 1, function(r) {
    cells <- vapply(r[-1], function(v) {
      if (is.na(v)) "—" else sprintf("%s%%", format(round(as.numeric(v), 1)))
    }, character(1))
    paste0("<tr><td>", esc(r[["Metabolite"]]), "</td>",
           paste0("<td>", cells, "</td>", collapse = ""), "</tr>")
  })

  condition_cv_section_html <- sprintf(
    '<h2>CV per IS, per condition</h2>
<p>
Expectation: a low CV in every condition, since the IS are added in the same amount to every
sample. A high CV in all conditions points to a general issue of the run (e.g. instrument
drift). A high CV in only some conditions is more concerning. It points to an issue specific to
those samples (e.g. a batch effect from how they were handled), since the IS values should
not depend on the condition. The table has the CV of each IS in each condition, with one row per IS and one column per
condition, and is saved as <code>%s</code>.
</p>
%s<div class="scroll-table"><table><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>',
    esc(condition_cv_out),
    if (length(no_condition_ids) > 0) paste0(
      "<p><strong>", length(no_condition_ids), " sample", if (length(no_condition_ids) == 1) " has" else "s have",
      " no value in the ", esc(conditions_col), " column and ", if (length(no_condition_ids) == 1) "is" else "are",
      " grouped as the condition 'NA':</strong> ",
      names_html(no_condition_ids), ".</p>\\n") else "",
    cond_header, paste(cond_rows, collapse = "")
  )
}

high_var_html <- if (length(high_var) > 0) {
  paste0("<ul>", paste0("<li>", esc(high_var), "</li>", collapse = ""),
         "</ul>")
} else if (all(is.na(cv_df\$CV))) {
  "<p><em>No IS has a CV. Therefore none can be compared with the CV cutoff.</em></p>"
} else paste0("<p><em>None — every IS ", if (any(is.na(cv_df\$CV))) "with a CV " else "",
              "is below the CV cutoff.</em></p>")

# Render the actual matched IS feature names as a readable comma-separated
# list (not a regex). Lists the names of the standards, not pattern metacharacters.
is_names_html <- if (length(is_feats) > 0) {
  paste(esc(is_feats), collapse = ", ")
} else "—"

# CV and MissingValuepercentage columns (both in percent): round + add "%". The saved TSV keeps the plain numbers.
cv_rows <- apply(cv_df, 1, function(r) {
  vals <- vapply(names(r), function(col) {
    if (col %in% c("CV", "MissingValuepercentage")) {
      v <- suppressWarnings(as.numeric(r[[col]]))
      if (!is.na(v)) sprintf("%s%%", format(round(v, 1))) else "—"
    } else {
      if (is.na(r[[col]])) "—" else as.character(r[[col]])
    }
  }, character(1))
  paste0("<tr>", paste0("<td>", esc(vals), "</td>", collapse = ""), "</tr>")
})
cv_table <- paste0(
  '<div class="scroll-table"><table><thead><tr>',
  paste0("<th>", esc(sub("^Metabolite\$", "IS", colnames(cv_df))), "</th>", collapse = ""),
  '</tr></thead><tbody>',
  paste(cv_rows, collapse = ""),
  '</tbody></table></div>')

# IS with a problem (not detected, one value only, identical values): one section with the names and the advice. The
# numbers stay in the CV table (a dash if no CV exists).
unusable_html <- ""
if (length(is_not_detected) + length(is_one_value) + length(is_identical) > 0) {
  unusable_html <- paste0(
    "<h2>IS to check</h2>\\n",
    if (length(is_not_detected) > 0) paste0(
      "<p><strong>Not detected in any sample (all values missing or 0):</strong> ", names_html(is_not_detected),
      ". No CV can be calculated. This points to a problem with the spike-in, the extraction or the measurement of this compound.</p>\\n") else "",
    if (length(is_one_value) > 0) paste0(
      "<p><strong>Only one value in all samples:</strong> ", info_html(is_one_value, one_value_info),
      ". A CV needs at least two values.</p>\\n") else "",
    if (length(is_identical) > 0) paste0(
      "<p><strong>Identical value in every sample:</strong> ", names_html(is_identical),
      ". The CV is 0. Exactly identical measured values are very unlikely. Please check whether they were filled in or changed earlier.</p>\\n") else "")
}
one_is_html <- if (length(is_feats) == 1) {
  "<p>With only one IS, a high CV cannot show whether the problem is specific to this compound or affects the whole run.</p>"
} else ""
is_extra_rows <- if (length(is_not_detected) > 0) {
  sprintf("\\n   <dt>IS not detected</dt><dd>%d</dd>", length(is_not_detected))
} else ""

# Box plots exist only for IS with values, say why for the others.
superplot_part <- if (nzchar(superplot_html)) superplot_html else if (length(is_not_detected) == length(is_feats)) {
  "<p><strong>No box plots, since no IS has values.</strong></p>"
} else "<p><em>No box plots returned.</em></p>"
if (nzchar(superplot_html) && length(is_not_detected) > 0) {
  superplot_part <- paste0(superplot_part, "<p><strong>No box plot for </strong>", names_html(is_not_detected),
                           "<strong>, since ", if (length(is_not_detected) == 1) "it has" else "they have", " no values.</strong></p>")
}

# Result of the database ID check of the auto-detected names.
is_base_html <- paste0(is_label_html, if (length(is_attempted) > 0) paste0(" (fallback, since no IS features were found by the ",
                                                                                is_attempted_html, ")") else "")
strategy_txt <- is_base_html
# IS with an entry in the database ID columns: a bold sentence under the run summary.
id_flag_note_html <- if (length(id_flagged) > 0) paste0(
  "<p><strong>IS with an entry in the database ID columns of the feature metadata: ",
  paste0("<code>", esc(id_flagged), "</code> (", vapply(id_found_cols, function(x) esc(paste(x, collapse = ", ")), character(1)), ")",
         collapse = ", "),
  ". ", if (length(id_flagged) == 1) "It is" else "They are", " used as IS. Please check that ",
  if (length(id_flagged) == 1) "it is a real IS" else "they are real IS", ".</strong></p>\\n") else ""

# Result of the IS removal for the run summary and the description.
files_html <- function(f) if (length(f) == 1) paste0("<code>", esc(f), "</code>") else
  paste0(paste0("<code>", esc(head(f, -1)), "</code>", collapse = ", "), " and <code>", esc(tail(f, 1)), "</code>")
outputs_note_html <- sprintf(
  "<p>The IS are not part of the output %s %s, since they are spike-ins added for this quality control and not measured metabolites of the samples. The removed IS are the rows of the CV table <code>%s</code>.</p>",
  if (length(se_out) == 1) "file" else "files", files_html(se_out), esc(cv_out))
# One bold sentence under the run summary if any IS has a problem, the details are in the section "IS to check".
problem_counts <- c(length(is_not_detected), length(is_one_value), length(is_identical))
problem_one   <- c("not detected in any sample", "only one value in all samples", "identical values in every sample")
problem_many  <- c("not detected", "with one value only", "with identical values")
keep <- problem_counts > 0
all_case_html <- if (ncol(se) == 1) {
  "<div class=\"notice\"><strong>The input has only one sample. No CV can be calculated. See \\"IS to check\\" below.</strong></div>\\n"
} else if (any(keep)) {
  n_problem <- sum(problem_counts)
  who  <- if (length(is_feats) == 1) "the IS" else if (n_problem == length(is_feats)) paste0("all ", length(is_feats), " IS") else
    paste0(n_problem, " of ", length(is_feats), " IS")
  what <- if (sum(keep) == 1) problem_one[keep] else paste0(problem_counts[keep], " ", problem_many[keep], collapse = ", ")
  flag <- paste0("Issue with ", who, ": ", what, ".",
    if (all(is.na(cv_df\$CV))) " No CV can be calculated." else "",
    " See \\"IS to check\\" below.")
  paste0("<div class=\"notice\"><strong>", esc(flag), "</strong></div>\\n")
} else ""

n_features_txt <- sprintf("%d feature%s after %d IS %s removed", nrow(se_no_is), if (nrow(se_no_is) == 1) "" else "s",
                          length(is_feats), if (length(is_feats) == 1) "was" else "were")
n_samples_txt  <- paste0(ncol(se), " sample", if (ncol(se) == 1) "" else "s")
run_meta_html <- sprintf(
'<dl>
   <dt>Input source</dt><dd>%s</dd>
   <dt>n features</dt><dd>%s</dd>
   <dt>n samples</dt><dd>%s</dd>
   <dt>IS detection strategy</dt><dd>%s</dd>
   <dt>IS used</dt><dd>%s</dd>
   <dt>CV cutoff</dt><dd>%s%%</dd>
   <dt>IS above cutoff</dt><dd>%d / %d</dd>%s
 </dl>',
 esc(source_label),
 n_features_txt, n_samples_txt,
 strategy_txt, is_names_html, format(cutoff_cv),
 length(high_var), nrow(cv_df), is_extra_rows)

# ── SESSION SUMMARY (packages used + citation reminder) ─────────────────────
session_info_text <- paste(capture.output(sessionInfo()), collapse = "\\n")
session_summary_html <- sprintf(
'<details>
  <summary>Full R session info (click to expand)</summary>
  <pre>%s</pre>
</details>
<p>Please cite MetaProViz, base R and the dependencies listed here and in the Citations section
if you use results or plots from this report.</p>',
esc(session_info_text))

html <- sprintf('<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Internal Standard (IS) Report</title>
<style>
  body{font-family:Arial,sans-serif;max-width:980px;margin:2em auto;padding:0 1em;line-height:1.5}
  h1{color:#1F3A5F;border-bottom:2px solid #1F3A5F;padding-bottom:.3em}
  h2{color:#1F3A5F;margin-top:2em}
  h3{color:#345;margin-top:1.4em}
  dl{display:grid;grid-template-columns:max-content 1fr;gap:.3em 1em}
  dt{font-weight:bold;color:#345}
  table{border-collapse:collapse;width:100%%;font-size:.9em}
  th,td{border:1px solid #ccc;padding:.35em .6em;text-align:left}
  th{background:#f0f4f8}
  img{max-width:100%%;border:1px solid #ddd;padding:4px;background:#fff;margin:.4em 0 1.2em}
  code{background:#f2f2f2;padding:.15em .35em;border-radius:3px}
  .notice{background:#f0f0f0;border-left:4px solid #888;border-radius:4px;padding:.7em 1em;margin:.9em 0}
  li{margin:.25em 0}
  .scroll-table{max-height:320px;overflow-y:auto;border:1px solid #cbd6e4;border-radius:6px}
  .scroll-table table{border:0;margin:0}
  .scroll-table thead th{position:sticky;top:0;background:#f0f4f8;z-index:1}
  details{margin:.6em 0}
  summary{cursor:pointer;color:#1F3A5F;font-weight:bold}
  details pre{max-height:400px;overflow:auto;font-size:.8em;background:#f7f9fb;border:1px solid #cbd6e4;border-radius:6px;padding:.8em}
</style></head><body>
<h1>Internal Standard (IS) Report</h1>
<p><em>Generated by the METAPROVIZ_INTERNALSTANDARD module.</em></p>

<h2>Run summary</h2>%s
%s%s
<h2>Citations</h2>
<p>
If you use results or plots from this analysis in a publication, please cite:
</p>
<ul>
  <li><strong>MetaProViz:</strong>
  <a href="https://doi.org/10.1038/s44320-026-00231-8">Schmidt et al., Integrated metabolomics data analysis to generate mechanistic hypotheses with MetaProViz, Molecular Systems Biology 2026.</a></li>
  <li><strong>Dependencies:</strong>
  <a href="https://doi.org/10.32614/R.manuals">R Core Team, R: A Language and Environment for Statistical Computing, R Foundation for Statistical Computing</a>, <a href="https://doi.org/10.18129/B9.bioc.SummarizedExperiment">Morgan et al., SummarizedExperiment: A container (S4 class) for matrix-like assays, Bioconductor</a>, <a href="https://doi.org/10.1038/nmeth.3252">Huber et al., Orchestrating high-throughput genomic analysis with Bioconductor, Nature Methods 2015</a> and <a href="https://doi.org/10.1007/978-3-319-24277-4">Wickham, ggplot2: Elegant Graphics for Data Analysis, Springer-Verlag New York 2016</a>.
  The other R packages used are listed with their versions in the session summary at the end of this report.</li>
  <li><strong>This module:</strong> <code>metaproviz/internalstandard</code>, an nf-core module. It has no publication of its own and is cited through the nf-core publication below. In the methods, name the module.</li>
  <li><strong>nf-core:</strong> <a href="https://doi.org/10.1038/s41587-020-0439-x">Ewels et al., The nf-core framework for community-curated bioinformatics pipelines, Nature Biotechnology 2020.</a></li>
  <li><strong>Nextflow:</strong> <a href="https://doi.org/10.1038/nbt.3820">Di Tommaso et al., Nextflow enables reproducible computational workflows, Nature Biotechnology 2017.</a></li>
</ul>

<h2>Description</h2>
<p>
This report evaluates IS (spike-in) features to assess the technical noise of the measurement,
independent of biology. IS are added in a fixed, known amount to every sample
before measurement. Unlike real metabolites, their values are expected to stay constant
across the run.
</p>
<p>
This module computes the coefficient of variation (CV) of each IS across all samples.
Pool and QC samples must be removed before this module (e.g. with <code>metaproviz/poolestimation</code>),
since their already-pooled nature would understate the true injection-to-injection variability. A high CV
indicates instrument drift or injection issues rather than biological variation.
</p>
<p>
CV is calculated as:
</p>
<p><code>CV(%%) = 100 × sd(x1, x2, ..., xn) / mean(x1, x2, ..., xn)</code></p>
<p>
<code>MetaProViz::pool_estimation()</code> uses the same formula, so this module reuses that function. The function
is made for pool samples, which are aliquots of the same mixture that are measured several times. A feature should
have the same value in all of them. An IS is added in the same amount to all samples and should therefore give the same
value in each sample. For this reason all samples are temporarily labelled as pool samples for the calculation.
In some cases the function cannot be used (a single IS, or an IS with identical values in every sample). Then the module
calculates the CV directly with the formula. An IS with no value or only one value has no CV and is listed under
"IS to check".
</p>
<p>
If condition metadata is available, the same CV calculation is additionally broken down per
condition (see the "CV per IS, per condition" table further down), and the box plots
below also color individual dots by condition.
</p>
%s

<h2>Plots</h2>
<p>
Each IS gets one box plot below, generated with MetaProViz\\'s
<code>viz_superplot()</code> function, showing its values across all samples used for the CV
calculation.
</p>
<p>
What to expect: dots of different colors should be interspersed throughout the box, not
clustered by color, since the IS is expected to behave the same regardless of
condition. If dots of one color cluster together at their own distinct position, that points to
a systematic shift for that specific condition (e.g. a batch effect). Without condition
metadata, all dots are shown in black instead.
</p>
%s
%s

<h2>CV per IS, overall</h2>
<p>
The CV of every IS, calculated across all samples together. The table has one row per IS.
<code>CV</code> is given in percent. <code>HighVar</code> shows whether the CV is above the cutoff of <code>%s%%</code>.
<code>MissingValuepercentage</code> is the percentage of samples in which the IS was not detected, which means that its
value is missing or 0. The table is saved as <code>%s</code>.
</p>
%s
%s
%s<h2>IS above the CV cutoff (%s%%)</h2>
%s

<h2>Run notes (log)</h2>
<p>A separate file, <code>%s</code>, records every step of this run,
including any warnings. Check it first if a result here looks unexpected.</p>

<h2>Session summary</h2>
<p>Full R session used to generate this report — every package, and its version:</p>
%s

</body></html>',
run_meta_html,
all_case_html,
id_flag_note_html,
outputs_note_html,
superplot_part,
condition_cv_section_html,
format(cutoff_cv), esc(cv_out),
cv_table, one_is_html, unusable_html, format(cutoff_cv), high_var_html,
esc(log_out), session_summary_html)

writeLines(html, report_out)
log_msg("INFO", "Written: ", report_out)

# Cleaned up not needed files
unlink("MetaProViz_Results", recursive = TRUE)
unlink("metaproviz-log", recursive = TRUE)
unlink("omnipathr-log", recursive = TRUE)
unlink("Rplots.pdf")

flush_log(log_out)
message(sprintf(
  "\\n✓ internal_standard complete — %d / %d IS flagged HighVar | %d warning(s)",
  length(high_var), nrow(cv_df), wc))
