#!/usr/bin/env Rscript
# =============================================================================
# metaproviz_processingimpute.R: missing value imputation step of the MetaProViz processing chain
#
# Approach: MetaProViz::mvi_imputation() replaces each missing value (NA or 0) with a percentage of the minimum of that
# feature in its condition (--mvi_percentage, default 50%, the half-minimum, Wei et al., 2018,
# https://doi.org/10.1038/s41598-017-19120-0). The Conditions column is required.
#
# Features removed from the output: the ones flagged TRUE in --filteredfeatures_col (unless --filter_flagged FALSE,
# --keep_features keeps single ones) and the ones with no measured value in a condition (unless --remove_unimputable
# FALSE), which cannot be imputed and would get value NA instead.
#
# Pool/QC samples must be removed before this module, they are not detected here.
#
# Outputs, all named <prefix>.<suffix>: the SE (or three tables if output_format is tsv or csv), imputed_values.tsv,
# removed_features.tsv, the HTML report and the log.
# =============================================================================


set.seed(42)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(magrittr)
  library(SummarizedExperiment); library(S4Vectors)
})

suppressPackageStartupMessages(library(MetaProViz))

# ── PARSE PARAMETERS FROM NEXTFLOW ──────────────────────────────────────────
as_bool <- function(x) toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES")
# Empty optional inputs arrive as "", "null" or "[]".
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

# Defaults of all options, overwritten by the ones in task.ext.args.
opt <- list(
  conditions_col       = "Conditions",
  biorep_col           = "Biological_Replicates",
  mvi_percentage       = "50",
  filter_flagged       = "TRUE",
  filteredfeatures_col = "FilteredFeaturesFlag",
  remove_unimputable   = "TRUE",
  keep_features        = character(0)
)
# R raw string: quotes and apostrophes inside the arguments (feature names like 5'-AMP) cannot break the R code.
args_opt <- parse_args(r"---($task.ext.args)---", list_opts = "keep_features")
invalid_opts <- setdiff(names(args_opt), names(opt))   # reported below, once the log exists
for (ao in setdiff(names(args_opt), invalid_opts)) opt[[ao]] <- args_opt[[ao]]

# The inputs of the process, empty for the shape that is not used.
se_path             <- nullify('$se_rds')
data_matrix_path    <- nullify('$data_matrix')
feature_matrix_path <- nullify('$feature_matrix')
sample_matrix_path  <- nullify('$sample_matrix')
conditions_col <- opt\$conditions_col
biorep_col     <- opt\$biorep_col
mvi_percentage <- as.numeric(opt\$mvi_percentage)
filter_flagged       <- as_bool(opt\$filter_flagged)
filteredfeatures_col <- opt\$filteredfeatures_col
remove_unimputable   <- as_bool(opt\$remove_unimputable)
keep_asked           <- trimws(opt\$keep_features); keep_asked <- keep_asked[nzchar(keep_asked)]

# Output files are named <prefix>.<suffix>.
prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
# output_format is a val input of the module. se: the SE as one RDS. tsv/csv: three files (data, feature and sample
# matrix), the same shape as the input matrices. An empty value means se.
output_format_in   <- trimws('$output_format')
output_format      <- if (output_format_in %in% c("", "null", "[]")) "se" else tolower(output_format_in)
se_out             <- if (output_format == "se") paste0(prefix, ".imputed.rds") else
  paste0(prefix, ".imputed", c(".data_matrix.", ".feature_matrix.", ".sample_matrix."), output_format)
values_out         <- paste0(prefix, ".imputed_values.tsv")
removed_out        <- paste0(prefix, ".removed_features.tsv")
report_out         <- paste0(prefix, ".report.html")
log_out            <- paste0(prefix, ".log")

# ── LOGGING ─────────────────────────────────────────────────────────────────
log_lines <- character(0); warn_count <- 0; err_count <- 0
log_msg <- function(level = "INFO", ...) {
  txt <- paste0("[", level, "] ", paste(..., sep = ""))
  log_lines <<- c(log_lines, txt); message(txt)
  if (level == "WARN")  warn_count <<- warn_count + 1
  if (level == "ERROR") err_count  <<- err_count  + 1
}
log_section <- function(title) {
  sep <- paste(rep("─", 60), collapse = "")
  log_lines <<- c(log_lines, "", sep, paste0("  ", title), sep)
  message(sep, "\\n  ", title, "\\n", sep)
}
flush_log <- function(path = log_out) {
  writeLines(c(log_lines, "",
    sprintf("SUMMARY: %d warning(s), %d error(s)", warn_count, err_count)), path)
}
# abort() writes the error to the log, saves the log and stops the run.
abort <- function(...) { log_msg("ERROR", ...); flush_log(); stop(paste0(...), call. = FALSE) }
# Escapes text for the HTML report.
esc <- function(x) gsub(">", "&gt;", gsub("<", "&lt;",
                    gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE), fixed = TRUE)

if (length(invalid_opts) > 0)
  abort("Invalid option", if (length(invalid_opts) > 1) "s" else "", ": ", paste0("--", invalid_opts, collapse = ", "),
        ". The options are: ", paste0("--", names(opt), collapse = ", "), ".")
if (!output_format %in% c("se", "tsv", "csv"))
  abort("output_format must be se, tsv or csv (empty means se), not '", output_format_in, "'.")

# ── METAPROVIZ VERSION ──────────────────────────────────────────────────────
# mvi_imputation() is part of MetaProViz from version 4.1.
if (!"mvi_imputation" %in% getNamespaceExports("MetaProViz"))
  abort("This module needs MetaProViz >= 4.1 (exports mvi_imputation()), found ",
        as.character(packageVersion("MetaProViz")), ".")

# ── LOAD INPUT ──────────────────────────────────────────────────────────────
log_section("Loading input data")

# Two acceptable input shapes, mutually exclusive: --se (a SummarizedExperiment
# .rds), or all three of --data_matrix/--feature_matrix/--sample_matrix (plain
# TSV/CSV, language-independent).
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

# how the report names the feature and sample metadata, depending on how the input was given
feat_meta_txt <- if (have_se) {
  sprintf("feature metadata (<code>rowData</code> of <code>%s</code>)", esc(basename(se_path)))
} else sprintf("feature metadata (<code>%s</code>)", esc(basename(feature_matrix_path)))
samp_meta_txt <- if (have_se) {
  sprintf("sample metadata (<code>colData</code> of <code>%s</code>)", esc(basename(se_path)))
} else sprintf("sample metadata (<code>%s</code>)", esc(basename(sample_matrix_path)))

if (have_se) {
  if (!file.exists(se_path)) abort("SE file not found: ", se_path)
  se <- tryCatch(readRDS(se_path),
                 error = function(e) abort("Failed to read SE: ",
                                           conditionMessage(e)))
  if (!is(se, "SummarizedExperiment"))
    abort("Object in ", se_path, " is not a SummarizedExperiment.")
  log_msg("INFO", "Loaded SE from ", basename(se_path), ": ", nrow(se),
          " features × ", ncol(se), " samples.")
} else {
  for (p in c(data_matrix_path, feature_matrix_path, sample_matrix_path))
    if (!file.exists(p)) abort("Input file not found: ", p)

  # TSV or CSV, by the file extension.
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

  # Build the SE from the three tables, matched by the IDs in their first column.
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

# Name of the input for the report.
source_label <- if (have_se) {
  basename(se_path)
} else {
  paste(basename(data_matrix_path), basename(feature_matrix_path),
        basename(sample_matrix_path), sep = " + ")
}


# ── FLAGGED FEATURES ────────────────────────────────────────────────────────
# The filter module only flags features. With --filter_flagged TRUE the features that are TRUE are removed here.
log_section("Flagged features")
flagged_ids   <- character(0)
kept_flagged_ids <- character(0)   # flagged TRUE, but kept because they are in --keep_features
keep_not_found   <- setdiff(keep_asked, rownames(se))
if (length(keep_not_found) > 0)
  log_msg("WARN", "The list --keep_features contains feature(s) that are not in the input: ", paste(keep_not_found, collapse = ", "), ".")
flag_col_found <- filteredfeatures_col %in% colnames(rowData(se))
if (filter_flagged && flag_col_found) {
  # The flags can be TRUE/FALSE, T/F, 1/0 or yes/no.
  raw_flags <- rowData(se)[[filteredfeatures_col]]
  flags_in  <- if (is.logical(raw_flags)) raw_flags else toupper(trimws(as.character(raw_flags))) %in% c("TRUE", "T", "1", "YES")
  flagged_ids <- rownames(se)[flags_in %in% TRUE]
  kept_flagged_ids <- intersect(flagged_ids, keep_asked)
  flagged_ids      <- setdiff(flagged_ids, keep_asked)
  if (length(kept_flagged_ids) > 0)
    log_msg("INFO", "Kept ", length(kept_flagged_ids), " feature(s) flagged TRUE in '", filteredfeatures_col,
            "', since they are in --keep_features: ", paste(kept_flagged_ids, collapse = ", "),
            ". The flag stays TRUE in the output SE.")
  if (length(flagged_ids) >= nrow(se))
    abort("All ", nrow(se), " features are flagged TRUE in the column '", filteredfeatures_col, "'. No feature would be left ",
          "to impute. Please check the filtering parameters of the metaproviz/processingfilter module.")
  if (length(flagged_ids) > 0) {
    se <- se[setdiff(rownames(se), flagged_ids), , drop = FALSE]
    log_msg("INFO", "Removed ", length(flagged_ids), " feature(s) flagged TRUE in '", filteredfeatures_col, "': ",
            paste(flagged_ids, collapse = ", "), ".")
  } else {
    log_msg("INFO", "No feature is flagged TRUE in '", filteredfeatures_col, "', all ", nrow(se), " features are kept.")
  }
} else if (filter_flagged) {
  log_msg("WARN", "The feature metadata has no column '", filteredfeatures_col, "', and no feature is removed. All ",
          nrow(se), " features are imputed.")
} else {
  log_msg("INFO", "--filter_flagged FALSE, all ", nrow(se), " features are kept.")
}

# A missing value is NA or 0.
missing_mask     <- is.na(assay(se, 1)) | assay(se, 1) == 0

# ── IMPUTATION ────────────────────────────────────────────────────────────
log_section("Running MetaProViz::mvi_imputation()")

# mvi_imputation() takes the samples in rows and the features in columns.
assay_df       <- t(assay(se, 1)) %>% as.data.frame(check.names = FALSE)
sample_info_df <- as.data.frame(colData(se), check.names = FALSE)
if (!conditions_col %in% colnames(sample_info_df))
  abort("Conditions column '", conditions_col, "' not found in colData. ",
        "Imputation takes the minimum per condition, so conditions are required: ",
        "set --conditions_col to a colData column that holds them.")

# Samples with NA, an empty value or only whitespace in the Conditions column are grouped as one condition called "NA",
# and a warning lists them. Numbers in the column are turned into text, so they are conditions like any other.
cond_chr <- as.character(sample_info_df[[conditions_col]])
no_cond  <- is.na(cond_chr) | trimws(cond_chr) == ""
no_condition_ids <- rownames(sample_info_df)[no_cond]
if (any(no_cond)) {
  cond_chr[no_cond] <- "NA"
  log_msg("WARN", length(no_condition_ids), " sample(s) have no value in the Conditions column '", conditions_col, "': ",
          paste(no_condition_ids, collapse = ", "), ". They are treated as the condition 'NA'.")
}
sample_info_df[[conditions_col]] <- cond_chr

# The conditions are required, the biological replicates are optional.
metadata_info <- c(Conditions = conditions_col)
if (biorep_col %in% colnames(sample_info_df))
  metadata_info <- c(metadata_info, Biological_Replicates = biorep_col)

imputed_df <- tryCatch(
  MetaProViz::mvi_imputation(
    data            = assay_df,
    metadata_sample = sample_info_df,
    metadata_info   = metadata_info,
    core            = FALSE,
    mvi_percentage  = mvi_percentage
  ),
  error = function(e) abort("mvi_imputation() failed: ", conditionMessage(e))
)
if (is.null(imputed_df)) abort("mvi_imputation() did not return a data frame.")

# ── BUILD OUTPUT SE ──────────────────────────────────────────────────────────
# imputed_df is samples x features, the assay is features x samples in the order of se.
imputed_assay <- t(as.matrix(imputed_df))[rownames(se), colnames(se), drop = FALSE]
se_imputed <- se
assay(se_imputed, 1) <- imputed_assay

# A feature with no measured value in any sample of a condition has no minimum there and stays NA. It is removed, unless
# --remove_unimputable is FALSE.
na_features <- rownames(se_imputed)[rowSums(is.na(assay(se_imputed, 1))) > 0]
na_removed  <- if (remove_unimputable) na_features else character(0)
na_kept     <- setdiff(na_features, na_removed)
if (length(na_removed) > 0)
  log_msg("WARN", "Removed ", length(na_removed), " feature(s) with no measured value in any sample of at least one condition, ",
          "because there is no minimum to impute from: ", paste(na_removed, collapse = ", "), ".")
if (length(na_kept) > 0)
  log_msg("WARN", length(na_kept), " feature(s) stay NA in at least one condition, because --remove_unimputable is FALSE: ",
          paste(na_kept, collapse = ", "), ". They should be removed before processingnormalize, which expects data without ",
          "missing values.")
# One row per removed feature, with the reason.
removed_df <- data.frame(
  Feature = c(flagged_ids, na_removed),
  Reason  = c(rep("flagged", length(flagged_ids)), rep("unimputable", length(na_removed))),
  stringsAsFactors = FALSE
)
if (length(na_removed) >= nrow(se_imputed))
  abort("All ", nrow(se_imputed), " features have no measured value in any sample of at least one condition. No feature would ",
        "be left after the imputation.")

# The cells that were missing, with the value they got (NA if none could be imputed).
imputed_idx <- which(missing_mask, arr.ind = TRUE)
imputed_values_df <- data.frame(
  Metabolite = rownames(se)[imputed_idx[, "row"]],
  Sample     = colnames(se)[imputed_idx[, "col"]],
  Condition  = sample_info_df[colnames(se)[imputed_idx[, "col"]], conditions_col],
  `Imputed value` = imputed_assay[imputed_idx],
  check.names = FALSE
)
imputed_values_df <- imputed_values_df[order(imputed_values_df\$Metabolite,
                                              imputed_values_df\$Sample), , drop = FALSE]
# The removed features leave the SE, their cells stay in the table of imputed values.
se_imputed <- se_imputed[!rownames(se_imputed) %in% na_removed, , drop = FALSE]

n_imputed       <- sum(!is.na(imputed_values_df[["Imputed value"]]))
n_missing_after <- sum(is.na(assay(se_imputed, 1)))
log_msg("INFO", "Imputed ", n_imputed, " missing value(s), ", n_missing_after, " stay NA.")

# ── WRITE OUTPUTS ─────────────────────────────────────────────────────────────
log_section("Writing outputs")
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
if (output_format != "se" && length(assays(se_imputed)) > 1)
  log_msg("WARN", "The input has ", length(assays(se_imputed)), " assays. Only the first one is written to the data matrix.")
write_out(se_imputed, se_out);                      log_msg("INFO", "Written: ", paste(se_out, collapse = ", "))
write.table(imputed_values_df, values_out,
            sep = "\\t", quote = FALSE, row.names = FALSE)
log_msg("INFO", "Written: ", values_out)
write.table(removed_df, removed_out, sep = "\\t", quote = FALSE, row.names = FALSE)
log_msg("INFO", "Written: ", removed_out)

# ── HTML REPORT ────────────────────────────────────────────────────────────
log_section("Building report")

# HTML table of a data frame, in a box that scrolls when it is long.
build_table_html <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  header <- paste0("<th>", esc(colnames(df)), "</th>", collapse = "")
  is_num <- vapply(df, is.numeric, logical(1))
  rows <- apply(df, 1, function(r) {
    cells <- vapply(seq_along(r), function(i) {
      if (is_num[i]) {
        v <- suppressWarnings(as.numeric(r[[i]]))
        if (!is.na(v)) format(round(v, 4)) else as.character(r[[i]])
      } else as.character(r[[i]])
    }, character(1))
    paste0("<tr>", paste0("<td>", esc(cells), "</td>", collapse = ""), "</tr>")
  })
  paste0('<div class="scroll-table"><table><thead><tr>', header,
         '</tr></thead><tbody>', paste(rows, collapse = ""), '</tbody></table></div>')
}
imputed_values_tbl  <- build_table_html(imputed_values_df)
imputed_values_html <- if (!is.null(imputed_values_tbl)) imputed_values_tbl else
  "<p><em>None, no missing values were found.</em></p>"

imputation_method_label <- if (isTRUE(all.equal(mvi_percentage, 50))) {
  "50% of the minimum (half-minimum), per condition"
} else {
  sprintf("%s%% of the minimum, per condition", format(mvi_percentage))
}

# Names for the report: at most 25, and a pointer to the log file, which has all of them, only if some are left out.
names_html <- function(x, max_n = 25) paste0(
  paste0("<code>", esc(head(x, max_n)), "</code>", collapse = " "),
  if (length(x) > max_n) paste0(" and ", length(x) - max_n, " more (see <code>", esc(log_out), "</code>)") else "")
# Feature names for the report: all of them, one per line, in a list that scrolls when it is long.
feat_names_html <- function(x) paste0('<div class="namelist">', paste0("<code>", esc(x), "</code>", collapse = "<br>"), "</div>")
graybox <- function(...) paste0('<div class="graybox">', paste0(...), "</div>")
sel <- function(param, value, default) {
  paste0("<code>", param, " ", esc(value), "</code> ",
         if (identical(as.character(value), as.character(default))) "(default)"
         else paste0("(default: ", esc(default), ")"))
}
pl <- function(n, one, many) if (n == 1) one else many
# Gray box under the imputed values table: the features that stay NA because --remove_unimputable is FALSE.
na_box_html <- if (length(na_kept) > 0) graybox(
  if (length(na_kept) > 0) paste0(
    "<p><strong>", length(na_kept), " feature", pl(length(na_kept), " stays", "s stay"), " NA in at least one condition:</strong></p><p>",
    feat_names_html(na_kept), "</p><p>", pl(length(na_kept), "It was", "They were"), " not removed, because <code>--remove_unimputable FALSE</code> is set, ",
    "but should be removed before continuing with <code>processingnormalize</code>, which expects data without missing values.</p>") else "") else ""
# Gray box under the conditions column: samples without a condition, or only one group of samples.
group_values <- unique(sample_info_df[[conditions_col]])
cond_box_html <- if (length(no_condition_ids) > 0 || length(group_values) == 1) graybox(
  if (length(no_condition_ids) == ncol(se)) paste0(
    "<p><strong>No sample has a value in the Conditions column.</strong> All samples are treated as ",
    "one condition called 'NA', and the minimum of each feature is taken over all samples.</p>") else paste0(
    if (length(no_condition_ids) > 0) paste0(
      "<p><strong>", length(no_condition_ids), " sample", pl(length(no_condition_ids), " has", "s have"), " no value in the Conditions column:</strong></p><p>", names_html(no_condition_ids), "</p><p>",
      pl(length(no_condition_ids), "This sample is", "These samples are"), " treated as one condition called 'NA'.</p>") else "",
    if (length(group_values) == 1) paste0(
      "<p><strong>All samples are in one condition (<code>", esc(group_values), "</code>).</strong> The minimum of each feature is taken over all samples.</p>") else "")) else ""

params_html <- paste0(
  "<ul>",
  "<li><p>", sel("--mvi_percentage", format(mvi_percentage), "50"),
  ": the imputed value as a percentage (%) of the feature minimum, between 0 and 100.</p></li>",
  "<li><p>", sel("--conditions_col", conditions_col, "Conditions"),
  ": the column of the ", samp_meta_txt, " that holds the conditions. The minimum is taken per condition.",
  "</p>", cond_box_html, "</li>",
  if (filter_flagged) paste0("<li><p>", sel("--filteredfeatures_col", filteredfeatures_col, "FilteredFeaturesFlag"),
                             ": the column of the ", feat_meta_txt, " that marks the features to remove with TRUE.</p>",
                             if (!flag_col_found) graybox("<p><strong>There is no such column in the feature metadata, and therefore no feature was removed.</strong></p>") else
                               if (length(flagged_ids) == 0) graybox("<p><strong>No feature is flagged TRUE, and all features are kept.</strong></p>") else
                                 graybox("<p><strong>", length(flagged_ids), " flagged feature", pl(length(flagged_ids), " was", "s were"), " removed:</strong></p><p>",
                                         feat_names_html(flagged_ids), "</p>"),
                             "</li>") else "",
  "<li><p>", sel("--filter_flagged", if (filter_flagged) "TRUE" else "FALSE", "TRUE"), ": ",
  if (filter_flagged) "the flagged features are removed before the imputation." else
    "the flagged features are not removed and are imputed like all others.", "</p></li>",
  "<li><p>", sel("--remove_unimputable", if (remove_unimputable) "TRUE" else "FALSE", "TRUE"), ": ",
  if (remove_unimputable) "the features with no measured value in at least one condition are removed." else
    "the features with no measured value in at least one condition are kept with NA values.", "</p>",
  if (length(na_removed) > 0) graybox(
    "<p><strong>", length(na_removed), " unimputable feature", pl(length(na_removed), " was", "s were"), " removed:</strong></p><p>",
    feat_names_html(na_removed), "</p>") else "",
  "</li>",
  if (length(keep_asked) > 0) paste0(
    "<li><p>", sel("--keep_features", paste(keep_asked, collapse = ","), "none"),
    ": features that are kept even if they are flagged TRUE.</p>",
    if (length(kept_flagged_ids) > 0) graybox(
      "<p><strong>", length(kept_flagged_ids), " flagged feature", pl(length(kept_flagged_ids), " was", "s were"), " kept:</strong></p><p>",
      feat_names_html(kept_flagged_ids), "</p>") else "",
    if (length(keep_not_found) > 0) graybox(
      "<p><strong>The list contains ", pl(length(keep_not_found), "a feature that is", "features that are"), " not in the input:</strong></p><p>",
      feat_names_html(keep_not_found), "</p>") else "",
    "</li>") else "",
  "</ul>")

# The feature count of the output, with the features that were removed before it.
n_fl <- length(flagged_ids); n_na <- length(na_removed)
features_txt <- paste0(nrow(se_imputed), " feature", pl(nrow(se_imputed), "", "s"),
  if (n_fl + n_na == 0) "" else if (n_na == 0) paste0(" after ", n_fl, " flagged feature", pl(n_fl, " was", "s were"), " removed") else
    if (n_fl == 0) paste0(" after ", n_na, " feature", pl(n_na, " with", "s with"), " no measured value in a condition ",
                          pl(n_na, "was", "were"), " removed") else
      paste0(" after ", n_fl + n_na, " features were removed (", n_fl, " flagged, ", n_na, " with no measured value in a condition)"))
run_meta_html <- sprintf(
'<dl>
   <dt>Input</dt><dd>%s</dd>
   <dt>n features</dt><dd>%s</dd>
   <dt>n samples</dt><dd>%d</dd>
   <dt>Imputation method</dt><dd>%s</dd>
   <dt>Imputed values</dt><dd>%d</dd>
 </dl>',
 esc(source_label),
 features_txt,
 ncol(se), imputation_method_label,
 n_imputed)

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
<title>Processing: Imputation Report</title>
<style>
  body{font-family:Arial,sans-serif;max-width:980px;margin:2em auto;padding:0 1em;line-height:1.5}
  h1{color:#1F3A5F;border-bottom:2px solid #1F3A5F;padding-bottom:.3em}
  h2{color:#1F3A5F;margin-top:2em}
  dl{display:grid;grid-template-columns:max-content 1fr;gap:.3em 1em}
  dt{font-weight:bold;color:#345}
  code{background:#f2f2f2;padding:.15em .35em;border-radius:3px}
  table{border-collapse:collapse;width:100%%;font-size:.9em}
  th,td{border:1px solid #ccc;padding:.35em .6em;text-align:left}
  th{background:#f0f4f8}
  .graybox{background:#f2f2f2;border-left:4px solid #999;padding:.1em 1em;margin:1em 0;border-radius:3px}
  .graybox code{background:#e0e0e0}
  .namelist{max-height:9.5em;overflow-y:auto}
  .scroll-table{max-height:320px;overflow-y:auto;border:1px solid #cbd6e4;border-radius:6px;margin:.4em 0}
  .scroll-table table{border:0;margin:0}
  .scroll-table thead th{position:sticky;top:0;background:#f0f4f8;z-index:1}
</style></head><body>
<h1>Processing: Missing-Value Imputation</h1>
<p><em>Generated by the METAPROVIZ_PROCESSINGIMPUTE module.</em></p>

<h2>Run summary</h2>%s

<h2>Citations</h2>
<p>
If you use results from this analysis in a publication, please cite:
</p>
<ul>
  <li><strong>MetaProViz:</strong>
  <a href="https://doi.org/10.1038/s44320-026-00231-8">Schmidt et al., Integrated metabolomics data analysis to generate mechanistic hypotheses with MetaProViz, Molecular Systems Biology 2026.</a></li>
  <li><strong>Dependencies:</strong>
  <a href="https://doi.org/10.32614/R.manuals">R Core Team, R: A Language and Environment for Statistical Computing, R Foundation for Statistical Computing</a>, <a href="https://doi.org/10.18129/B9.bioc.SummarizedExperiment">Morgan et al., SummarizedExperiment: A container (S4 class) for matrix-like assays, Bioconductor</a> and <a href="https://doi.org/10.1038/nmeth.3252">Huber et al., Orchestrating high-throughput genomic analysis with Bioconductor, Nature Methods 2015</a>.
  The other R packages used are listed with their versions in the session summary at the end of this report.</li>
  <li><strong>Imputation method:</strong>
  <a href="https://doi.org/10.1038/s41598-017-19120-0">Wei et al., Missing Value Imputation Approach for Mass Spectrometry-based Metabolomics Data, Scientific Reports 2018.</a></li>
  <li><strong>This module:</strong> <code>metaproviz/processingimpute</code>, an nf-core module. It has no publication of its own and is cited through the nf-core publication below. In the methods, name the module.</li>
  <li><strong>nf-core:</strong> <a href="https://doi.org/10.1038/s41587-020-0439-x">Ewels et al., The nf-core framework for community-curated bioinformatics pipelines, Nature Biotechnology 2020.</a></li>
  <li><strong>Nextflow:</strong> <a href="https://doi.org/10.1038/nbt.3820">Di Tommaso et al., Nextflow enables reproducible computational workflows, Nature Biotechnology 2017.</a></li>
</ul>

<h2>Description</h2>
<p>
This is the second module in the processing chain (<code>processingfilter</code>, <code>processingimpute</code>,
<code>processingnormalize</code> and <code>processingoutliers</code>). It continues from the output of
<code>processingfilter</code>. The data should therefore contain no pool or QC samples, and by default the features
flagged in the data are removed before the imputation.
</p>
<p>
This module replaces each missing value (NA or 0) with a percentage of the minimum observed value of that feature,
taken separately within each condition (<code>--mvi_percentage</code>, default 50%%). At 50%% this is the half-minimum
approach (Wei et al., 2018), and any other percentage is a generalized fraction of the minimum. A feature with no
measured value in any sample of a condition has no minimum there and cannot be imputed. By default, these features are removed, because the normalization in the next module, <code>processingnormalize</code>,
needs data without missing values. To keep them with NA values, use <code>--remove_unimputable FALSE</code>.
</p>
<p><strong>The parameters selected for this run:</strong></p>
%s

<h2>Imputed values</h2>
<p>
Every imputed value with its feature (<code>Metabolite</code>), its <code>Sample</code>, the <code>Condition</code> of that
sample and the value that replaced the missing one (<code>Imputed value</code>). Saved as <code>%s</code>. The value is NA where the feature has no measured value
in the condition and cannot be imputed.
</p>
%s
%s

<h2>Run notes (log)</h2>
<p>A separate file, <code>%s</code>, records every step of this run,
including any warnings. Check it first if a result here looks unexpected.</p>

<h2>Session summary</h2>
<p>Full R session used to generate this report — every package, and its version:</p>
%s

</body></html>',
run_meta_html, params_html, esc(values_out), imputed_values_html, na_box_html, esc(log_out), session_summary_html)

writeLines(html, report_out)
log_msg("INFO", "Written: ", report_out)

# MetaProViz and OmnipathR write their own log files and plot dumps to the working directory. They are no outputs of
# this module and are deleted here, at the very end.
unlink("metaproviz-log", recursive = TRUE)
unlink("omnipathr-log", recursive = TRUE)
unlink(Sys.glob("Rplots*.pdf"))

flush_log()
message(sprintf(
  "\\n✓ processing_impute complete — %d missing value(s) imputed | %d warning(s) | %d error(s)",
  n_imputed, warn_count, err_count))
