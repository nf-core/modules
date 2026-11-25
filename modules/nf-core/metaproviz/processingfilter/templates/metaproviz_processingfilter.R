#!/usr/bin/env Rscript
# =============================================================================
# metaproviz_processingfilter.R: feature filtering step of the MetaProViz processing chain
#
# Approach: MetaProViz::feature_filtering() flags the features to remove. No feature is removed from the data, the flags
# are added to the feature metadata as one TRUE/FALSE column (TRUE = remove suggestion).
#
# Filtering options (--featurefilt):
#   Modified  detection rate per condition (Yang et al., 2015), the default
#   Standard  detection rate across all samples (Bijlsma et al., 2006)
#   Manual    features from a list (--features_to_flag) or a TRUE/FALSE column (--filteredfeatures_col), no cutoff
#
# If something an option needs is missing, the module warns and moves one step down instead of failing:
# Manual (needs the list or the column) -> Modified (needs a Conditions column with at least two conditions) -> Standard.
#
# Pool/QC samples must be removed before this module (for example with the metaproviz/poolestimation module), they are
# not detected here and would be used like all other samples.
#
# Outputs, all named <prefix>.<suffix>: the SE (or three tables if output_format is tsv or csv), flagged_features.tsv,
# the HTML report and the log.
# =============================================================================

set.seed(42)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(magrittr)
  library(SummarizedExperiment); library(S4Vectors)
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
  conditions_col       = "Conditions",
  filteredfeatures_col = "FilteredFeaturesFlag",
  biorep_col         = "Biological_Replicates",
  featurefilt        = "Modified",
  cutoff_featurefilt = "0.8",
  features_to_flag   = character(0)
)
# R raw string: quotes and apostrophes inside the arguments (feature names like 5'-AMP) cannot break the R code.
args_opt <- parse_args(r"---($task.ext.args)---", list_opts = c("features_to_flag"))
for (ao in names(args_opt)) {
  if (!ao %in% names(opt)) stop(paste("Invalid option:", ao))
  opt[[ao]] <- args_opt[[ao]]
}

se_path             <- nullify('$se_rds')
data_matrix_path    <- nullify('$data_matrix')
feature_matrix_path <- nullify('$feature_matrix')
sample_matrix_path  <- nullify('$sample_matrix')
conditions_col       <- opt\$conditions_col
filteredfeatures_col <- opt\$filteredfeatures_col
biorep_col         <- opt\$biorep_col
featurefilt        <- opt\$featurefilt
cutoff_featurefilt <- as.numeric(opt\$cutoff_featurefilt)
features_to_flag_list <- opt\$features_to_flag
features_to_flag_str  <- paste(features_to_flag_list, collapse = ",")

prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
# output_format is a val input of the module. se: the SE as one RDS. tsv/csv: three files (data, feature and sample
# matrix), the same shape as the input matrices. An empty value means se.
output_format_in   <- trimws('$output_format')
output_format      <- if (output_format_in %in% c("", "null", "[]")) "se" else tolower(output_format_in)
se_out             <- if (output_format == "se") paste0(prefix, ".flagged.rds") else
  paste0(prefix, ".flagged", c(".data_matrix.", ".feature_matrix.", ".sample_matrix."), output_format)
flagged_features_out <- paste0(prefix, ".flagged_features.tsv")
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
abort <- function(...) { log_msg("ERROR", ...); flush_log(); stop(paste0(...), call. = FALSE) }
esc <- function(x) gsub(">", "&gt;", gsub("<", "&lt;",
                    gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE), fixed = TRUE)

if (!output_format %in% c("se", "tsv", "csv"))
  abort("output_format must be se, tsv or csv (empty means se), not '", output_format_in, "'.")

# ── LOAD INPUT ──────────────────────────────────────────────────────────────
log_section("Loading input data")

# Two acceptable input shapes, mutually exclusive: --se (a SummarizedExperiment
# .rds), or all three of --data_matrix/--feature_matrix/--sample_matrix (plain
# TSV/CSV, language-independent).
have_se  <- !is.null(se_path) && nzchar(se_path)
# how the report names the feature and sample metadata, depending on how the input was given
feat_meta_txt <- if (have_se) {
  sprintf("feature metadata (<code>rowData</code> of <code>%s</code>)", esc(basename(se_path)))
} else sprintf("feature metadata (<code>%s</code>)", esc(basename(feature_matrix_path)))
# the feature metadata of the output: rowData of the SE file, or the feature_matrix table
out_feat_meta_txt <- if (output_format == "se") {
  sprintf("feature metadata (<code>rowData</code> of <code>%s</code>)", esc(se_out))
} else sprintf("feature metadata (<code>%s</code>)", esc(se_out[2]))
samp_meta_txt <- if (have_se) {
  sprintf("sample metadata (<code>colData</code> of <code>%s</code>)", esc(basename(se_path)))
} else sprintf("sample metadata (<code>%s</code>)", esc(basename(sample_matrix_path)))
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

n_features_input <- nrow(se); n_samples_input <- ncol(se)

# ── FEATURE FILTERING ────────────────────────────────────────────────────────
log_section("Running MetaProViz::feature_filtering()")

assay_df        <- t(assay(se, 1)) %>% as.data.frame(check.names = FALSE)
sample_info_df  <- as.data.frame(colData(se), check.names = FALSE)
feature_info_df <- as.data.frame(rowData(se), check.names = FALSE)

if (!featurefilt %in% c("Standard", "Modified", "Manual"))
  abort("--featurefilt must be 'Standard', 'Modified' or 'Manual', got '", featurefilt, "'.")

# Fallback chain Manual -> Modified -> Standard, a warning each time and no failure.
featurefilt_requested <- featurefilt

# Manual with a list of names: the TRUE/FALSE column is built from the list (TRUE = flag). The list takes precedence
# over a column of the same name in the feature metadata, which is kept under another name (see BUILD OUTPUT SE).
use_list       <- FALSE
list_asked     <- character(0)
list_not_found <- character(0)
if (nzchar(features_to_flag_str)) {
  if (featurefilt != "Manual") {
    log_msg("WARN", "--features_to_flag is only used with --featurefilt Manual and is ignored.")
  } else {
    use_list       <- TRUE
    list_asked     <- unique(trimws(features_to_flag_list)); list_asked <- list_asked[nzchar(list_asked)]
    list_not_found <- setdiff(list_asked, rownames(se))
    log_msg("INFO", "Features to flag (", length(list_asked), "): ", paste(list_asked, collapse = ", "))
    if (length(list_not_found) > 0)
      log_msg("WARN", length(list_not_found), " feature(s) of --features_to_flag are not in the data: ",
              paste(list_not_found, collapse = ", "), ".")
    feature_info_df[[filteredfeatures_col]] <- rownames(se) %in% list_asked
  }
}

fallback_notes        <- character(0)
fallback_sent         <- character(0)   # the same, as report sentences that say what is missing

if (featurefilt == "Manual" && !filteredfeatures_col %in% colnames(feature_info_df)) {
  reason <- sprintf("there is no column '%s' in rowData", filteredfeatures_col)
  log_msg("WARN", "--featurefilt Manual takes the features to remove from a TRUE/FALSE column, but ",
          reason, ". Using Modified instead.")
  fallback_notes <- c(fallback_notes, sprintf("Manual to Modified, because %s", reason))
  fallback_sent  <- c(fallback_sent, sprintf(paste0(
    "<code>--featurefilt Manual</code> cannot be used because there is no column <code>%s</code> in the %s."),
    esc(filteredfeatures_col), feat_meta_txt))
  featurefilt <- "Modified"
}

# Samples with NA, an empty value or only whitespace in the Conditions column are all grouped together as one condition
# called "NA", and a warning lists them. Numbers in the column are turned into text, so they are conditions like any other.
has_conditions   <- conditions_col %in% colnames(sample_info_df)
no_condition_ids <- character(0)
if (has_conditions) {
  cond_chr <- as.character(sample_info_df[[conditions_col]])
  no_cond  <- is.na(cond_chr) | trimws(cond_chr) == ""
  if (any(no_cond)) {
    no_condition_ids <- rownames(sample_info_df)[no_cond]
    cond_chr[no_cond] <- "NA"
    log_msg("WARN", length(no_condition_ids), " sample(s) have no value in the Conditions column '", conditions_col,
            "': ", paste(no_condition_ids, collapse = ", "), ". They are treated as the condition 'NA'.")
  }
  sample_info_df[[conditions_col]] <- cond_chr
}
n_conditions   <- if (has_conditions) length(unique(sample_info_df[[conditions_col]])) else 0L
if (featurefilt == "Modified" && n_conditions < 2) {
  reason <- if (!has_conditions) {
    sprintf("there is no column '%s' in colData", conditions_col)
  } else {
    sprintf("the column '%s' has fewer than two different values", conditions_col)
  }
  log_msg("WARN", "--featurefilt Modified checks every condition separately, but ",
          reason, ". Using Standard (one check across all samples) instead.")
  fallback_notes <- c(fallback_notes, sprintf("Modified to Standard, because %s", reason))
  what <- if (!has_conditions) {
    paste0("there is no column <code>", esc(conditions_col), "</code> in the ", samp_meta_txt)
  } else {
    paste0("the column <code>", esc(conditions_col), "</code> has fewer than two different values")
  }
  fallback_sent  <- c(fallback_sent, paste0("<code>--featurefilt Modified</code> cannot be used because ", what, "."))
  featurefilt <- "Standard"
}

# Manual: read the TRUE/FALSE column (TRUE = remove). TRUE/FALSE, T/F, 1/0 and yes/no are accepted, empty means keep.
# Anything else stops the run.
metadata_feature <- NULL
n_flag_true <- n_flag_false <- n_flag_empty <- 0L
if (featurefilt == "Manual") {
  raw_flags <- feature_info_df[[filteredfeatures_col]]
  chr       <- toupper(trimws(as.character(raw_flags)))
  flags_in  <- if (is.logical(raw_flags)) raw_flags else {
    out <- rep(NA, length(chr))
    out[chr %in% c("TRUE", "T", "1", "YES")]  <- TRUE
    out[chr %in% c("FALSE", "F", "0", "NO")]  <- FALSE
    out
  }
  blank <- is.na(raw_flags) | chr == ""
  bad   <- !blank & is.na(flags_in)
  if (any(bad))
    abort("The column '", filteredfeatures_col, "' must hold TRUE or FALSE (T/F, 1/0 or yes/no also work), ",
          "but ", sum(bad), " feature(s) have other values, for example: ",
          paste(head(unique(as.character(raw_flags)[bad]), 5), collapse = ", "), ".")
  n_flag_true  <- sum(flags_in %in% TRUE)
  n_flag_false <- sum(flags_in %in% FALSE)
  n_flag_empty <- sum(blank)
  if (n_flag_empty > 0)
    log_msg("WARN", n_flag_empty, " feature(s) have no value in '", filteredfeatures_col,
            "' and are kept.")
  if (n_flag_true == 0)
    log_msg("WARN", if (use_list) "None of the features of --features_to_flag is in the data, so no feature will be flagged." else
            paste0("No feature is TRUE in '", filteredfeatures_col, "', so no feature will be flagged."))
  metadata_feature <- data.frame(flag = flags_in, row.names = rownames(se))
  names(metadata_feature) <- filteredfeatures_col
}

metadata_info <- if (has_conditions) c(Conditions = conditions_col) else NULL
if (biorep_col %in% colnames(sample_info_df))
  metadata_info <- c(metadata_info, Biological_Replicates = biorep_col)
if (featurefilt == "Manual")
  metadata_info <- c(metadata_info, FilteredFeatures = filteredfeatures_col)

filt_res <- tryCatch(
  MetaProViz::feature_filtering(
    data                = assay_df,
    metadata_sample     = sample_info_df,
    metadata_info       = metadata_info,
    core                = FALSE,
    featurefilt         = featurefilt,
    cutoff_featurefilt  = cutoff_featurefilt,
    metadata_feature    = metadata_feature
  ),
  error = function(e) abort("feature_filtering() failed: ", conditionMessage(e))
)

filtered_df <- filt_res[["DF"]]
flags       <- filt_res[["FilteredFeatures"]]
if (is.null(filtered_df)) abort("feature_filtering() did not return DF.")
if (is.null(flags)) abort("feature_filtering() did not return FilteredFeatures.")
removed <- flags[["FeatureID"]][flags[["FilteredFeatures"]]]
# From here on, `featurefilt` is the method MetaProViz applied.
if (!is.null(filt_res[["featurefilt_used"]])) featurefilt <- filt_res[["featurefilt_used"]]

if (featurefilt == "Manual") {
  log_msg("INFO", "Flagged ", length(removed), " / ", n_features_input,
          " feature(s) as TRUE in '", filteredfeatures_col, "'.")
} else {
  log_msg("INFO", "Flagged ", length(removed), " / ", n_features_input,
          " feature(s) below the detection-rate cutoff.")
}

# % missing (NA or 0, as feature_filtering() counts it) of each flagged feature: one column per condition for Modified,
# one overall column for Standard and Manual.
flagged_features_df <- if (length(removed) == 0) {
  data.frame(Metabolite = character(0), check.names = FALSE)
} else if (featurefilt == "Modified") {
  long_df <- assay_df[, removed, drop = FALSE] %>%
    tibble::rownames_to_column("Sample") %>%
    tidyr::pivot_longer(-Sample, names_to = "Metabolite", values_to = "Value") %>%
    dplyr::mutate(Condition = sample_info_df[Sample, conditions_col])
  df <- long_df %>%
    dplyr::group_by(Metabolite, Condition) %>%
    dplyr::summarise(`% missing` = round(100 * mean(is.na(Value) | Value == 0), 1),
                      .groups = "drop") %>%
    tidyr::pivot_wider(id_cols = Metabolite, names_from = Condition, values_from = `% missing`)
  df[match(removed, df\$Metabolite), , drop = FALSE]
} else {
  pct_missing <- vapply(removed, function(feat) {
    v <- assay_df[[feat]]
    100 * mean(is.na(v) | v == 0)
  }, numeric(1))
  df <- data.frame(Metabolite = removed, `% missing` = round(pct_missing, 1),
                    check.names = FALSE)
  df[order(-df\$`% missing`), , drop = FALSE]
}

# ── BUILD OUTPUT SE ───────────────────────────────────────────────────────
# The SE is not filtered. All features stay, with their values unchanged, and the flags (TRUE = remove) are added
# as one column of the feature metadata, named like --filteredfeatures_col.
kept_features <- colnames(filtered_df)
stopifnot(setequal(c(kept_features, removed), rownames(se)))
# If the input already has a column with this name and the flags of this run are not that column (a list was given, or
# Modified or Standard was used), the input column is kept as <name>_original. With Manual and only a column, the
# column is the input and is written back under its name.
old_column_kept_as <- NULL
if (filteredfeatures_col %in% colnames(rowData(se)) && (use_list || featurefilt != "Manual")) {
  old_column_kept_as <- paste0(filteredfeatures_col, "_original")
  n_taken <- 2L
  while (old_column_kept_as %in% colnames(rowData(se))) {
    old_column_kept_as <- paste0(filteredfeatures_col, "_original", n_taken); n_taken <- n_taken + 1L
  }
  log_msg("WARN", "The feature metadata already has a column '", filteredfeatures_col, "'. It is kept as '",
          old_column_kept_as, "', and the flags of this run are written to '", filteredfeatures_col, "'.")
}
se_flagged <- se
if (!is.null(old_column_kept_as)) {
  rd_in <- rowData(se_flagged)
  colnames(rd_in)[colnames(rd_in) == filteredfeatures_col] <- old_column_kept_as
  rowData(se_flagged) <- rd_in
}
rowData(se_flagged)[[filteredfeatures_col]] <- flags[["FilteredFeatures"]][match(rownames(se), flags[["FeatureID"]])]

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
if (output_format != "se") {
  if (length(assays(se_flagged)) > 1)
    log_msg("WARN", "The input has ", length(assays(se_flagged)), " assays. Only the first one is written to the data matrix.")
}
write_out(se_flagged, se_out);                      log_msg("INFO", "Written: ", paste(se_out, collapse = ", "))
write.table(flagged_features_df, flagged_features_out,
            sep = "\\t", quote = FALSE, row.names = FALSE)
log_msg("INFO", "Written: ", flagged_features_out)

# ── HTML REPORT ────────────────────────────────────────────────────────────
log_section("Building report")

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
flagged_features_tbl  <- build_table_html(flagged_features_df)
flagged_features_html <- if (!is.null(flagged_features_tbl)) flagged_features_tbl else
  if (featurefilt == "Manual") "<p><em>None, no feature is flagged.</em></p>" else
    "<p><em>None, every feature passed the detection-rate cutoff.</em></p>"

allowed_pct <- format(round(100 * (1 - cutoff_featurefilt), 1))
missing_def <- paste0("Zeros and NAs in the data both count as missing values, meaning the feature was not detected in ",
                      "that sample.")
flagged_features_lead <- if (featurefilt == "Modified") {
  paste0("The table shows the share of samples (%) in which each flagged feature is missing, one column ",
         "per condition.")
} else if (featurefilt == "Manual") {
  "The table shows the share of samples (%) in which each flagged feature is missing."
} else {
  paste0("The table shows the share of samples (%) in which each flagged feature is missing, across all ",
         "samples.")
}
flagged_features_note <- if (featurefilt == "Manual") {
  paste0(
    "These features are flagged for removal only because ", if (use_list) "they are in <code>--features_to_flag</code>" else
      paste0("they are TRUE in the feature column <code>", esc(filteredfeatures_col), "</code>"), ". No cutoff is used.")
} else if (featurefilt == "Modified") {
  paste0(
    missing_def, " A feature is flagged for the whole dataset as soon as its share of missing values is above the ",
    "missingness allowed by <code>--cutoff_featurefilt = ", esc(format(cutoff_featurefilt)), "</code> (at most ", allowed_pct, "% missing) in even ",
    "one condition, not necessarily in all of them.")
} else {
  paste0(
    missing_def, " A feature is flagged when its share of missing values is above the missingness allowed by ",
    "<code>--cutoff_featurefilt = ", esc(format(cutoff_featurefilt)), "</code> (at most ", allowed_pct, "% missing).")
}

if (!is.null(old_column_kept_as))
  flagged_features_note <- paste0(flagged_features_note, " The input already had a column <code>", esc(filteredfeatures_col),
    "</code>. It is kept as <code>", esc(old_column_kept_as), "</code>, and the flags of this run are in <code>",
    esc(filteredfeatures_col), "</code>.")

fallback_from <- sub(" to .*", "", fallback_notes)
filter_method_label <- featurefilt
fallback_box_html <- if (length(fallback_notes) == 0) "" else paste0(
  "<div class=\"graybox\"><p><strong>Fallback:</strong> ",
  paste0("<code>--featurefilt ", fallback_from, "</code>", collapse = " and "), " could not be used, therefore <code>--featurefilt ",
  featurefilt, "</code> is used as fallback.</p></div>")
# the cutoff row is only shown for the options that use a cutoff
cutoff_row <- if (featurefilt == "Manual") "" else
  sprintf("\n   <dt>Detection-rate cutoff</dt><dd>%s</dd>", format(cutoff_featurefilt))

run_meta_html <- sprintf(
'<dl>
   <dt>Input</dt><dd>%s</dd>
   <dt>n features</dt><dd>%s</dd>
   <dt>n samples</dt><dd>%s</dd>
   <dt>Filter method</dt><dd>%s</dd>%s
 </dl>',
 esc(source_label),
 sprintf("%d feature%s, %d flagged for removal", n_features_input, if (n_features_input == 1) "" else "s", length(removed)),
 paste0(ncol(se), " sample", if (ncol(se) == 1) "" else "s"),
 esc(filter_method_label), cutoff_row)
run_meta_html <- paste0(run_meta_html, fallback_box_html)

# ── DESCRIPTION (explains only the filtering option that was used) ──────────
# tolerated(n): the largest number of missing values out of n samples that still keeps a feature, with the same
# comparison as MetaProViz, so the examples below match what the module does.
tolerated_at <- function(n, cutoff) max(which(!((0:n) > (1 - cutoff) * n)) - 1L)
tolerated    <- function(n) tolerated_at(n, cutoff_featurefilt)
pl        <- function(k) if (k == 1) "value" else "values"
max_missing_pct <- format(round(100 * (1 - cutoff_featurefilt), 1))
cutoff_txt      <- format(cutoff_featurefilt)

# a parameter as the user set it: "--name value" followed by "(default)" or, if changed, "(default: X)"
sel <- function(param, value, default) {
  paste0("<code>", param, " ", esc(value), "</code> ",
         if (identical(as.character(value), as.character(default))) "(default)"
         else paste0("(default: ", esc(default), ")"))
}
cutoff_is_default <- isTRUE(all.equal(cutoff_featurefilt, 0.8))

# how a high cutoff hits small groups (a condition for Modified, the dataset for Standard), always with the default 0.8
small_unit <- if (featurefilt == "Modified") "a condition" else "a dataset"
example5 <- paste0(
  "Look out for small sample sizes. For example, with the default 0.8, ", small_unit, " with 5 samples tolerates ",
  if (tolerated_at(5, 0.8) == 0) "no missing value, and a single one flags the feature." else
    sprintf("up to %d missing %s, one more flags the feature.", tolerated_at(5, 0.8), pl(tolerated_at(5, 0.8))))

in_this_run <- if (featurefilt == "Modified") {
  sizes <- table(sample_info_df[[conditions_col]]); lo <- min(sizes); hi <- max(sizes)
  if (lo == hi) {
    sprintf("In this run each condition has %d samples: up to %d missing %s per condition are tolerated, one more flags the feature for the whole dataset.",
            lo, tolerated(lo), pl(tolerated(lo)))
  } else {
    tl <- tolerated(lo); th <- tolerated(hi)
    tol_txt <- if (tl == th) sprintf("up to %d missing %s", tl, pl(tl)) else if (th - tl == 1)
      sprintf("up to %d or %d missing values", tl, th) else sprintf("between %d and %d missing values", tl, th)
    sprintf("In this run the conditions have %d to %d samples: %s per condition are tolerated, one more flags the feature for the whole dataset.",
            lo, hi, tol_txt)
  }
} else if (featurefilt == "Standard") {
  sprintf("In this run %d samples are used: up to %d missing %s are tolerated, one more flags the feature.",
          ncol(se), tolerated(ncol(se)), pl(tolerated(ncol(se))))
} else ""

intro_html <- paste0(
  "<p>This module is the first of the processing chain (<code>processingfilter</code>, <code>processingimpute</code>, ",
  "<code>processingnormalize</code> and <code>processingoutliers</code>). It flags metabolite features for removal, ",
  "either those that were detected in too few samples or those on a custom list. The flags are a TRUE/FALSE column ",
  "in the feature metadata of the output, and no feature is removed.</p>",
  "<p>Pool and QC samples must be removed before this module (for example with <code>metaproviz/poolestimation</code>), ",
  "since they are technical quality-control mixtures and not biological samples. A pool contains the features of all ",
  "study samples, which makes the features look better detected than they are.</p>",
  "<p><strong>The parameters selected for this run:</strong></p>")

# --featurefilt: what the option does
mode_text <- switch(featurefilt,
  Modified = paste0(
    "checks the detection rate separately within each condition (Yang et al., 2015). The flag is not per ",
    "condition: if a feature fails the check in even one condition, it is flagged for the entire dataset. It ",
    "needs at least two conditions in <code>--conditions_col</code>."),
  Standard = paste0(
    "checks the detection rate once, across all samples and independent of the conditions ",
    "(Bijlsma et al., 2006)."),
  Manual = if (use_list) paste0(
    "flags exactly the features given in <code>--features_to_flag</code>.") else paste0(
    "flags exactly the features that are marked with TRUE in <code>--filteredfeatures_col</code>."))

# After a fallback the item describes the option that was used, and a gray box says what was missing for the selected one.
fallback_item_box <- if (length(fallback_notes) == 0) "" else paste0(
  "<div class=\"graybox\"><p><strong>Fallback from ", paste(fallback_from, collapse = ", then "), ":</strong> ",
  paste(fallback_sent, collapse = " "), "</p></div>")
item_featurefilt <- paste0(
  "<li><p>", if (length(fallback_notes) == 0) sel("--featurefilt", featurefilt, "Modified") else
    paste0("<code>--featurefilt ", featurefilt, "</code>"), ": ",
  "<strong>", featurefilt, "</strong> ", mode_text, "</p>", fallback_item_box, "</li>")

item_cutoff <- if (featurefilt %in% c("Modified", "Standard")) {
  where <- if (featurefilt == "Modified") "of the samples of a condition" else "of all samples"
  price_text <- if (featurefilt == "Modified") paste0(
    "Any cutoff above 0 has a price. A metabolite that is really absent in one condition, for example because a ",
    "pathway is switched off there, fails the check as well") else paste0(
    "Any cutoff above 0 has a price. A metabolite that is really absent in one condition, for example because a pathway ",
    "is switched off there, is flagged as well if that condition makes up more samples than the cutoff allows to ",
    "be missing")
  paste0(
    "<li><p>", sel("--cutoff_featurefilt", cutoff_txt, "0.8"), ": the minimum required detection rate, a number ",
    "between 0 and 1 (", cutoff_txt, " = ", format(round(100 * cutoff_featurefilt, 1)), "%). A feature is flagged ",
    "once it is missing in more than the allowed share of samples, here ", max_missing_pct, "% ", where, ". ",
    "<strong>", in_this_run, "</strong></p>",
    "<p>Features with too many missing values should be removed, because missing values make a comparison between ",
    "groups less reliable. What counts as too many is a choice. ", example5, "</p>",
    "<p>", price_text, " (see &quot;Noise or biology?&quot; below).</p></li>")
} else {
  ""
}

# The conditions column is only listed if Modified was used.
# Names for the report: all of them up to 25, beyond that a pointer to the log file, which has all of them.
names_html <- function(x, max_n = 25) paste0(
  paste0("<code>", esc(head(x, max_n)), "</code>", collapse = " "),
  if (length(x) > max_n) paste0(" and ", length(x) - max_n, " more (see <code>", esc(log_out), "</code>)") else "")
item_conditions_col <- if (featurefilt == "Modified") {
  paste0("<li><p>", sel("--conditions_col", conditions_col, "Conditions"),
         ": the column of the ", samp_meta_txt, " that holds the conditions.</p>",
         if (length(no_condition_ids) > 0) paste0(
           "<div class=\"graybox\"><p><strong>", length(no_condition_ids), " sample", if (length(no_condition_ids) == 1) " has" else "s have",
           " no value in it and ", if (length(no_condition_ids) == 1) "is" else "are", " treated as the condition 'NA':</strong> ",
           names_html(no_condition_ids), ".</p></div>") else "",
         "</li>")
} else ""
item_features_to_flag <- if (use_list) {
  not_found_box <- if (length(list_not_found) == 0) "" else paste0(
    "<div class=\"graybox\"><p><strong>", length(list_not_found), if (length(list_not_found) == 1) " name is" else " names are",
    " not in the data:</strong> ", names_html(list_not_found), ".</p></div>")
  paste0("<li><p><code>--features_to_flag</code>: a list of feature names. The listed features are flagged for removal ",
         "(TRUE), all other features are kept (FALSE).</p>", not_found_box, "</li>")
} else ""
item_filtered_col <- if (featurefilt != "Manual" || use_list) {
  ""
} else {
  paste0("<li><p>", sel("--filteredfeatures_col", filteredfeatures_col, "FilteredFeaturesFlag"),
         ": the column of the ", feat_meta_txt, " that marks the features to remove. TRUE (or T, 1, yes) flags ",
         "a feature for removal, FALSE (or F, 0, no, or an empty value) keeps it.</p></li>")
}

params_html <- paste0(
  "<ul>", item_conditions_col, item_filtered_col, item_features_to_flag, item_featurefilt, item_cutoff,
  "</ul>",
  "<p>The filtering is done by <code>MetaProViz::feature_filtering()</code>, see <a href=\"https://saezlab.github.io/MetaProViz/reference/feature_filtering.html\">",
  "https://saezlab.github.io/MetaProViz/reference/feature_filtering.html</a>.</p>")

description_html <- paste0(intro_html, params_html)

# ── NOISE OR BIOLOGY (depends on the filtering option) ───────────────────────
noise_html <- if (featurefilt == "Manual") {
  paste0(
    "<p>A feature that is missing in many samples is not automatically bad data. It can be technical noise or a ",
    "real biological difference. <code>--featurefilt Manual</code> does not check this. It flags exactly the features ",
    "that were selected. The decision is therefore made before the run, when the features are selected. The table above shows ",
    "in how many samples each flagged feature is missing, which helps to check the selection.</p><p>",
  "Example from the intracellular cell-line example data of MetaProViz (HK2 healthy kidney cells and the ",
  "786-O, 786-M1A and 786-M2A renal cancer cell lines, 12 samples each): SAICAR, an intermediate of de novo purine synthesis, is ",
  "missing in all 12 HK2 samples and in 8 of 12 786-M2A samples, but detected in every 786-M1A and 786-O sample. ",
  "That pattern follows the cell line. If the measurement itself is sound, it points to a systematic difference ",
  "between the cell lines (the metabolite may be absent there, or present only at a level the method cannot ",
  "detect) rather than to random technical noise.",
    " With <code>--featurefilt Manual</code> such a feature stays in the data as long as it is marked FALSE.</p>")
} else {
  paste0(
    "<p>A feature that is missing in many samples is not automatically bad data. Technical noise, for example a ",
    "signal close to the detection limit or a badly integrated peak, tends to drop out here and there, in samples ",
    "of every condition. A real biological difference follows the conditions: the metabolite is present in some ",
    "and absent in others. The cutoff only counts missing values",
    if (featurefilt == "Standard") " across all samples and does not look at the conditions" else "",
    ". <code>--featurefilt ", featurefilt, "</code> therefore cannot tell the two apart. Look at the flagged features ",
    "before accepting the result.</p><p>",
  "Example from the intracellular cell-line example data of MetaProViz (HK2 healthy kidney cells and the ",
  "786-O, 786-M1A and 786-M2A renal cancer cell lines, 12 samples each): SAICAR, an intermediate of de novo purine synthesis, is ",
  "missing in all 12 HK2 samples and in 8 of 12 786-M2A samples, but detected in every 786-M1A and 786-O sample. ",
  "That pattern follows the cell line. If the measurement itself is sound, it points to a systematic difference ",
  "between the cell lines (the metabolite may be absent there, or present only at a level the method cannot ",
  "detect) rather than to random technical noise.",
    if (featurefilt == "Modified") paste0(
      " With the default cutoff of 0.8, Modified flags SAICAR anyway, since in HK2 it is missing in more samples ",
      "than the cutoff allows. The per-condition columns of the table above show whether missing values cluster in one ",
      "condition (likely biology) or are spread over all conditions (likely technical).")
    else paste0(
      " With the default cutoff of 0.8, Standard flags SAICAR as well, since across all 48 samples it is missing ",
      "in 20, which is more than the cutoff allows. Standard cannot show whether missing values cluster in one ",
      "condition. For that, use <code>--featurefilt Modified</code>, whose table has one column per condition."),
    " To keep a feature like SAICAR, use <code>--featurefilt Manual</code>.</p>")
}

# ── CITATIONS (the tool, the method that was used, the software around it, the workflow) ──
method_ref_html <- switch(featurefilt,
  Modified = paste0(
    "<li><strong>Modified filtering (<code>--featurefilt Modified</code>):</strong> ",
    "<a href='https://doi.org/10.3389/fmolb.2015.00004'>Yang et al., A data preprocessing strategy for ",
    "metabolomics to reduce the mask effect in data analysis, Frontiers in Molecular Biosciences 2015.</a></li>"),
  Standard = paste0(
    "<li><strong>Standard filtering (<code>--featurefilt Standard</code>):</strong> ",
    "<a href='https://doi.org/10.1021/ac051495j'>Bijlsma et al., Large-scale human metabolomics studies: a ",
    "strategy for data (pre-) processing and validation, Analytical Chemistry 2006.</a></li>"),
  "")

citations_html <- paste0(
  "<p>If you use results from this analysis in a publication, please cite:</p><ul>",
  "<li><strong>MetaProViz:</strong> ",
  "<a href='https://doi.org/10.1038/s44320-026-00231-8'>Schmidt et al., Integrated metabolomics data analysis ",
  "to generate mechanistic hypotheses with MetaProViz, Molecular Systems Biology 2026.</a></li>",
  method_ref_html,
  "<li><strong>Dependencies:</strong> ",
  "<a href='https://doi.org/10.32614/R.manuals'>R Core Team, R: A Language and Environment for Statistical Computing, ",
  "R Foundation for Statistical Computing</a>, ",
  "<a href='https://doi.org/10.18129/B9.bioc.SummarizedExperiment'>Morgan et al., SummarizedExperiment: A ",
  "container (S4 class) for matrix-like assays, Bioconductor</a> and ",
  "<a href='https://doi.org/10.1038/nmeth.3252'>Huber et al., Orchestrating high-throughput genomic analysis with Bioconductor, Nature Methods 2015</a>. ",
  "The other R packages used are listed with their ",
  "versions in the session summary at the end of this report.</li>",
  "<li><strong>This module:</strong> <code>metaproviz/processingfilter</code>, an nf-core module. It ",
  "has no publication of its own and is cited through the nf-core publication below. In the methods, ",
  "name the module.</li>",
  "<li><strong>nf-core:</strong> <a href='https://doi.org/10.1038/s41587-020-0439-x'>Ewels et al., The ",
  "nf-core framework for community-curated bioinformatics pipelines, Nature Biotechnology ",
  "2020.</a></li>",
  "<li><strong>Nextflow:</strong> <a href='https://doi.org/10.1038/nbt.3820'>Di Tommaso et al., ",
  "Nextflow enables reproducible computational workflows, Nature Biotechnology 2017.</a></li>",
  "</ul>")

# ── SESSION SUMMARY (packages used + citation reminder) ─────────────────────
session_info_text <- paste(capture.output(sessionInfo()), collapse = "\\n")
session_summary_html <- sprintf(
'<details>
  <summary>Full R session info (click to expand)</summary>
  <pre>%s</pre>
</details>
<p>Please cite MetaProViz, base R and the dependencies listed here and in the Citations section
if you use results from this report.</p>',
esc(session_info_text))

html <- sprintf('<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Processing: Feature Filtering Report</title>
<style>
  body{font-family:Arial,sans-serif;max-width:980px;margin:2em auto;padding:0 1em;line-height:1.5}
  h1{color:#1F3A5F;border-bottom:2px solid #1F3A5F;padding-bottom:.3em}
  h2{color:#1F3A5F;margin-top:2em}
  dl{display:grid;grid-template-columns:max-content 1fr;gap:.3em 1em}
  dt{font-weight:bold;color:#345}
  code{background:#f2f2f2;padding:.15em .35em;border-radius:3px}
  .graybox{background:#f2f2f2;border-left:4px solid #999;padding:.1em 1em;margin:1em 0;border-radius:3px}
  li{margin:.25em 0}
  table{border-collapse:collapse;width:100%%;font-size:.9em}
  th,td{border:1px solid #ccc;padding:.35em .6em;text-align:left}
  th{background:#f0f4f8}
  .scroll-table{max-height:320px;overflow-y:auto;border:1px solid #cbd6e4;border-radius:6px;margin:.4em 0}
  .scroll-table table{border:0;margin:0}
  .scroll-table thead th{position:sticky;top:0;background:#f0f4f8;z-index:1}
</style></head><body>
<h1>Processing: Feature Filtering</h1>
<p><em>Generated by the METAPROVIZ_PROCESSINGFILTER module.</em></p>

<h2>Run summary</h2>%s

<h2>Citations</h2>
%s

<h2>Description</h2>
%s

<h2>Features flagged for removal</h2>
<p>%s</p>
%s
<p>
%s This table is saved as <code>%s</code>.
</p>
<p>
No feature is removed from the output. The features listed here are TRUE in the added column
<code>%s</code> of the %s, all other features are FALSE. The rest of the data stays unchanged.
</p>

<h2>Noise or biology?</h2>
%s

<h2>Run notes (log)</h2>
<p>A separate file, <code>%s</code>, records every step of this run,
including any warnings. Check it first if a result here looks unexpected.</p>

<h2>Session summary</h2>
<p>Full R session used to generate this report — every package, and its version:</p>
%s

</body></html>',
run_meta_html, citations_html, description_html, flagged_features_lead, flagged_features_html,
flagged_features_note, esc(flagged_features_out), esc(filteredfeatures_col), out_feat_meta_txt, noise_html, esc(log_out), session_summary_html)

writeLines(html, report_out)
log_msg("INFO", "Written: ", report_out)

# MetaProViz and OmnipathR write their own log files and plot dumps to the working directory. They are no outputs of
# this module and are deleted here, at the very end.
unlink("metaproviz-log", recursive = TRUE)
unlink("omnipathr-log", recursive = TRUE)
unlink(Sys.glob("Rplots*.pdf"))

flush_log()
message(sprintf(
  "\\n✓ processing_filter complete — %d / %d feature(s) flagged for removal | %d warning(s) | %d error(s)",
  length(removed), n_features_input, warn_count, err_count))
