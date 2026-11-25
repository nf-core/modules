#!/usr/bin/env Rscript
# =============================================================================
# metaproviz_poolestimation.R: run MetaProViz::pool_estimation on a SummarizedExperiment
# to assess metabolite detection quality from pool samples (mixtures of all
# experimental samples that are measured several times during the run).
#
# Pool samples can be identified in three ways (top precedence first):
#   1. Explicit list            --pool_samples POOL1,POOL2,...
#   2. Metadata-column lookup   --pool_metadata_col Conditions
#                               --pool_metadata_value Pool  (default value: Pool, upper and lower case are ignored)
#   3. Pattern on sample names  --pool_pattern "pool|qc"  (default; case-insensitive)
# If the list or the column finds no pool sample, the next way is tried, and the pattern is the last one.
#
# The format of the SE output is the val input `output_format` of the module (not an option in ext.args, since it
# decides which outputs exist): se (default, also if empty) or tsv or csv.
#
# Options (set via task.ext.args):
#   --cutoff_cv        CV cutoff for flagging high-variance features, passed to
#                      MetaProViz::pool_estimation(). Default: 30.
#   --filter_high_var  TRUE: drop HighVar features from the SE (default).
#                      FALSE: keep them and only report.
#
# Outputs (all named <prefix>.<suffix>; prefix = task.ext.prefix or meta.id):
#   <prefix>.pool_filtered.rds  SE without the pool samples, with HighVar features dropped (if filtering is on).
#                               With the val input output_format tsv/csv: <prefix>.pool_filtered.{data,feature,sample}_matrix.<fmt>
#                               instead (three tables, the same shape as the input matrices)
#   <prefix>.pools.txt          names of the pool samples, one per line
#   <prefix>.cv.tsv             full CV table, one row per feature, with the HighVar flag
#   <prefix>.plots.rds          named list of ggplot objects
#   <prefix>.report.html        self-contained HTML report (plots embedded)
#   <prefix>.log                INFO/WARN log
# =============================================================================

set.seed(42)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(magrittr)
  library(ggplot2)
  library(SummarizedExperiment)
  library(S4Vectors)
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
  pool_pattern        = "pool|qc",
  pool_samples        = character(0),
  pool_metadata_col   = "",
  pool_metadata_value = "Pool",
  cutoff_cv           = "30",
  filter_high_var     = "TRUE"
)
# R raw string: quotes and apostrophes inside the arguments (feature names like 5'-AMP) cannot break the R code.
args_opt <- parse_args(r"---($task.ext.args)---", list_opts = c("pool_samples"))
for (ao in names(args_opt)) {
  if (!ao %in% names(opt)) stop(paste("Invalid option:", ao))
  opt[[ao]] <- args_opt[[ao]]
}

se_path             <- nullify('$se_rds')
data_matrix_path    <- nullify('$data_matrix')
feature_matrix_path <- nullify('$feature_matrix')
sample_matrix_path  <- nullify('$sample_matrix')
pool_pattern        <- opt\$pool_pattern
pool_samples_list <- opt\$pool_samples
pool_samples_str    <- paste(pool_samples_list, collapse = ",")
pool_metadata_col   <- opt\$pool_metadata_col
pool_metadata_value <- opt\$pool_metadata_value
cutoff_cv           <- as.numeric(opt\$cutoff_cv)
filter_high_var     <- toupper(opt\$filter_high_var) == "TRUE"

prefix              <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
# output_format is a val input of the module. se: the SE as one RDS. tsv/csv: three files (data, feature and sample
# matrix), the same shape as the input matrices. An empty value means se.
output_format_in    <- trimws('$output_format')
output_format       <- if (output_format_in %in% c("", "null", "[]")) "se" else tolower(output_format_in)
se_out              <- if (output_format == "se") paste0(prefix, ".pool_filtered.rds") else
  paste0(prefix, ".pool_filtered", c(".data_matrix.", ".feature_matrix.", ".sample_matrix."), output_format)
# the output names in running text
se_out_txt          <- if (length(se_out) == 1) se_out else paste0(paste(head(se_out, -1), collapse = ", "), " and ", tail(se_out, 1))
pools_out           <- paste0(prefix, ".pools.txt")
cv_out              <- paste0(prefix, ".cv.tsv")
plots_out           <- paste0(prefix, ".plots.rds")
report_out          <- paste0(prefix, ".report.html")
log_out             <- paste0(prefix, ".log")

# ── LOGGING ─────────────────────────────────────────────────────────────────

log_lines  <- character(0)
warn_count <- 0; err_count <- 0

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
               sprintf("SUMMARY: %d warning(s), %d error(s)",
                       warn_count, err_count)), path)
}
abort <- function(...) { log_msg("ERROR", ...); flush_log(); stop(paste0(...), call. = FALSE) }
if (!output_format %in% c("se", "tsv", "csv"))
  abort("output_format must be se, tsv or csv (empty means se), not '", output_format_in, "'.")

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

# Used in the report/log wherever the input source needs naming; se_path is
# NULL when the module was given the 3 flat files instead of an .rds.
source_label <- if (have_se) {
  basename(se_path)
} else {
  paste(basename(data_matrix_path), basename(feature_matrix_path),
        basename(sample_matrix_path), sep = " + ")
}

# Data layout for MetaProViz: samples x features (rownames = sample IDs)
assay_df <- t(assay(se, 1)) %>% as.data.frame(check.names = FALSE)
sample_info_df <- as.data.frame(colData(se), check.names = FALSE)

# ── IDENTIFY POOL SAMPLES ───────────────────────────────────────────────────

log_section("Identifying pool samples")

pool_ids <- character(0)
strategy <- "none"
attempted <- character(0)   # ways of finding the pools that were tried and found no pool sample

# 1) Explicit list
if (nzchar(pool_samples_str)) {
  asked <- trimws(pool_samples_list)
  asked <- asked[nzchar(asked)]
  pool_ids <- intersect(asked, rownames(assay_df))
  missing  <- setdiff(asked, rownames(assay_df))
  if (length(missing) > 0)
    log_msg("WARN", "Pool sample(s) not found in SE: ",
            paste(missing, collapse = ", "))
  strategy <- "explicit_list"
  if (length(pool_ids) == 0) attempted <- c(attempted, list = "explicit list")
}

# 2) Metadata-column lookup
if (length(pool_ids) == 0 && nzchar(pool_metadata_col)) {
  if (!pool_metadata_col %in% colnames(sample_info_df)) {
    log_msg("WARN", "pool_metadata_col '", pool_metadata_col,
            "' not in colData. Available: ",
            paste(colnames(sample_info_df), collapse = ", "))
    attempted <- c(attempted, column = sprintf("metadata column '%s' with the value '%s'", pool_metadata_col, pool_metadata_value))
  } else {
    hit <- which(tolower(as.character(sample_info_df[[pool_metadata_col]])) == tolower(pool_metadata_value))
    pool_ids <- rownames(sample_info_df)[hit]
    strategy <- "metadata_column"
    if (length(pool_ids) == 0)
      attempted <- c(attempted, column = sprintf("metadata column '%s' with the value '%s'", pool_metadata_col, pool_metadata_value))
  }
}

# 3) Pattern on rownames (default fallback)
if (length(pool_ids) == 0 && nzchar(pool_pattern)) {
  hit <- grepl(pool_pattern, rownames(assay_df), ignore.case = TRUE)
  pool_ids <- rownames(assay_df)[hit]
  strategy <- paste0("pattern '", pool_pattern, "'")
}

# Wording of the detection strategy for messages and the report.
not_given_label <- if (strategy == "metadata_column") {
  sprintf("do not have the value '%s' in the metadata column '%s'", pool_metadata_value, pool_metadata_col)
} else "are not part of the given explicit list"
strategy_label <- if (strategy == "explicit_list") "explicit list" else if (strategy == "metadata_column") {
  sprintf("metadata column '%s' with the value '%s'", pool_metadata_col, pool_metadata_value)
} else strategy

if (length(pool_ids) == 0) {
  abort("No pool samples were found. Please check the pool detection parameters of this module.")
}
if (length(pool_ids) < 2) {
  abort("Only ", length(pool_ids), " pool sample found (", paste(pool_ids, collapse = ", "),
        ", found by ", strategy_label, "). At least two pool samples are needed to compute a CV. ",
        "Please check the pool detection parameters of this module.")
}

if (length(pool_ids) >= ncol(se)) {
  abort("All ", ncol(se), " samples were identified as pool samples (found by ", strategy_label, "). Study samples are needed ",
        "next to the pool samples, since the PCA compares the two and the output holds the study samples. ",
        "Please check the pool detection parameters of this module.")
}

# If the list or the column found no pool sample, the pattern was used instead. Say so in the log and the report.
attempted_txt <- ""
if (length(attempted) > 0 && startsWith(strategy, "pattern")) {
  attempted_txt <- paste(attempted, collapse = " and the ")
  log_msg("WARN", "No pool samples were found by the ", attempted_txt, ". The pool detection falls back to the ",
          strategy_label, ".")
}

log_msg("INFO", "Pool detection strategy: ", strategy)
log_msg("INFO", "Pool samples (", length(pool_ids), "): ",
        paste(head(pool_ids, 12), collapse = ", "),
        if (length(pool_ids) > 12) " ..." else "")

# An explicit list or a metadata column takes precedence over the name pattern. Warn if samples that match the
# pattern are left in the data, so a list that misses some pools does not go unnoticed.
pool_leftover_ids <- character(0)
if (!startsWith(strategy, "pattern") && nzchar(pool_pattern)) {
  pool_leftover_ids <- setdiff(rownames(assay_df)[grepl(pool_pattern, rownames(assay_df), ignore.case = TRUE)], pool_ids)
  if (length(pool_leftover_ids) > 0) {
    log_msg("WARN", length(pool_leftover_ids), " sample(s) match the pool pattern '", pool_pattern, "' as well, but ",
            not_given_label, ": ", paste(pool_leftover_ids, collapse = ", "),
            ". In this run they are treated as normal samples and are kept in ", se_out_txt, ".")
  }
}

# ── RUN POOL ESTIMATION ─────────────────────────────────────────────────────

log_section("Running MetaProViz::pool_estimation")

# MetaProViz::pool_estimation() finds the pool samples by the value "Pool" in a column named Conditions of the sample
# table. This module finds the pools in other ways (sample name, explicit list, metadata column), so it passes them on
# in this form. It does this on a copy of the sample table. The input and the output data are not changed.
#   Conditions column exists, pool rows already say "Pool":  nothing changes
#   Conditions column exists, pools found by name only:      the pool rows are set to "Pool" for the call
#   no Conditions column:                                    one is created for the call ("Sample" for all rows,
#                                                            "Pool" for the pool rows)
#   a non-pool sample already says "Pool":                   MetaProViz would take it as a pool, so it is set to
#                                                            "Sample" for the call, with a warning

pool_label <- "Pool"

if (!"Conditions" %in% colnames(sample_info_df)) {
  log_msg("INFO", "No 'Conditions' column found in colData — creating one.")
  sample_info_df\$Conditions <- "Sample"
} else {
  sample_info_df\$Conditions_original <- sample_info_df\$Conditions
}

is_pool_row    <- rownames(sample_info_df) %in% pool_ids
no_condition   <- !is_pool_row & (is.na(sample_info_df\$Conditions) | trimws(as.character(sample_info_df\$Conditions)) == "")
if (any(no_condition)) {
  log_msg("WARN", sum(no_condition), " sample(s) have no value in the Conditions column: ",
          paste(rownames(sample_info_df)[no_condition], collapse = ", "), ". They are therefore treated as normal samples and kept in ", se_out_txt, ".")
}
already_pool   <- (sample_info_df\$Conditions == pool_label & !is_pool_row) %in% TRUE
if (any(already_pool, na.rm = TRUE)) {
  log_msg("WARN", sum(already_pool), " sample(s) have '", pool_label, "' in their Conditions column but were not ",
          "detected as pool samples: ", paste(rownames(sample_info_df)[already_pool], collapse = ", "),
          ". MetaProViz would treat them as pools, and their Conditions are set to 'Sample' while it runs. ",
          "The output keeps the original Conditions.")
  sample_info_df\$Conditions[already_pool] <- "Sample"
}

sample_info_df\$Conditions[is_pool_row] <- pool_label
log_msg("INFO", "Stamped Conditions = '", pool_label, "' on ",
        sum(is_pool_row), " pool sample row(s).")

# Conditions going into pool_estimation()
log_msg("INFO", "Conditions levels going into pool_estimation: ",
        paste(sprintf("%s=%d", names(table(sample_info_df\$Conditions)),
                      as.integer(table(sample_info_df\$Conditions))),
              collapse = ", "))

pe <- tryCatch(
  MetaProViz::pool_estimation(
    data            = assay_df,
    metadata_sample = sample_info_df,
    metadata_info   = c(PoolSamples = pool_label, Conditions = "Conditions"),
    cutoff_cv       = cutoff_cv,
    print_plot      = FALSE,
    save_plot       = NULL
  ),
  error = function(e) abort("pool_estimation() failed: ", conditionMessage(e))
)

# pool_estimation() writes its own result files and logs to the working directory. They are no outputs of this
# module and are deleted here.
unlink("MetaProViz_Results", recursive = TRUE)
unlink("metaproviz-log", recursive = TRUE)
unlink("omnipathr-log", recursive = TRUE)
unlink("Rplots.pdf")

# ── EXTRACT RESULTS ─────────────────────────────────────────────────────────

log_section("Extracting results")

cv_df <- pe[["DF"]][["CV"]]
if (is.null(cv_df))
  abort("pool_estimation() did not return DF\$CV — unexpected result shape.")
log_msg("INFO", "CV table: ", nrow(cv_df), " feature(s), columns: ",
        paste(colnames(cv_df), collapse = ", "))

high_var <- cv_df %>% filter(HighVar == TRUE) %>% pull(Metabolite)
log_msg("INFO", "High-variance features at CV > ", cutoff_cv,
        "%: ", length(high_var), " / ", nrow(cv_df))

plots <- pe[["Plot"]]
if (is.null(plots) || length(plots) == 0) {
  log_msg("WARN", "pool_estimation() did not return any plots.")
} else {
  log_msg("INFO", "Plots returned: ", paste(names(plots), collapse = ", "))
}

# ── OPTIONAL FILTERING ──────────────────────────────────────────────────────

log_section("Filtering high-variance features (optional)")

if (filter_high_var && length(high_var) > 0) {
  keep <- setdiff(rownames(se), high_var)
  log_msg("INFO", "Removing ", length(high_var),
          " high-variance feature(s); keeping ", length(keep), ".")
  se_filtered <- se[keep, ]
  metadata(se_filtered)\$pool_estimation <- list(
    cutoff_cv      = cutoff_cv,
    pool_samples   = pool_ids,
    high_var_feats = high_var,
    strategy       = strategy
  )
} else {
  se_filtered <- se
  if (filter_high_var) {
    log_msg("INFO", "No HighVar features to remove — SE unchanged.")
  } else {
    log_msg("INFO", "filter_high_var=FALSE — SE unchanged.")
  }
  metadata(se_filtered)\$pool_estimation <- list(
    cutoff_cv    = cutoff_cv,
    pool_samples = pool_ids,
    strategy     = strategy
  )
}

# ── SPLIT OFF THE POOL SAMPLES ──────────────────────────────────────────────

log_section("Separating the pool samples")

# The pool samples were only needed for this quality control. They are mixtures of the study samples and not
# biological samples, so they are taken out of the output. Their names are written to a separate file.
se_filtered <- se_filtered[, setdiff(colnames(se_filtered), pool_ids), drop = FALSE]
log_msg("INFO", "Pool samples removed from the output: ", length(pool_ids), ", samples kept: ", ncol(se_filtered), ".")
log_msg("INFO", "Pool sample names: ", paste(pool_ids, collapse = ", "))

# ── WRITE OUTPUTS ───────────────────────────────────────────────────────────

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
  if (length(assays(se_filtered)) > 1)
    log_msg("WARN", "The input has ", length(assays(se_filtered)), " assays. Only the first one is written to the data matrix.")
  log_msg("INFO", "The metadata of the SE (cutoff, pool samples and pool detection strategy) is not part of the tables, ",
          "it is in this log and in the report.")
}
write_out(se_filtered, se_out);                     log_msg("INFO", "Written: ", paste(se_out, collapse = ", "))
writeLines(pool_ids, pools_out);                    log_msg("INFO", "Written: ", pools_out)

# The CV is written with 10 significant digits.
cv_tsv_df <- cv_df
cv_tsv_df\$CV <- signif(as.numeric(cv_tsv_df\$CV), 10)
write.table(cv_tsv_df, cv_out,
            sep = "\\t", quote = FALSE, row.names = FALSE, na = "NA")
log_msg("INFO", "Written: ", cv_out)

saveRDS(plots, plots_out)
log_msg("INFO", "Written: ", plots_out, " (a named list of ggplots)")

# ── HTML REPORT ─────────────────────────────────────────────────────────────
# Self-contained: ggplots >>> PNG >>> base64 >>> embedded <img>.

log_section("Building HTML report")

html_escape <- function(x) {
  x <- gsub("&", "&amp;",  x, fixed = TRUE)
  x <- gsub("<", "&lt;",   x, fixed = TRUE)
  x <- gsub(">", "&gt;",   x, fixed = TRUE)
  x
}
# Names for the report: at most 25, and a pointer to the file that lists all of them only if not all are shown.
names_html <- function(x, ref_file, max_n = 25) paste0(
  paste0("<code>", html_escape(head(x, max_n)), "</code>", collapse = " "),
  if (length(x) > max_n) paste0(" and ", length(x) - max_n, " more (see <code>", html_escape(ref_file), "</code>)") else "")
# the SE output file(s) in the report: <code>a</code>, <code>b</code> and <code>c</code>
files_html <- function(f) if (length(f) == 1) paste0("<code>", html_escape(f), "</code>") else
  paste0(paste0("<code>", html_escape(head(f, -1)), "</code>", collapse = ", "), " and <code>", html_escape(tail(f, 1)), "</code>")
plot_to_b64 <- function(p, width = 8, height = 5, dpi = 100) {
  tmp <- tempfile(fileext = ".png")
  ggplot2::ggsave(tmp, plot = p, width = width, height = height,
                  dpi = dpi, units = "in", bg = "white")
  raw <- readBin(tmp, what = "raw", n = file.info(tmp)\$size)
  unlink(tmp)
  paste0("data:image/png;base64,",
         base64enc::base64encode(raw))
}
df_to_html_table <- function(df, max_rows = 5, source_file = NULL) {
  # Only the top `max_rows` rows are shown, the full data is in the TSV next to the HTML. The table sits in a
  # scrollable container with a sticky header.
  rows_shown <- head(df, max_rows)
  th <- paste0("<th>", html_escape(colnames(rows_shown)), "</th>", collapse = "")
  body_rows <- apply(rows_shown, 1, function(r) {
    paste0("<tr>",
           paste0("<td>", html_escape(as.character(r)), "</td>", collapse = ""),
           "</tr>")
  })
  more <- if (nrow(df) > max_rows) {
    src_hint <- if (!is.null(source_file))
      sprintf(" (full data in %s)", html_escape(source_file)) else ""
    sprintf("<p class='more-rows'><em>Showing first %d of %d row(s)%s.</em></p>",
            max_rows, nrow(df), src_hint)
  } else {
    ""
  }
  paste0(
    '<div class="scroll-table"><table><thead><tr>', th,
    "</tr></thead><tbody>",
    paste(body_rows, collapse = ""), "</tbody></table></div>", more
  )
}

# Build embedded plot section
plot_html <- ""
if (!is.null(plots) && length(plots) > 0) {
  for (nm in names(plots)) {
    p <- plots[[nm]]
    if (!inherits(p, "ggplot")) next
    src <- tryCatch(plot_to_b64(p),
                    error = function(e) {
                      log_msg("WARN", "Could not render plot '", nm, "': ",
                              conditionMessage(e)); NA_character_
                    })
    if (!is.na(src)) {
      plot_html <- paste0(
        plot_html,
        sprintf('<h3>%s</h3>\\n<img src="%s" alt="%s" />\\n',
                html_escape(nm), src, html_escape(nm))
      )
    }
  }
}

# The ways that found no pool sample, written with the names of their parameters.
attempted_html <- paste(vapply(names(attempted), function(k) if (k == "list") {
  "explicit list <code>--pool_samples</code>"
} else {
  sprintf("metadata column <code>--pool_metadata_col %s</code> with the value <code>--pool_metadata_value %s</code>",
          html_escape(pool_metadata_col), html_escape(pool_metadata_value))
}, character(1)), collapse = " and the ")
# The detection strategy written with the names of its parameters.
strategy_html <- if (strategy == "explicit_list") {
  "explicit list <code>--pool_samples</code>"
} else if (strategy == "metadata_column") {
  sprintf("metadata column <code>--pool_metadata_col %s</code> with the value <code>--pool_metadata_value %s</code>",
          html_escape(pool_metadata_col), html_escape(pool_metadata_value))
} else sprintf("pattern <code>--pool_pattern %s</code>", html_escape(pool_pattern))
not_given_html <- if (strategy == "metadata_column") {
  sprintf("do not have the value <code>--pool_metadata_value %s</code> in the metadata column <code>--pool_metadata_col %s</code>",
          html_escape(pool_metadata_value), html_escape(pool_metadata_col))
} else "are not part of the given explicit list <code>--pool_samples</code>"
# CV and MissingValuepercentage columns (both in percent): rounded to one decimal, the unit is named in the text above
# the table. The saved TSV keeps the plain numbers.
cv_display_df <- as.data.frame(lapply(cv_df, function(x) x), stringsAsFactors = FALSE)
for (col in intersect(c("CV", "MissingValuepercentage"), colnames(cv_display_df))) {
  v <- suppressWarnings(as.numeric(cv_df[[col]]))
  cv_display_df[[col]] <- ifelse(is.na(v), "—", format(round(v, 1)))
}
cv_display_df[] <- lapply(cv_display_df, function(x) ifelse(is.na(x), "—", as.character(x)))
high_var_html <- if (length(high_var) > 0) {
  # more than 10 features: the list scrolls
  paste0(if (length(high_var) > 10) '<div class="scroll-list">' else "<div>",
         "<ul>", paste0("<li>", html_escape(high_var), "</li>", collapse = ""), "</ul></div>")
} else "<p><em>None — all features below the CV cutoff.</em></p>"

n_features_txt <- sprintf("%d feature%s, %d flagged as HighVar", nrow(se), if (nrow(se) == 1) "" else "s", length(high_var))
n_samples_txt  <- sprintf("%d sample%s, %d pool sample%s", ncol(se_filtered), if (ncol(se_filtered) == 1) "" else "s",
                          length(pool_ids), if (length(pool_ids) == 1) "" else "s")
filter_txt     <- sprintf("%s, total of %d feature%s kept", if (filter_high_var) "TRUE" else "FALSE", nrow(se_filtered),
                          if (nrow(se_filtered) == 1) "" else "s")
run_meta_html <- sprintf(
  "<dl>
     <dt>Input source</dt><dd>%s</dd>
     <dt>n features</dt><dd>%s</dd>
     <dt>n samples</dt><dd>%s</dd>
     <dt>Pool detection strategy</dt><dd>%s</dd>
     <dt>CV cutoff</dt><dd>%s%%</dd>
     <dt>Filter HighVar from SE</dt><dd>%s</dd>
   </dl>",
  html_escape(source_label),
  n_features_txt, n_samples_txt,
  if (nzchar(attempted_txt)) paste0(strategy_html, " (fallback, since no pool samples were found by the ",
                                    attempted_html, ")") else strategy_html,
  format(cutoff_cv),
  filter_txt
)

# Sentence end after the default cutoff of 30%: nothing more for the default, or the cutoff set with --cutoff_cv.
cutoff_note_html <- if (cutoff_cv == 30) "." else paste0(". This run used <code>--cutoff_cv ", format(cutoff_cv), "</code>.")

# Why the HighVar features should be removed, then that the module does it by default, or that this run kept them.
filter_note_html <- paste0(
  "Features above the cutoff should be removed from the output data to keep only reliable features for the downstream ",
  "analysis. ",
  if (filter_high_var) "This module removes them by default (<code>--filter_high_var TRUE</code>)." else
    paste0("Since keeping them was selected (<code>--filter_high_var FALSE</code>), they are flagged as HighVar and kept in ",
           "this run. Use the CV table and the plots to check whether their removal is wanted."))

# ── POOL SAMPLES IN THE REPORT ──────────────────────────────────────────────
pool_section_html <- paste0(
  "<h2>Pool samples</h2>\\n<p>", length(pool_ids), " pool ", if (length(pool_ids) == 1) "sample was" else "samples were",
  " found by ", strategy_html, ":</p>\\n<p>", names_html(pool_ids, pools_out),
  "</p>\\n<p>The pool samples are needed only for this quality control. As a mixture of all study samples, they do not ",
  "represent the biology of a single sample or condition and add no information for the downstream analysis. They are ",
  "therefore removed from the output ", if (length(se_out) == 1) "" else "files ", files_html(se_out), ". Their names are saved in <code>", html_escape(pools_out),
  "</code>.</p>\\n",
  if (length(pool_leftover_ids) > 0) paste0(
    "<p><strong>Additionally matching the pool pattern:</strong> ",
    names_html(pool_leftover_ids, log_out), ". These samples match the pattern <code>--pool_pattern ",
    html_escape(pool_pattern), "</code> as well, but ", not_given_html, ". In this run they are therefore treated as normal samples and are kept in ",
    files_html(se_out), ".</p>\\n") else "")

# ── SESSION SUMMARY ─────────────────────────────────────────────────────────
session_info_text <- paste(capture.output(sessionInfo()), collapse = "\\n")
session_summary_html <- sprintf(
'<details>
  <summary>Full R session info (click to expand)</summary>
  <pre>%s</pre>
</details>
<p>Please cite MetaProViz, base R and the dependencies listed here and in the Citations section
if you use results or plots from this report.</p>',
html_escape(session_info_text))

html <- sprintf('<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>Pool Estimation Report — %s</title>
<style>
  body  { font-family: Arial, sans-serif; max-width: 980px; margin: 2em auto;
          color: #1e1e1e; line-height: 1.5; padding: 0 1em; }
  h1    { color: #1F3A5F; border-bottom: 2px solid #1F3A5F; padding-bottom: .3em; }
  h2    { color: #1F3A5F; margin-top: 2em; }
  h3    { color: #345; margin-top: 1.4em; }
  dl    { display: grid; grid-template-columns: max-content 1fr; gap: .3em 1em; }
  dt    { font-weight: bold; color: #345; }
  table { border-collapse: collapse; width: 100%%; font-size: .9em; }
  th,td { border: 1px solid #ccc; padding: .35em .6em; text-align: left; }
  th    { background: #f0f4f8; }
  tr:nth-child(even) td { background: #fafafa; }
  img   { max-width: 100%%; border: 1px solid #ddd; padding: 4px;
          background: #fff; margin: .4em 0 1.2em; }
  code  { background: #f2f2f2; padding: .15em .35em; border-radius: 3px; }
    pre   { background: #f7f9fc; border: 1px solid #d9e2ec; border-radius: 6px;
      padding: .7em .9em; overflow-x: auto; }
    li    { margin: .25em 0; }
  ul    { margin-top: .3em; }
  /* Scrollable table container with sticky header */
  .scroll-table       { max-height: 320px; overflow-y: auto;
                        border: 1px solid #cbd6e4; border-radius: 6px;
                        margin: .4em 0 .6em; }
  .scroll-table table { border: 0; margin: 0; }
  .scroll-table thead th { position: sticky; top: 0; z-index: 1;
                           background: #f0f4f8; }
  /* Scrollable list of the high-variance features (more than 10) */
  .scroll-list        { max-height: 19em; overflow-y: auto; border: 1px solid #cbd6e4; border-radius: 6px;
                        padding: 0 .6em; margin: .4em 0 .6em; }
  .scroll-list ul     { margin: .4em 0; }
  .more-rows { color: #555; font-size: .9em; margin: .2em 0 1.2em; }
</style></head><body>

<h1>Pool Estimation Report</h1>
<p><em>Generated by the METAPROVIZ_POOLESTIMATION module.</em></p>

<h2>Run summary</h2>
%s

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
  <li><strong>This module:</strong> <code>metaproviz/poolestimation</code>, an nf-core module. It has no publication of its own and is cited through the nf-core publication below. In the methods, name the module.</li>
  <li><strong>nf-core:</strong> <a href="https://doi.org/10.1038/s41587-020-0439-x">Ewels et al., The nf-core framework for community-curated bioinformatics pipelines, Nature Biotechnology 2020.</a></li>
  <li><strong>Nextflow:</strong> <a href="https://doi.org/10.1038/nbt.3820">Di Tommaso et al., Nextflow enables reproducible computational workflows, Nature Biotechnology 2017.</a></li>
</ul>

<h2>Description</h2>
<p>
This report uses pooled quality-control (QC) samples to assess technical data
quality. Pool samples are a homogeneous mixture of all study samples and should
be stable across repeated measurements.
</p>
<p>
To assess this stability, the coefficient of variation (CV) of each feature
across the pool samples is calculated:
</p>
<p><code>CV(%%) = 100 × sd(x1, x2, ..., xn) / mean(x1, x2, ..., xn)</code></p>
<p>
This module computes the CVs through <code>MetaProViz::pool_estimation()</code>. The result is a feature-level
CV table and the plots shown in this report.
</p>
<p>
A low CV means that a feature is measured consistently across the pool samples, and a high CV points to low
detection accuracy. As a practical rule of thumb, pool-sample CV values are typically expected to be below 30%% (the
default cutoff of this module)%s %s
</p>

%s<h2>Plots</h2>
<p>
<strong>How to read the plots:</strong>
</p>
<ul>
  <li><strong>PCA plot:</strong> Pool samples are expected to cluster near the coordinate origin,
  reflecting their role as a homogeneous mixture of all samples. This pattern is
  most interpretable when all sample groups are included in the PCA view.</li>
  <li><strong>Histogram:</strong> x-axis = CV and y-axis = frequency. The cutoff line in this run is %s%%. Most features are expected at lower CV values (left side),
  with frequencies tapering toward higher CV values.</li>
  <li><strong>Violin plot:</strong> Distribution of metabolite CV values split into two
  groups relative to the cutoff (%s%%): features below cutoff vs features above
  cutoff (high-variance candidates).</li>
</ul>
<p>
Use these diagnostics together with the high-variance feature list and full CV
table to check the features that are flagged as HighVar.
</p>
%s

<h2>High-variance features (CV > %s%%)</h2>
%s

<h2>Full CV table</h2>
<p>
Each row shows the CV (%%) of a feature across the pool samples, its HighVar flag and its missing values (%%) among
the pool samples. HighVar is TRUE if the CV is above the cutoff of %s%%, which marks the feature as too variable, and
FALSE if not. The table is saved as <code>%s</code>.
</p>
%s

<h2>Run notes (log)</h2>
<p>A separate file, <code>%s</code>, records every step of this run,
including any warnings. Check it first if a result here looks unexpected.</p>

<h2>Session summary</h2>
<p>Full R session used to generate this report — every package, and its version:</p>
%s

</body></html>',
  html_escape(source_label),
  run_meta_html,
  cutoff_note_html,
  filter_note_html,
  pool_section_html,
  format(cutoff_cv),
  format(cutoff_cv),
  if (nzchar(plot_html)) plot_html else "<p><em>No plots returned.</em></p>",
  format(cutoff_cv),
  high_var_html,
  format(cutoff_cv),
  html_escape(cv_out),
  # The full CV table, in a scrollable container (.scroll-table).
  df_to_html_table(cv_display_df, max_rows = 5000,
                   source_file = cv_out),
  html_escape(log_out),
  session_summary_html
)
writeLines(html, report_out)
log_msg("INFO", "Written: ", report_out)

flush_log(log_out)
log_msg("INFO", "Written: ", log_out)

message(sprintf(
  "\\n✓ pool_estimation complete — %d / %d features flagged HighVar | %d warning(s) | %d error(s)",
  length(high_var), nrow(cv_df), warn_count, err_count
))
