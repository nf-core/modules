#!/usr/bin/env Rscript
# =============================================================================
# metaproviz_processingnormalize.R: total ion count (TIC) normalization step of the MetaProViz processing chain.
#
# Approach: MetaProViz::tic_norm() divides each sample by its total ion count (NA counted as 0), which corrects for
# differences in injection amount or instrument sensitivity (Wulff and Mitchell, 2018,
# https://doi.org/10.4236/abb.2018.98022). The Conditions column is optional, it only colors the RLA
# plots. Without it all samples form one group.
#
# A sample with a total ion count of 0 (every value NA or 0) would be divided by 0 and is removed, with a warning. The
# data should be imputed first (metaproviz/processingimpute), NA or 0 values are reported.
#
# Pool/QC samples must be removed before this module, they are not detected here.
#
# Outputs, all named <prefix>.<suffix>: the SE (or three tables if output_format is tsv or csv), removed_samples.tsv,
# missing_values.tsv, plots.rds, the HTML report and the log.
# =============================================================================


set.seed(42)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(magrittr)
  library(ggplot2); library(SummarizedExperiment); library(S4Vectors)
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
  conditions_col = "Conditions"
)
# R raw string: quotes and apostrophes inside the arguments (feature names like 5'-AMP) cannot break the R code.
args_opt <- parse_args(r"---($task.ext.args)---", list_opts = character(0))
invalid_opts <- setdiff(names(args_opt), names(opt))   # reported below, once the log exists
for (ao in setdiff(names(args_opt), invalid_opts)) opt[[ao]] <- args_opt[[ao]]

# The inputs of the process, empty for the shape that is not used.
se_path             <- nullify('$se_rds')
data_matrix_path    <- nullify('$data_matrix')
feature_matrix_path <- nullify('$feature_matrix')
sample_matrix_path  <- nullify('$sample_matrix')
conditions_col <- opt\$conditions_col

# Output files are named <prefix>.<suffix>.
prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
# output_format is a val input of the module. se: the SE as one RDS. tsv/csv: three files (data, feature and sample
# matrix), the same shape as the input matrices. An empty value means se.
output_format_in   <- trimws('$output_format')
output_format      <- if (output_format_in %in% c("", "null", "[]")) "se" else tolower(output_format_in)
se_out             <- if (output_format == "se") paste0(prefix, ".normalized.rds") else
  paste0(prefix, ".normalized", c(".data_matrix.", ".feature_matrix.", ".sample_matrix."), output_format)
removed_out        <- paste0(prefix, ".removed_samples.tsv")
missing_out        <- paste0(prefix, ".missing_values.tsv")
plots_out          <- paste0(prefix, ".plots.rds")
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
png_to_b64 <- function(p, w = 8, h = 5) {
  tmp <- tempfile(fileext = ".png")
  ggplot2::ggsave(tmp, plot = p, width = w, height = h, dpi = 100,
                  units = "in", bg = "white")
  raw <- readBin(tmp, "raw", n = file.info(tmp)\$size); unlink(tmp)
  paste0("data:image/png;base64,", base64enc::base64encode(raw))
}

if (length(invalid_opts) > 0)
  abort("Invalid option", if (length(invalid_opts) > 1) "s" else "", ": ", paste0("--", invalid_opts, collapse = ", "),
        ". The options are: ", paste0("--", names(opt), collapse = ", "), ".")
if (!output_format %in% c("se", "tsv", "csv"))
  abort("output_format must be se, tsv or csv (empty means se), not '", output_format_in, "'.")

# ── METAPROVIZ VERSION ──────────────────────────────────────────────────────
# tic_norm() is part of MetaProViz from version 4.1.
if (!"tic_norm" %in% getNamespaceExports("MetaProViz"))
  abort("This module needs MetaProViz >= 4.1 (exports tic_norm()), found ",
        as.character(packageVersion("MetaProViz")), ".")

# ── LOAD INPUT ─────────────────────────────────────────────────────────────────
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


# ── TIC NORMALIZATION ────────────────────────────────────────────────────────
log_section("Running MetaProViz::tic_norm()")

# tic_norm() takes the samples in rows and the features in columns.
assay_df       <- t(assay(se, 1)) %>% as.data.frame(check.names = FALSE)
sample_info_df <- as.data.frame(colData(se), check.names = FALSE)

# Each sample is divided by its total ion count (NA counted as 0, as tic_norm() does). A sample with a total of 0 would be
# divided by 0 and is removed, with a warning. If no sample is left, the run stops.
sample_totals   <- rowSums(replace(as.matrix(assay_df), is.na(as.matrix(assay_df)), 0))
zero_total_ids  <- rownames(assay_df)[sample_totals == 0]
if (length(zero_total_ids) == nrow(assay_df))
  abort("All ", nrow(assay_df), " samples have a total ion count of 0 (every value is NA or 0). Nothing is left to normalize.")
removed_df <- data.frame(Sample = zero_total_ids, stringsAsFactors = FALSE)
if (length(zero_total_ids) > 0) {
  log_msg("WARN", length(zero_total_ids), " sample(s) have a total ion count of 0 (every value is NA or 0) and were removed, ",
          "since the normalization would divide by 0: ", paste(zero_total_ids, collapse = ", "),
          ". They are not in the output.")
  keep_samples   <- !rownames(assay_df) %in% zero_total_ids
  assay_df       <- assay_df[keep_samples, , drop = FALSE]
  sample_info_df <- sample_info_df[keep_samples, , drop = FALSE]
  se             <- se[, keep_samples]
}

# Conditions only color the plots, so a missing column is not fatal: warn and give every sample the same placeholder
# condition (in the copy handed to MetaProViz only, the output SE is built from the original). Samples without a value
# (NA or empty) are one condition of their own, called "NA". Numbers are turned into text.
has_conditions   <- conditions_col %in% colnames(sample_info_df)
no_condition_ids <- character(0)
if (!has_conditions) {
  log_msg("WARN", "Conditions column '", conditions_col, "' not found in colData. ",
          "Continuing without conditions: all samples are treated as one group, ",
          "so the plots are not colored by condition.")
  if (!nzchar(conditions_col)) conditions_col <- "Conditions"
  sample_info_df[[conditions_col]] <- "All samples"
} else {
  cond_chr <- as.character(sample_info_df[[conditions_col]])
  no_cond  <- is.na(cond_chr) | trimws(cond_chr) == ""
  no_condition_ids <- rownames(sample_info_df)[no_cond]
  if (any(no_cond)) {
    cond_chr[no_cond] <- "NA"
    log_msg("WARN", length(no_condition_ids), " sample(s) have no value in the Conditions column '", conditions_col, "': ",
            paste(no_condition_ids, collapse = ", "), ". They are treated as the condition 'NA'.")
  }
  sample_info_df[[conditions_col]] <- cond_chr
}

# tic_norm() sets NA to 0 before it normalizes, and zeros cannot be shown in the RLA plot (log of 0). The data should be
# imputed first (metaproviz/processingimpute), so NA or 0 values are reported here.
missing_mask <- is.na(assay_df) | assay_df == 0
n_missing    <- sum(as.matrix(missing_mask))
n_na         <- sum(is.na(as.matrix(assay_df)))
missing_features <- colnames(assay_df)[colSums(as.matrix(missing_mask)) > 0]
# One row per feature with missing values, counted in the input.
missing_df <- data.frame(
  Feature = missing_features,
  `NA values in input` = vapply(missing_features, function(f) sum(is.na(assay_df[[f]])), numeric(1)),
  `0 values in input`  = vapply(missing_features, function(f) sum(assay_df[[f]] == 0, na.rm = TRUE), numeric(1)),
  check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
if (n_missing > 0)
  log_msg("WARN", n_missing, " missing value(s) (NA or 0) in ", length(missing_features), " feature(s): ",
          paste(missing_features, collapse = ", "), ". tic_norm() sets NA to 0",
          if (n_na > 0) paste0(" (", n_na, " NA value(s) are 0 in the output)") else "",
          ". Missing values should be imputed first, see the metaproviz/processingimpute module.")

metadata_info <- c(Conditions = conditions_col)

tic_res <- tryCatch(
  MetaProViz::tic_norm(
    data            = assay_df,
    metadata_sample = sample_info_df,
    metadata_info   = metadata_info,
    tic             = TRUE
  ),
  error = function(e) abort("tic_norm() failed: ", conditionMessage(e))
)

normalized_df <- tic_res[["DF"]][["data_tic"]]
plots         <- tic_res[["Plot"]]; if (is.null(plots)) plots <- list()
if (is.null(normalized_df)) abort("tic_norm() did not return DF\$data_tic.")

# ── BUILD OUTPUT SE ──────────────────────────────────────────────────────────
# normalized_df is samples x features, the assay is features x samples in the order of se.
normalized_assay <- t(as.matrix(normalized_df))[rownames(se), colnames(se), drop = FALSE]
se_normalized <- se
assay(se_normalized, 1) <- normalized_assay

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
if (output_format != "se" && length(assays(se_normalized)) > 1)
  log_msg("WARN", "The input has ", length(assays(se_normalized)), " assays. Only the first one is written to the data matrix.")
write_out(se_normalized, se_out);                   log_msg("INFO", "Written: ", paste(se_out, collapse = ", "))
write.table(removed_df, removed_out, sep = "\\t", quote = FALSE, row.names = FALSE)
# The file has column names without spaces, the report table keeps the longer ones.
missing_file_df <- setNames(missing_df, c("Feature", "NA_values", "Zero_values"))
write.table(missing_file_df, missing_out, sep = "\\t", quote = FALSE, row.names = FALSE)
saveRDS(plots, plots_out)
log_msg("INFO", "Written: ", removed_out, ", ", missing_out, ", ", plots_out)

# ── HTML REPORT ────────────────────────────────────────────────────────────
log_section("Building report")

plot_titles <- c(
  RLA_BeforeticNorm = "RLA plot, before TIC normalization",
  RLA_AfterticNorm  = "RLA plot, after TIC normalization",
  norm_plots        = "RLA before/after, side by side"
)
plot_html <- ""
plot_order <- c(intersect(c("RLA_BeforeticNorm", "RLA_AfterticNorm", "norm_plots"), names(plots)),
                setdiff(names(plots), c("RLA_BeforeticNorm", "RLA_AfterticNorm", "norm_plots")))
for (nm in plot_order) {
  p <- plots[[nm]]
  if (!inherits(p, "ggplot")) next
  title <- if (nm %in% names(plot_titles)) plot_titles[[nm]] else nm
  src <- tryCatch(png_to_b64(p),
                  error = function(e) {
                    log_msg("WARN", "Could not render plot '", nm, "': ",
                            conditionMessage(e)); NA_character_
                  })
  if (!is.na(src)) {
    plot_html <- paste0(plot_html,
      sprintf('<h3>%s</h3>\\n<img src="%s" alt="%s" />\\n', esc(title), src, esc(title)))
  }
}

# Names for the report: at most 25, and a pointer to the log file, which has all of them, only if some are left out.
names_html <- function(x, max_n = 25) paste0(
  paste0("<code>", esc(head(x, max_n)), "</code>", collapse = " "),
  if (length(x) > max_n) paste0(" and ", length(x) - max_n, " more (see <code>", esc(log_out), "</code>)") else "")
graybox <- function(...) paste0('<div class="graybox">', paste0(...), "</div>")
pl <- function(n, one, many) if (n == 1) one else many
sel <- function(param, value, default) {
  paste0("<code>", param, " ", esc(value), "</code> ",
         if (identical(as.character(value), as.character(default))) "(default)"
         else paste0("(default: ", esc(default), ")"))
}
# Gray box under the conditions column: no such column, samples without a condition, or only one group of samples.
group_values <- unique(sample_info_df[[conditions_col]])
cond_box_html <- if (!has_conditions) graybox(
  "<p><strong>There is no such column in the sample metadata, and therefore all samples are treated as one group.</strong></p>") else
  if (length(no_condition_ids) > 0 || length(group_values) == 1) graybox(
    if (length(no_condition_ids) == ncol(se)) paste0(
      "<p><strong>No sample has a value in the Conditions column.</strong> All samples are treated as one condition called ",
      "'NA' and share one color in the plots.</p>") else paste0(
      if (length(no_condition_ids) > 0) paste0(
        "<p><strong>", length(no_condition_ids), " sample", pl(length(no_condition_ids), " has", "s have"),
        " no value in the Conditions column:</strong></p><p>", names_html(no_condition_ids), "</p><p>",
        pl(length(no_condition_ids), "This sample is", "These samples are"), " treated as one condition called 'NA'.</p>") else "",
      if (length(group_values) == 1) paste0(
        "<p><strong>All samples are in one condition (<code>", esc(group_values), "</code>).</strong> They share one color in the plots.</p>") else "")) else ""
params_html <- paste0(
  "<ul>",
  "<li><p>", sel("--conditions_col", conditions_col, "Conditions"),
  ": the column of the ", samp_meta_txt, " that holds the conditions. The plots are colored by them.</p>", cond_box_html, "</li>",
  "</ul>")

missing_row <- if (n_missing == 0) "none" else paste0(
  n_missing, " in ", length(missing_features), " feature", if (length(missing_features) == 1) "" else "s")

# HTML table of a data frame, in a box that scrolls when it is long.
build_table_html <- function(df) {
  header <- paste0("<th>", esc(colnames(df)), "</th>", collapse = "")
  rows <- apply(df, 1, function(r)
    paste0("<tr>", paste0("<td>", esc(as.character(r)), "</td>", collapse = ""), "</tr>"))
  paste0('<div class="scroll-table"><table><thead><tr>', header,
         '</tr></thead><tbody>', paste(rows, collapse = ""), '</tbody></table></div>')
}
# Gray box in the missing values section: the samples that were removed because their total ion count is 0.
removed_box_html <- if (length(zero_total_ids) > 0) graybox(
  "<p><strong>", length(zero_total_ids), " sample", pl(length(zero_total_ids), " with a total ion count of 0 was", "s with a total ion count of 0 were"),
  " removed:</strong></p><p>", names_html(zero_total_ids), "</p><p>Every value of ", pl(length(zero_total_ids), "this sample is", "these samples is"),
  " NA or 0, which would turn the normalization into a division by 0.</p>") else ""


missing_section_html <- paste0(
  "<h2>Missing values</h2>\\n<p>",
  "If missing values are left, the normalization counts NA as 0, and values of 0 are left out of the RLA plots. They should ",
  "be imputed first, see the <code>metaproviz/processingimpute</code> module.</p>\\n",
  if (n_missing == 0) paste0("<p><em>None, no missing values were found", if (length(zero_total_ids) > 0) " in the remaining samples" else "", ".</em></p>\\n") else
    paste0("<p>", n_missing, " missing value", if (n_missing == 1) "" else "s", " (NA or 0) in ", length(missing_features),
           " feature", if (length(missing_features) == 1) "" else "s", ". The table below counts the NA and the 0 values in the input, ",
           "before the normalization sets NA to 0. It is saved as <code>", esc(missing_out), "</code>.</p>\\n", build_table_html(missing_df), "\\n"),
  removed_box_html)

run_meta_html <- sprintf(
'<dl>
   <dt>Input</dt><dd>%s</dd>
   <dt>n features</dt><dd>%d</dd>
   <dt>n samples</dt><dd>%s</dd>
   <dt>Missing values (NA or 0)</dt><dd>%s</dd>
 </dl>',
 esc(source_label), nrow(se),
 paste0(ncol(se), " sample", if (ncol(se) == 1) "" else "s",
        if (length(zero_total_ids) > 0) paste0(" after ", length(zero_total_ids), " sample",
                                               if (length(zero_total_ids) == 1) " was" else "s were", " removed") else ""),
 missing_row)

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
<title>Processing: TIC Normalization Report</title>
<style>
  body{font-family:Arial,sans-serif;max-width:980px;margin:2em auto;padding:0 1em;line-height:1.5}
  h1{color:#1F3A5F;border-bottom:2px solid #1F3A5F;padding-bottom:.3em}
  h2{color:#1F3A5F;margin-top:2em}
  h3{color:#345;margin-top:1.4em}
  dl{display:grid;grid-template-columns:max-content 1fr;gap:.3em 1em}
  dt{font-weight:bold;color:#345}
  img{max-width:100%%;border:1px solid #ddd;padding:4px;background:#fff;margin:.4em 0 1.2em}
  code{background:#f2f2f2;padding:.15em .35em;border-radius:3px}
  .graybox{background:#f2f2f2;border-left:4px solid #999;padding:.1em 1em;margin:1em 0;border-radius:3px}
  .graybox code{background:#e0e0e0}
  table{border-collapse:collapse;width:100%%;font-size:.9em}
  th,td{border:1px solid #ccc;padding:.35em .6em;text-align:left}
  th{background:#f0f4f8}
  .scroll-table{max-height:320px;overflow-y:auto;border:1px solid #cbd6e4;border-radius:6px;margin:.4em 0}
  .scroll-table table{border:0;margin:0}
  .scroll-table thead th{position:sticky;top:0;background:#f0f4f8;z-index:1}
</style></head><body>
<h1>Processing: TIC Normalization</h1>
<p><em>Generated by the METAPROVIZ_PROCESSINGNORMALIZE module.</em></p>

<h2>Run summary</h2>%s

<h2>Citations</h2>
<p>
If you use results or plots from this analysis in a publication, please cite:
</p>
<ul>
  <li><strong>MetaProViz:</strong>
  <a href="https://doi.org/10.1038/s44320-026-00231-8">Schmidt et al., Integrated metabolomics data analysis to generate mechanistic hypotheses with MetaProViz, Molecular Systems Biology 2026.</a></li>
  <li><strong>Dependencies:</strong>
  <a href="https://doi.org/10.32614/R.manuals">R Core Team, R: A Language and Environment for Statistical Computing, R Foundation for Statistical Computing</a>, <a href="https://doi.org/10.18129/B9.bioc.SummarizedExperiment">Morgan et al., SummarizedExperiment: A container (S4 class) for matrix-like assays, Bioconductor</a>, <a href="https://doi.org/10.1038/nmeth.3252">Huber et al., Orchestrating high-throughput genomic analysis with Bioconductor, Nature Methods 2015</a>, <a href="https://doi.org/10.1007/978-3-319-24277-4">Wickham, ggplot2: Elegant Graphics for Data Analysis, Springer-Verlag New York 2016</a> and <a href="https://doi.org/10.32614/CRAN.package.gridExtra">Auguie, gridExtra: Miscellaneous Functions for "Grid" Graphics, CRAN</a>.
  The other R packages used are listed with their versions in the session summary at the end of this report.</li>
  <li><strong>Normalization method:</strong>
  <a href="https://doi.org/10.4236/abb.2018.98022">Wulff and Mitchell, A Comparison of Various Normalization Methods for LC/MS Metabolomics Data, Advances in Bioscience and Biotechnology 2018.</a></li>
  <li><strong>This module:</strong> <code>metaproviz/processingnormalize</code>, an nf-core module. It has no publication of its own and is cited through the nf-core publication below. In the methods, name the module.</li>
  <li><strong>nf-core:</strong> <a href="https://doi.org/10.1038/s41587-020-0439-x">Ewels et al., The nf-core framework for community-curated bioinformatics pipelines, Nature Biotechnology 2020.</a></li>
  <li><strong>Nextflow:</strong> <a href="https://doi.org/10.1038/nbt.3820">Di Tommaso et al., Nextflow enables reproducible computational workflows, Nature Biotechnology 2017.</a></li>
</ul>

<h2>Description</h2>
<p>
This is the third module in the processing chain (<code>processingfilter</code>, <code>processingimpute</code>,
<code>processingnormalize</code> and <code>processingoutliers</code>). It continues from the output of
<code>processingimpute</code>. The data should therefore contain no pool or QC samples, and missing values should be
imputed first.
</p>
<p>
This module uses
<a href="https://saezlab.github.io/MetaProViz/reference/tic_norm.html"><code>MetaProViz::tic_norm()</code></a>
for the total ion count (TIC) normalization. It divides each sample by the sum of its detected intensities and multiplies the result by the median of these sums over
all samples, which keeps the values on the scale of the input (Wulff and Mitchell, 2018). A sample with a total ion count of 0 is removed from the output data, because dividing by 0 is not possible. Since
the total signal of a sample also depends on the measurement, for example on the amount that was injected or on the
sensitivity of the instrument, the division by the total corrects for technical differences between the samples.
</p>
<p><strong>The parameters selected for this run:</strong></p>
%s

%s

<h2>Plots</h2>
<p>
The RLA (relative log abundance) plots are boxplots with one box per sample. They show how the values of each sample are
distributed relative to a reference. For a direct visual comparison, look at the plot "RLA before/after, side by side",
which shows the before and the after plot next to each other. If the TIC normalization worked as intended, the "after"
boxplots look tighter and more centered around 0 than the "before" ones, since the RLA values of a well-normalized sample
center on 0. Samples that
still stand out afterwards may deserve a closer look. The boxes are colored by condition. Without a Conditions column,
or if all samples have the same condition, all samples form one group and share one color.
</p>
%s

<h2>Run notes (log)</h2>
<p>A separate file, <code>%s</code>, records every step of this run,
including any warnings. Check it first if a result here looks unexpected.</p>

<h2>Session summary</h2>
<p>Full R session used to generate this report — every package, and its version:</p>
%s

</body></html>',
run_meta_html, params_html, missing_section_html,
if (nzchar(plot_html)) plot_html else "<p><em>No plots returned.</em></p>",
esc(log_out), session_summary_html)

writeLines(html, report_out)
log_msg("INFO", "Written: ", report_out)

# MetaProViz and OmnipathR write their own log files and plot dumps to the working directory. They are no outputs of
# this module and are deleted here, at the very end.
unlink("metaproviz-log", recursive = TRUE)
unlink("omnipathr-log", recursive = TRUE)
unlink(Sys.glob("Rplots*.pdf"))

flush_log()
message(sprintf(
  "\\n✓ processing_normalize complete | %d warning(s) | %d error(s)",
  warn_count, err_count))
