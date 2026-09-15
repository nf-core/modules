#!/usr/bin/env Rscript
# =============================================================================
# internal_standard.R — technical-QC check on the internal-standard spike-ins
#
# Approach: re-use MetaProViz::pool_estimation on the SE subsetted to just the
# IS features, with every sample flagged as "Pool" (because for IS we want the
# CV across all samples — IS should be constant, so a high CV reveals LC-MS
# technical issues, not biology).
#
# Picking the IS features (top precedence first):
#   1. --internal_standards  "valine-d8,hippuric acid-d5,..."   (explicit list)
#   2. --pattern             "regex"                            (user regex)
#   3. Default auto-detect on the common IS naming conventions:
#         -dN suffix         (deuterated, e.g. valine-d8)
#         -13CN suffix       (carbon-13 labelled)
#         ISTD anywhere
#         IS_ prefix / _IS suffix / bare 'IS'
#
# Pool/QC samples are always excluded before the CV is computed when any are
# found (pool detection mirrors the POOL_ESTIMATION module). 
# Pools are homogeneous mixtures whose tight IS values would mask drift seen only 
# in the real samples. If no pools are found, all samples are used.
# =============================================================================

set.seed(42)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(magrittr)
  library(ggplot2); library(SummarizedExperiment); library(S4Vectors)
})

suppressPackageStartupMessages(library(MetaProViz))

# ── ARG PARSING ─────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args); if (is.na(idx) || idx == length(args)) return(default)
  v <- args[idx + 1]; if (startsWith(v, "--")) return(default); v
}

se_path             <- get_arg("--se")
data_matrix_path    <- get_arg("--data_matrix")
feature_matrix_path <- get_arg("--feature_matrix")
sample_matrix_path  <- get_arg("--sample_matrix")
is_list_str    <- get_arg("--internal_standards", "")
user_pattern   <- get_arg("--pattern",             "")
cutoff_cv      <- as.numeric(get_arg("--cutoff_cv", "30"))

# Pool/QC samples are ALWAYS dropped before the IS CV is computed (when any are
# found). Pools are homogeneous mixtures, so their internal-standard values are
# tighter than real samples and would pull the overall CV down, masking drift
# that is only visible in the biological samples. Pool samples are detected with
# the same precedence as the POOL_ESTIMATION module: explicit list >>> metadata
# column+value >>> name pattern. If no pools are found, all samples are used.
pool_samples_str    <- get_arg("--pool_samples",        "")
pool_metadata_col   <- get_arg("--pool_metadata_col",   "")
pool_metadata_value <- get_arg("--pool_metadata_value", "Pool")
pool_pattern        <- get_arg("--pool_pattern",        "pool")

# All output filenames are `<prefix>.<suffix>`, per nf-core naming
# convention (https://nf-co.re/docs/specifications/components/modules/naming-conventions):
# "output file names SHOULD consist of only ${prefix} and the file-format
# suffix" — needed so filenames don't collide across samples when this
# module runs on many samples in one pipeline.
prefix              <- get_arg("--prefix", "internal_standard")
cv_out              <- paste0(prefix, ".cv.tsv")
high_var_out        <- paste0(prefix, ".high_var.txt")
condition_cv_out    <- paste0(prefix, ".condition_cv.tsv")
plots_out           <- paste0(prefix, ".plots.rds")
report_out          <- paste0(prefix, ".report.html")
log_out             <- paste0(prefix, ".log")

# Default regex — case-insensitive PCRE. Covers common IS labels.
DEFAULT_PATTERN <- "-d\\d+|-13C\\d*|ISTD|^IS_|_IS$|^IS$"

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
  message(sep, "\n  ", t, "\n", sep)
}
flush_log <- function(p = log_out) {
  writeLines(c(ll, "",
    sprintf("SUMMARY: %d warning(s), %d error(s)", wc, ec)), p)
}
abort <- function(...) { log_msg("ERROR", ...); flush_log(); stop(paste(...), call. = FALSE) }
esc <- function(x) gsub(">", "&gt;", gsub("<", "&lt;",
                    gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE), fixed = TRUE)

placeholder_outputs <- function(reason_html) {
  # Always emit all declared outputs so Nextflow doesn't fail the process.
  # (condition_cv and versions.yml are both non-optional in main.nf; a
  # missing declared output fails the Nextflow process even when this R
  # script itself exits 0 — versions.yml is written separately by main.nf's
  # own script: block after this script returns, so it's covered, but
  # condition_cv has to be written here explicitly.)
  writeLines(character(0),                      high_var_out)
  write.table(data.frame(Metabolite=character(0), CV=numeric(0),
                         HighVar=logical(0)),
              cv_out, sep="\t", quote=FALSE,
              row.names=FALSE, na="NA")
  # Same header-only shape used when real condition metadata is absent
  # (see condition_cv_df's `else` branch further down).
  write.table(data.frame(Standard = character(0)),
              condition_cv_out, sep="\t", quote=FALSE,
              row.names=FALSE, na="NA")
  saveRDS(list(), plots_out)
  writeLines(sprintf(
    '<!DOCTYPE html><html><head><meta charset="utf-8"><title>IS Report</title>
<style>body{font-family:Arial,sans-serif;max-width:980px;margin:2em auto;
padding:0 1em;line-height:1.5}h1{color:#1F3A5F;border-bottom:2px solid #1F3A5F;
padding-bottom:.3em}p{margin:.5em 0}.note{background:#fff8e1;border:1px solid
#ffe082;border-radius:6px;padding:.7em 1em;color:#614e00}
</style></head><body>
<h1>Internal Standard Report</h1>
<p><em>Generated by the nf_metabolism pipeline (METAPROVIZ_INTERNALSTANDARD module).</em></p>
<div class="note">%s</div></body></html>', reason_html),
    report_out)
}

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
    sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
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

# ── ALWAYS EXCLUDE POOL SAMPLES (when any are found) ─────────────────────────
pool_n        <- 0L
pool_strategy <- "n/a"
{
  log_section("Excluding pool samples")
  sinfo_full <- as.data.frame(colData(se), check.names = FALSE)
  pool_ids   <- character(0)

  # 1) Explicit list
  if (nzchar(pool_samples_str)) {
    asked <- trimws(strsplit(pool_samples_str, ",")[[1]]); asked <- asked[nzchar(asked)]
    pool_ids <- intersect(asked, colnames(se))
    pool_strategy <- "explicit list (--pool_samples)"
  }
  # 2) Metadata column + value
  if (length(pool_ids) == 0 && nzchar(pool_metadata_col) &&
      pool_metadata_col %in% colnames(sinfo_full) && nzchar(pool_metadata_value)) {
    pool_ids <- rownames(sinfo_full)[
      as.character(sinfo_full[[pool_metadata_col]]) == pool_metadata_value]
    pool_strategy <- sprintf("metadata column '%s' == '%s'",
                             pool_metadata_col, pool_metadata_value)
  }
  # 3) Name pattern (default fallback)
  if (length(pool_ids) == 0 && nzchar(pool_pattern)) {
    pool_ids <- colnames(se)[grepl(pool_pattern, colnames(se), ignore.case = TRUE)]
    pool_strategy <- sprintf("pattern '%s' on sample names", pool_pattern)
  }

  pool_n <- length(pool_ids)
  if (pool_n > 0 && pool_n < ncol(se)) {
    log_msg("INFO", "Excluding ", pool_n, " pool sample(s) via ", pool_strategy,
            "; keeping ", ncol(se) - pool_n, " for the IS CV.")
    se <- se[, setdiff(colnames(se), pool_ids), drop = FALSE]
  } else if (pool_n >= ncol(se)) {
    log_msg("WARN", "Pool detection matched ALL samples (", pool_n,
            ") — refusing to drop everything; using all samples instead.")
    pool_n <- 0L; pool_strategy <- "all samples matched, exclusion skipped"
  } else {
    log_msg("INFO", "No pool samples found — using all samples for the IS CV.")
    pool_strategy <- "no pool samples found"
  }
}

# ── IDENTIFY IS FEATURES ────────────────────────────────────────────────────
log_section("Identifying internal-standard features")

feat_names <- rownames(se)
is_feats <- character(0)
strategy <- "none"

# 1) Explicit list
if (nzchar(is_list_str)) {
  asked <- trimws(strsplit(is_list_str, ",")[[1]])
  asked <- asked[nzchar(asked)]
  is_feats <- intersect(asked, feat_names)
  missing  <- setdiff(asked, feat_names)
  if (length(missing) > 0)
    log_msg("WARN", "Requested IS feature(s) not found in SE: ",
            paste(missing, collapse = ", "))
  strategy <- "Explicit list provided by --internal_standards"
}

# 2) User pattern
if (length(is_feats) == 0 && nzchar(user_pattern)) {
  hit <- grepl(user_pattern, feat_names, ignore.case = TRUE, perl = TRUE)
  is_feats <- feat_names[hit]
  strategy <- paste0("User-supplied pattern (--internal_standard_pattern)")
}

# 3) Default pattern
if (length(is_feats) == 0) {
  hit <- grepl(DEFAULT_PATTERN, feat_names, ignore.case = TRUE, perl = TRUE)
  is_feats <- feat_names[hit]
  strategy <- "Auto-detected from feature names (deuterated, 13C, IS / ISTD naming)"
}

if (length(is_feats) == 0) {
  log_msg("WARN", "No internal-standard features found by any strategy.")
  placeholder_outputs(paste0(
    "No internal-standard features were detected in this SE — neither in the ",
    "explicit list (if any), nor by the user pattern (if any), nor by the ",
    "default pattern (deuterated -dN, 13C, ISTD, IS_, _IS). The step is ",
    "harmless to skip; pass --internal_standards \"...\" or ",
    "--internal_standard_pattern \"...\" to point at the right names if your ",
    "data does include IS."))
  flush_log(log_out)
  message("\n✓ internal_standard complete (no IS features) | ",
          wc, " warning(s)")
  quit(save = "no", status = 0)
}

log_msg("INFO", "Detection strategy: ", strategy)
log_msg("INFO", "IS features (", length(is_feats), "): ",
        paste(head(is_feats, 12), collapse = ", "),
        if (length(is_feats) > 12) " ..." else "")

# ── SUBSET SE TO IS FEATURES + STAMP COND=POOL ──────────────────────────────
log_section("Preparing data for MetaProViz::pool_estimation")

se_is <- se[is_feats, , drop = FALSE]
assay_df       <- t(assay(se_is, 1)) %>% as.data.frame(check.names = FALSE)
sample_info_df <- as.data.frame(colData(se_is), check.names = FALSE)

# For pool_estimation: all samples are "pools" (we want CV across everything).
# Preserve any original Conditions column.
if ("Conditions" %in% colnames(sample_info_df)) {
  sample_info_df$Conditions_original <- sample_info_df$Conditions
}
sample_info_df$Conditions <- "Pool"

# ── RUN pool_estimation ─────────────────────────────────────────────────────
log_section("Running MetaProViz::pool_estimation (limited to IS features)")

pe <- tryCatch(
  MetaProViz::pool_estimation(
    data            = assay_df,
    metadata_sample = sample_info_df,
    metadata_info   = c(PoolSamples = "Pool", Conditions = "Conditions"),
    cutoff_cv       = cutoff_cv,
    print_plot      = FALSE,
    save_plot       = NULL
  ),
  error = function(e) abort("pool_estimation() failed on IS subset: ",
                            conditionMessage(e))
)

cv_df <- pe[["DF"]][["CV"]]
if (is.null(cv_df)) abort("pool_estimation did not return DF$CV.")
plots <- pe[["Plot"]]; if (is.null(plots)) plots <- list()

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
superplot_meta$Conditions <- NULL
superplot_meta$AllSamples <- "IS"
has_real_conditions <- "Conditions_original" %in% colnames(superplot_meta)
superplot_metadata_info <- if (has_real_conditions) {
  c(Conditions = "AllSamples", Superplot = "Conditions_original")
} else {
  c(Conditions = "AllSamples")
}

superplot_res <- tryCatch(
  MetaProViz::viz_superplot(
    data            = assay_df,
    metadata_sample = superplot_meta,
    metadata_info   = superplot_metadata_info,
    plot_type       = "Box",
    print_plot      = FALSE,
    # save_plot = NULL hits a bug in viz_superplot() (it still tries to use
    # an internal `folder` variable that only gets set when save_plot is
    # non-NULL). Workaround: give it a throwaway path and let it keep its
    # default save_plot = "svg" — we only use the returned Plot object
    # anyway, the files it writes here are never read.
    path            = tempdir()
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
  long_df <- assay_df %>%
    tibble::rownames_to_column("Sample") %>%
    tidyr::pivot_longer(-Sample, names_to = "Standard", values_to = "Value") %>%
    dplyr::mutate(Condition = sample_info_df[Sample, "Conditions_original"]) %>%
    dplyr::filter(!is.na(Value))

  long_df %>%
    dplyr::group_by(Standard, Condition) %>%
    dplyr::summarise(
      CV = if (dplyr::n() >= 2) 100 * sd(Value) / mean(Value) else NA_real_,
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(id_cols = Standard, names_from = Condition, values_from = CV)
} else {
  data.frame(Standard = character(0))
}

# ── WRITE OUTPUTS ──────────────────────────────────────────────────────────
log_section("Writing outputs")

high_var <- cv_df %>% filter(HighVar == TRUE) %>% pull(Metabolite)
log_msg("INFO", "IS with CV > ", cutoff_cv, "%: ",
        length(high_var), " / ", nrow(cv_df))

write.table(cv_df, cv_out,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
writeLines(if (length(high_var) > 0) high_var else character(0),
           high_var_out)
write.table(condition_cv_df, condition_cv_out,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
saveRDS(list(PoolEstimationStyle = plots, Superplot = superplot_plots),
        plots_out)
log_msg("INFO", "Written: ", cv_out, ", ", high_var_out, ", ",
        condition_cv_out, ", ", plots_out)

# ── HTML REPORT ────────────────────────────────────────────────────────────
log_section("Building report")

png_to_b64 <- function(p, w = 8, h = 5) {
  tmp <- tempfile(fileext = ".png")
  ggplot2::ggsave(tmp, plot = p, width = w, height = h, dpi = 100,
                  units = "in", bg = "white")
  raw <- readBin(tmp, "raw", n = file.info(tmp)$size); unlink(tmp)
  paste0("data:image/png;base64,", base64enc::base64encode(raw))
}

# One box plot per standard. Name + CV% go on the plot's own x-axis label.
superplot_html <- ""
if (length(superplot_plots) > 0) {
  cv_lookup <- setNames(cv_df$CV, cv_df$Metabolite)
  for (nm in names(superplot_plots)) {
    p <- superplot_plots[[nm]]
    if (!inherits(p, "ggplot")) next
    cv_val <- cv_lookup[[nm]]
    cv_str <- if (!is.null(cv_val) && !is.na(cv_val)) format(round(cv_val, 1)) else NA
    axis_label <- if (!is.na(cv_str)) sprintf("%s (%s%%)", nm, cv_str) else nm
    plot_title <- if (!is.na(cv_str)) {
      sprintf("Internal Standard: %s (%s%%)", nm, cv_str)
    } else {
      sprintf("Internal Standard: %s", nm)
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
        sprintf('<img src="%s" alt="%s" />\n', src, esc(axis_label))
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
  cond_cols <- setdiff(colnames(condition_cv_df), "Standard")
  cond_header <- paste0("<th>", esc(c("Standard", cond_cols)), "</th>", collapse = "")
  cond_rows <- apply(condition_cv_df, 1, function(r) {
    cells <- vapply(r[-1], function(v) {
      if (is.na(v)) "—" else sprintf("%s%%", format(round(as.numeric(v), 1)))
    }, character(1))
    paste0("<tr><td>", esc(r[["Standard"]]), "</td>",
           paste0("<td>", cells, "</td>", collapse = ""), "</tr>")
  })

  condition_cv_section_html <- sprintf(
    '<h2>CV per standard, per condition</h2>
<p>
Expectation: low CV across every condition, since internal standards should behave
technically the same regardless of biological condition. High CV in all conditions points to
a general, run-wide issue (e.g. instrument drift). High CV in only some conditions is more
concerning: it points to a condition-specific issue (e.g. a batch effect from how those samples
were handled), since a real internal standard should not care which condition a sample belongs
to. The table is also saved as <code>%s</code>.
</p>
<div class="scroll-table"><table><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>',
    esc(condition_cv_out), cond_header, paste(cond_rows, collapse = "")
  )
}

high_var_html <- if (length(high_var) > 0) {
  paste0("<ul>", paste0("<li>", esc(high_var), "</li>", collapse = ""),
         "</ul>")
} else "<p><em>None — every internal standard is below the CV cutoff.</em></p>"

# Render the actual matched IS feature names as a readable comma-separated
# list (not a regex). Lists the names of the standards, not pattern metacharacters.
is_names_html <- if (length(is_feats) > 0) {
  paste(esc(is_feats), collapse = ", ")
} else "—"

# CV column: round + add "%" to match the per-condition CV table's style,
# instead of full-precision decimals with no unit.
cv_rows <- apply(cv_df, 1, function(r) {
  vals <- vapply(names(r), function(col) {
    if (col == "CV") {
      v <- suppressWarnings(as.numeric(r[[col]]))
      if (!is.na(v)) sprintf("%s%%", format(round(v, 1))) else as.character(r[[col]])
    } else {
      as.character(r[[col]])
    }
  }, character(1))
  paste0("<tr>", paste0("<td>", esc(vals), "</td>", collapse = ""), "</tr>")
})
cv_table <- paste0(
  '<div class="scroll-table"><table><thead><tr>',
  paste0("<th>", esc(colnames(cv_df)), "</th>", collapse = ""),
  '</tr></thead><tbody>',
  paste(cv_rows, collapse = ""),
  '</tbody></table></div>')

# Pool samples are always excluded when found; state what happened.
pool_mode_html <- if (pool_n > 0) {
  sprintf("Excluded — %d pool sample(s) dropped via %s; CV computed on the remaining %d sample(s).",
          pool_n, esc(pool_strategy), ncol(se))
} else {
  sprintf("No pool samples found (%s); CV computed across all %d sample(s).",
          esc(pool_strategy), ncol(se))
}

run_meta_html <- sprintf(
'<dl>
   <dt>Input source</dt><dd>%s</dd>
   <dt>n features (input → IS)</dt><dd>%d → %d</dd>
   <dt>n samples (input → used)</dt><dd>%d → %d</dd>
   <dt>Pool samples in CV</dt><dd>%s</dd>
   <dt>Detection strategy</dt><dd>%s</dd>
   <dt>Internal standards used</dt><dd>%s</dd>
   <dt>CV cutoff</dt><dd>%s%%</dd>
   <dt>IS above cutoff</dt><dd>%d / %d</dd>
 </dl>',
 esc(source_label),
 nrow(se), length(is_feats), n_samples_input, ncol(se),
 pool_mode_html,
 esc(strategy), is_names_html, format(cutoff_cv),
 length(high_var), nrow(cv_df))

html <- sprintf('<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Internal Standard Report</title>
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
  li{margin:.25em 0}
  .scroll-table{max-height:320px;overflow-y:auto;border:1px solid #cbd6e4;border-radius:6px}
  .scroll-table table{border:0;margin:0}
  .scroll-table thead th{position:sticky;top:0;background:#f0f4f8;z-index:1}
</style></head><body>
<h1>Internal Standard Report</h1>
<p><em>Generated by the nf_metabolism pipeline (METAPROVIZ_INTERNALSTANDARD module).</em></p>

<h2>Run summary</h2>%s

<h2>Citations</h2>
<p>
If you use results or plots from this analysis in a publication, please cite:
</p>
<ul>
  <li><strong>MetaProViz:</strong> Please cite:
  <a href="https://doi.org/10.1038/s44320-026-00231-8">Schmidt et al., Integrated metabolomics data analysis to generate mechanistic hypotheses with MetaProViz, Molecular Systems Biology 2026.</a>.</li>
  <li><strong>Dependencies:</strong> None for this module.</li>
</ul>

<h2>Description</h2>
<p>
This report evaluates internal-standard (spike-in) features to assess LC-MS technical noise,
independent of biology. Internal standards are added in a fixed, known amount to every sample
before measurement, so unlike real metabolites, their values are expected to stay constant
across the run.
</p>
<p>
This module computes each internal standard\'s coefficient of variation (CV) across all real
samples. Pool and QC samples are excluded, since their already-pooled nature would understate
true injection-to-injection variability. A high CV indicates instrument drift or injection
issues rather than biological variation.
</p>
<p>
CV is calculated as:
</p>
<p><code>CV(%%) = 100 × sd(x1, x2, ..., xn) / mean(x1, x2, ..., xn)</code></p>
<p>
This module reuses <code>MetaProViz::pool_estimation()</code> internally, applied only to the
internal-standard features, with every sample temporarily marked <code>Pool</code> so that the
same pool-CV computation logic applies to this different feature subset and framing.
</p>
<p>
If condition metadata is available, the same CV calculation is additionally broken down per
condition (see the "CV per standard, per condition" table further down), and the box plots
below also color individual dots by condition.
</p>

<h2>Plots</h2>
<p>
Each internal standard gets one box plot below, generated with MetaProViz\'s
<code>viz_superplot()</code> function, showing its values across all samples used for the CV
calculation.
</p>
<p>
What to expect: dots of different colors should be interspersed throughout the box, not
clustered by color, since the internal standard is expected to behave the same regardless of
condition. If dots of one color cluster together at their own distinct position, that points to
a systematic shift for that specific condition (e.g. a batch effect). Without condition
metadata, all dots are shown in black instead.
</p>
%s
%s

<h2>CV per standard, overall</h2>
<p>
Every internal standard\'s CV, computed across all real samples together. The table is also
saved as <code>%s</code>.
</p>
%s

<h2>Internal standards above the CV cutoff (%s%%)</h2>
%s

<h2>Run notes (log)</h2>
<p>A separate file, <code>%s</code>, records every step of this run,
including any warnings. Check it first if a result here looks unexpected.</p>

</body></html>',
run_meta_html,
if (nzchar(superplot_html)) superplot_html else "<p><em>No box plots returned.</em></p>",
condition_cv_section_html,
esc(cv_out),
cv_table, format(cutoff_cv), high_var_html,
esc(log_out))

writeLines(html, report_out)
log_msg("INFO", "Written: ", report_out)

flush_log(log_out)
message(sprintf(
  "\n✓ internal_standard complete — %d / %d IS flagged HighVar | %d warning(s)",
  length(high_var), nrow(cv_df), wc))
