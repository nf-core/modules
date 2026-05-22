#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(readr)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(purrr)
    library(ggplot2)
    library(DOTSeq)
    library(SummarizedExperiment)
})

################################################################################
## Parse parameters                                                           ##
################################################################################

option_list <- list(
    make_option("--output_prefix",       type = "character", default = NULL),
    make_option("--count_file",          type = "character", default = NULL),
    make_option("--sample_file",         type = "character", default = NULL),
    make_option("--flattened_gtf",       type = "character", default = NULL),
    make_option("--flattened_bed",       type = "character", default = NULL),
    make_option("--contrast_variable",   type = "character", default = NULL),
    make_option("--reference_level",     type = "character", default = NULL),
    make_option("--target_level",        type = "character", default = NULL),
    make_option("--sample_id_col",       type = "character", default = "run"),
    make_option("--strategy_col",        type = "character", default = "strategy"),
    make_option("--replicate_col",       type = "character", default = "replicate"),
    make_option("--sample_name_regex",   type = "character", default = NULL,
                help = "Regex applied to count-table column names; the first capture group is kept (matches DOTSeq's vignette pattern)."),
    make_option("--modules",             type = "character", default = "DOU,DTE",
                help = "Which DOTSeq modules to run [default: %default]"),
    make_option("--min_count",           type = "integer",   default = 1L),
    make_option("--stringent",           type = "character", default = "TRUE",
                help = "TRUE / FALSE / NULL [default: %default]"),
    make_option("--dispersion_modeling", type = "character", default = "auto"),
    make_option("--nullweight",          type = "double",    default = 500),
    make_option("--contrasts_method",    type = "character", default = "revpairwise"),
    make_option("--generate_plots",      type = "logical",   default = TRUE),
    make_option("--alpha",               type = "double",    default = 0.05,
                help = "Padj cut-off for the DTE p-value distribution plot"),
    make_option("--top_hits",            type = "integer",   default = 25L),
    make_option("--cores",               type = "integer",   default = 1L)
)

# Defaults wired in by the Nextflow template; task.ext.args (if any) layers
# on top so users can override anything via `--key value`.
nf_defaults <- c(
    paste0("--output_prefix=",     ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')),
    paste0("--count_file=",        '$counts'),
    paste0("--sample_file=",       '$samplesheet'),
    paste0("--flattened_gtf=",     '$flattened_gtf'),
    paste0("--flattened_bed=",     '$flattened_bed'),
    paste0("--contrast_variable=", '$contrast_variable'),
    paste0("--reference_level=",   '$reference'),
    paste0("--target_level=",      '$target'),
    paste0("--cores=",             '$task.cpus')
)

ext_args_raw <- '$task.ext.args'
ext_argv <- if (identical(ext_args_raw, "null") || !nzchar(trimws(ext_args_raw))) {
    character(0)
} else {
    strsplit(ext_args_raw, "\\\\s+", perl = TRUE)[[1]] |> (\\(x) x[nzchar(x)])()
}

opt <- parse_args(OptionParser(option_list = option_list), args = c(nf_defaults, ext_argv))

# DOTSeq accepts TRUE, FALSE, or NULL for `stringent` (three filter modes);
# optparse won't natively parse a tri-state into a logical, so we round-trip.
stringent_val <- switch(toupper(opt\$stringent),
    "TRUE"  = TRUE,
    "FALSE" = FALSE,
    "NULL"  = NULL,
    stop("`--stringent` must be one of TRUE, FALSE, NULL")
)
modules <- trimws(strsplit(opt\$modules, ",")[[1]])

walk(c("count_file", "sample_file", "flattened_gtf", "flattened_bed"), \\(x) {
    if (!file.exists(opt[[x]])) stop("Missing input file: ", x, " = ", opt[[x]])
})

prefix <- opt\$output_prefix

################################################################################
## Read inputs and normalise the sample sheet                                 ##
################################################################################

# featureCounts emits a `# Program:...` banner as the first line; `comment`
# strips it. Returns a tibble; DOTSeqDataSetsFromFeatureCounts() wants a
# vanilla data.frame so coerce.
cnt <- read_tsv(opt\$count_file, comment = "#", show_col_types = FALSE,
                progress = FALSE) |> as.data.frame()

# featureCounts column names often carry the full BAM path; the vignette uses
# `gsub(".*(SRR[0-9]+).*", "\\1", names(cnt))` to keep just the run accession.
# Expose the same regex via the CLI so users can adapt to their own naming.
if (!is.null(opt\$sample_name_regex) && nzchar(opt\$sample_name_regex)) {
    names(cnt) <- gsub(opt\$sample_name_regex, "\\\\1", names(cnt))
}

cond <- read_csv(opt\$sample_file, show_col_types = FALSE, progress = FALSE) |>
    as.data.frame() |>
    rename_with(\\(nm) tolower(trimws(nm)))

# DOTSeq's parse_condition_table() insists on columns named exactly
# `run, strategy, replicate, condition` (lower-case). Allow the user to point
# at differently-named columns via task.ext.args and rename in-place.
col_map <- c(
    run       = tolower(opt\$sample_id_col),
    strategy  = tolower(opt\$strategy_col),
    replicate = tolower(opt\$replicate_col),
    condition = tolower(opt\$contrast_variable)
)
missing_cols <- col_map[!col_map %in% names(cond)]
if (length(missing_cols) > 0) {
    stop(sprintf("Sample sheet missing column(s): %s. Have: %s",
                 paste(missing_cols, collapse = ", "),
                 paste(names(cond),  collapse = ", ")))
}
for (req in names(col_map)) {
    src <- col_map[[req]]
    if (src != req) names(cond)[names(cond) == src] <- req
}

# Subset to the two contrast levels and put `reference` first so it becomes
# the implicit baseline in the DESeq2 / glmmTMB design.
cond <- cond |>
    filter(.data\$condition %in% c(opt\$reference_level, opt\$target_level)) |>
    mutate(
        condition = factor(.data\$condition, levels = c(opt\$reference_level, opt\$target_level)),
        strategy  = factor(.data\$strategy)
    )

if (nrow(cond) == 0) stop("No samples remain after filtering on the contrast levels.")

################################################################################
## DOTSeq: DOU + DTE                                                          ##
################################################################################

d <- DOTSeqDataSetsFromFeatureCounts(
    count_table     = cnt,
    condition_table = cond,
    flattened_gtf   = opt\$flattened_gtf,
    flattened_bed   = opt\$flattened_bed,
    min_count       = opt\$min_count,
    stringent       = stringent_val,
    baseline        = opt\$reference_level,
    verbose         = FALSE
)

d <- DOTSeq(
    datasets            = d,
    modules             = modules,
    target              = opt\$target_level,
    baseline            = opt\$reference_level,
    min_count           = opt\$min_count,
    stringent           = stringent_val,
    dispersion_modeling = opt\$dispersion_modeling,
    nullweight          = opt\$nullweight,
    contrasts_method    = opt\$contrasts_method,
    parallel            = list(n = opt\$cores, autopar = TRUE),
    verbose             = FALSE
)

# testDOU() and the DTE wrap-up in DOTSeq both lift rownames into an `orf_id`
# column and clear the rownames before returning, so we just coerce to tibble.
get_contrasts_df <- function(x, type) {
    res <- tryCatch(getContrasts(x, type = type), error = \\(e) NULL)
    if (is.null(res)) NULL else as_tibble(as.data.frame(res))
}

dou_interaction <- get_contrasts_df(getDOU(d), "interaction")
dou_strategy    <- get_contrasts_df(getDOU(d), "strategy")
dte_interaction <- get_contrasts_df(getDTE(d), "interaction")
dte_strategy    <- get_contrasts_df(getDTE(d), "strategy")

################################################################################
## Write result tables                                                        ##
##                                                                            ##
## DTE interaction is written as `translation.dotseq.results.tsv` because it  ##
## is the per-ORF differential translation efficiency contrast.               ##
################################################################################

# Always emit the file (even empty) so downstream Nextflow channels stay
# consistent across runs with different significance counts.
empty_safe <- function(df) if (is.null(df)) tibble() else df

write_tsv(empty_safe(dte_interaction), paste0(prefix, ".translation.dotseq.results.tsv"))
write_tsv(empty_safe(dou_interaction), paste0(prefix, ".dou.dotseq.results.tsv"))

if (!is.null(dou_strategy) && nrow(dou_strategy) > 0) {
    write_tsv(dou_strategy, paste0(prefix, ".dou_strategy.dotseq.results.tsv"))
}
if (!is.null(dte_strategy) && nrow(dte_strategy) > 0) {
    write_tsv(dte_strategy, paste0(prefix, ".dte_strategy.dotseq.results.tsv"))
}

saveRDS(d, file = paste0(prefix, ".DOTSeqDataSets.rds"))

################################################################################
## Plots                                                                      ##
##                                                                            ##
## Volcano / composite / venn / heatmap come from DOTSeq's native plotDOT().  ##
## The p-value distribution plot is a plain ggplot on top of the package's   ##
## own DTE padj column.                                                       ##
################################################################################

if (opt\$generate_plots) {

    if (!is.null(dte_interaction) && nrow(dte_interaction) > 0) {
        pdist <- dte_interaction |> filter(!is.na(padj))
        if (nrow(pdist) > 0) {
            pdist_plot <- ggplot(pdist, aes(x = padj)) +
                geom_histogram(bins = 40, fill = "#3498db", colour = "white", alpha = 0.85) +
                geom_vline(xintercept = opt\$alpha, linetype = "dashed", colour = "#e74c3c") +
                labs(
                    x = "Adjusted p-value (DTE interaction)",
                    y = "Count",
                    title = sprintf("DTE p-value distribution (alpha = %s)", opt\$alpha)
                ) +
                theme_bw(base_size = 13)
            ggsave(paste0(prefix, ".interaction_p_distribution.png"),
                   plot = pdist_plot, width = 8, height = 6, dpi = 100)
        }
    }

    # plotDOT defaults to `force_new_device = TRUE` which unconditionally
    # resets the active graphics device, killing the png() we just opened.
    # Disable that so the PNG actually captures the plot.
    safe_plot_dot <- function(plot_type, fname, results_df = NULL, data = NULL,
                              annotation_params = list()) {
        if (is.null(results_df) || nrow(results_df) == 0) return(invisible(NULL))
        tryCatch({
            png(fname, width = 900, height = 800, res = 110)
            plotDOT(
                plot_type         = plot_type,
                results           = results_df,
                data              = data,
                id_mapping        = FALSE,
                plot_params       = list(top_hits = opt\$top_hits),
                annotation_params = annotation_params,
                force_new_device  = FALSE
            )
            dev.off()
        }, error = \\(e) {
            while (length(dev.list()) > 0) dev.off()
            message(sprintf("plotDOT(%s) failed: %s", plot_type, conditionMessage(e)))
        })
    }

    # plotDOT wants both DOU and DTE columns on a single results frame keyed by orf_id.
    plotdot_df <- if (!is.null(dou_interaction) && !is.null(dte_interaction)) {
        dou_interaction |> inner_join(dte_interaction, by = "orf_id", suffix = c("_dou", "_dte"))
    } else NULL

    safe_plot_dot("volcano",   paste0(prefix, ".volcano.png"),   plotdot_df, getDOU(d))
    safe_plot_dot("composite", paste0(prefix, ".composite.png"), plotdot_df, getDOU(d))
    safe_plot_dot("venn",      paste0(prefix, ".venn.png"),      plotdot_df)

    # The heatmap pairs mORFs with a chosen short-ORF class within each gene;
    # try uORF first (the package default) and fall back to dORF if no
    # significant gene has both. tryCatch in safe_plot_dot makes either a
    # no-op if the data don't support it.
    safe_plot_dot("heatmap", paste0(prefix, ".heatmap.png"),
                  plotdot_df, getDOU(d), list(sorf_type = "uORF"))
    if (!file.exists(paste0(prefix, ".heatmap.png"))) {
        safe_plot_dot("heatmap", paste0(prefix, ".heatmap.png"),
                      plotdot_df, getDOU(d), list(sorf_type = "dORF"))
    }
}

################################################################################
## Session info + versions                                                    ##
################################################################################

sink(paste0(prefix, ".R_sessionInfo.log"))
print(sessionInfo())
sink()

writeLines(
    c(
        '"${task.process}":',
        paste0("    bioconductor-dotseq: ", packageVersion("DOTSeq")),
        paste0("    r-optparse: ",          packageVersion("optparse")),
        paste0("    r-readr: ",             packageVersion("readr")),
        paste0("    r-dplyr: ",             packageVersion("dplyr"))
    ),
    "versions.yml"
)
