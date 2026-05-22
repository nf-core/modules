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

is_set <- function(x) !is.null(x) && nzchar(trimws(x))

# DOTSeq's `stringent` is tri-state TRUE / FALSE / NULL; normalise via switch.
# Without the upfront is_set() guard, switch(toupper(NULL), ...) silently
# returns NULL and bypasses the `stop()` fallback.
if (!is_set(opt\$stringent)) stop("--stringent must be TRUE / FALSE / NULL.")
stringent_val <- switch(toupper(opt\$stringent),
    "TRUE"  = TRUE,
    "FALSE" = FALSE,
    "NULL"  = NULL,
    stop("--stringent must be TRUE / FALSE / NULL (got: ", opt\$stringent, ")")
)
modules <- trimws(strsplit(opt\$modules, ",")[[1]])

# An empty output_prefix would produce dotfiles that the emit globs miss.
walk(c("contrast_variable", "reference_level", "target_level", "output_prefix"),
     \\(x) if (!is_set(opt[[x]])) stop("Missing required option: --", x))

prefix <- opt\$output_prefix

################################################################################
## Read inputs and normalise the sample sheet                                 ##
################################################################################

# `comment = "#"` strips featureCounts' `# Program:...` banner.
cnt <- read_tsv(opt\$count_file, comment = "#", show_col_types = FALSE,
                progress = FALSE) |> as.data.frame()

# featureCounts column names often carry the full BAM path; the vignette
# strips them to the run accession via `gsub(".*(SRR[0-9]+).*", "\\1", ...)`.
if (is_set(opt\$sample_name_regex)) {
    names(cnt) <- gsub(opt\$sample_name_regex, "\\\\1", names(cnt))
}

# meta.yml documents CSV or TSV; pick the delim by file extension.
sample_ext <- tolower(tools::file_ext(sub("\\\\.gz\$", "", basename(opt\$sample_file))))
cond <- read_delim(opt\$sample_file, delim = if (sample_ext == "csv") "," else "\t",
                   show_col_types = FALSE, progress = FALSE) |>
    as.data.frame() |>
    rename_with(\\(nm) tolower(trimws(nm)))

# DOTSeq insists on lower-case `run, strategy, replicate, condition` columns;
# rename from user-specified column names if necessary, refusing collisions.
col_map <- c(
    run       = tolower(opt\$sample_id_col),
    strategy  = tolower(opt\$strategy_col),
    replicate = tolower(opt\$replicate_col),
    condition = tolower(opt\$contrast_variable)
)
missing_cols <- setdiff(col_map, names(cond))
if (length(missing_cols) > 0) {
    stop("Sample sheet missing column(s): ", paste(missing_cols, collapse = ", "))
}
nm <- names(cond)
for (req in names(col_map)) {
    src <- col_map[[req]]
    if (src == req) next
    if (req %in% nm) stop("Cannot rename '", src, "' to '", req, "': '", req, "' column already present.")
    nm[nm == src] <- req
}
names(cond) <- nm

if (anyDuplicated(cond\$run)) {
    cond <- distinct(cond, run, .keep_all = TRUE)
}

if (!opt\$reference_level %in% cond\$condition) {
    stop("--reference_level '", opt\$reference_level, "' not in `condition`. Have: ",
         paste(unique(cond\$condition), collapse = ", "))
}
if (!opt\$target_level %in% cond\$condition) {
    stop("--target_level '", opt\$target_level, "' not in `condition`. Have: ",
         paste(unique(cond\$condition), collapse = ", "))
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

# DOTSeq's design `~ condition * strategy` requires both Ribo and RNA samples.
if (nlevels(droplevels(cond\$strategy)) < 2) {
    stop(sprintf("Strategy column must contain both Ribo and RNA after filtering; got: %s",
                 paste(unique(as.character(cond\$strategy)), collapse = ", ")))
}

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
# column and clear the rownames, so we just coerce to tibble. Interaction
# contrasts are the module's headline output: let real errors propagate
# rather than catching them and writing an empty TSV that looks like a
# successful "no significant ORFs". Strategy contrasts can legitimately be
# absent, so tryCatch is fine there.
dou_d <- getDOU(d)
dte_d <- getDTE(d)
contrasts_tibble <- function(res) if (is.null(res)) NULL else as_tibble(as.data.frame(res))
try_contrasts <- function(x, type) tryCatch(contrasts_tibble(getContrasts(x, type = type)), error = \\(e) NULL)

dou_interaction <- contrasts_tibble(getContrasts(dou_d, "interaction"))
dte_interaction <- contrasts_tibble(getContrasts(dte_d, "interaction"))
dou_strategy    <- try_contrasts(dou_d, "strategy")
dte_strategy    <- try_contrasts(dte_d, "strategy")

################################################################################
## Write result tables                                                        ##
##                                                                            ##
## DTE interaction is written as `translation.dotseq.results.tsv` because it  ##
## is the per-ORF differential translation efficiency contrast.               ##
################################################################################

# Always emit the mandatory tables (even empty) so Nextflow channels are
# consistent; strategy tables are optional and only written when populated.
empty_safe <- function(df) if (is.null(df)) tibble() else df
write_tsv(empty_safe(dte_interaction), paste0(prefix, ".translation.dotseq.results.tsv"))
write_tsv(empty_safe(dou_interaction), paste0(prefix, ".dou.dotseq.results.tsv"))
write_optional <- function(df, suffix) {
    if (!is.null(df) && nrow(df) > 0) write_tsv(df, paste0(prefix, ".", suffix))
}
write_optional(dou_strategy, "dou_strategy.dotseq.results.tsv")
write_optional(dte_strategy, "dte_strategy.dotseq.results.tsv")

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

    # plotDOT's default `force_new_device = TRUE` would reset our png() device;
    # disable it. Return success so the heatmap fallback can distinguish a real
    # plot from one that left a 0-byte file behind.
    safe_plot_dot <- function(plot_type, fname, results_df = NULL, data = NULL,
                              annotation_params = list()) {
        if (is.null(results_df) || nrow(results_df) == 0) return(invisible(FALSE))
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
            invisible(TRUE)
        }, error = \\(e) {
            while (length(dev.list()) > 0) dev.off()
            if (file.exists(fname)) unlink(fname)
            message(sprintf("plotDOT(%s) failed: %s", plot_type, conditionMessage(e)))
            invisible(FALSE)
        })
    }

    plotdot_df <- if (!is.null(dou_interaction) && !is.null(dte_interaction)) {
        dou_interaction |> inner_join(dte_interaction, by = "orf_id", suffix = c("_dou", "_dte"))
    } else NULL

    safe_plot_dot("volcano",   paste0(prefix, ".volcano.png"),   plotdot_df, dou_d)
    safe_plot_dot("composite", paste0(prefix, ".composite.png"), plotdot_df, dou_d)
    safe_plot_dot("venn",      paste0(prefix, ".venn.png"),      plotdot_df)

    # Heatmap needs mORF + sorf_type pairs per gene; try uORF (package default)
    # then dORF.
    heatmap_path <- paste0(prefix, ".heatmap.png")
    if (!safe_plot_dot("heatmap", heatmap_path, plotdot_df, dou_d, list(sorf_type = "uORF"))) {
        safe_plot_dot("heatmap", heatmap_path, plotdot_df, dou_d, list(sorf_type = "dORF"))
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
