#!/usr/bin/env Rscript

# Author: Samuel Ahuno
# Date: 2026-04-04
# Purpose: Composite plot of modkit localize results — overlay samples colored by condition
#          with configurable smoothing (loess default, rolling_mean, binned)

################################################
## Functions                                  ##
################################################

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#' @return named list of options and values similar to optparse

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})
    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[ ( ! parsed_args %in%  c('', 'null')) & ! is.na(parsed_args)]
}

################################################
## Libraries                                  ##
################################################

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(zoo)
})

################################################
## Parse options from ext.args                ##
################################################

opt <- list(
    prefix         = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    samplesheet    = '$samplesheet',
    smooth         = 'loess',
    loess_span     = 0.1,
    rolling_window = 50,
    bin_size       = 50
)

# Override defaults with ext.args
args_opt <- parse_args('$task.ext.args')

if ('smooth' %in% names(args_opt))         opt\$smooth         <- args_opt\$smooth
if ('loess_span' %in% names(args_opt))     opt\$loess_span     <- as.numeric(args_opt\$loess_span)
if ('rolling_window' %in% names(args_opt)) opt\$rolling_window <- as.integer(args_opt\$rolling_window)
if ('bin_size' %in% names(args_opt))       opt\$bin_size       <- as.integer(args_opt\$bin_size)

cat("Options:\\n")
cat(paste0("  prefix: ", opt\$prefix, "\\n"))
cat(paste0("  samplesheet: ", opt\$samplesheet, "\\n"))
cat(paste0("  smooth: ", opt\$smooth, "\\n"))
cat(paste0("  loess_span: ", opt\$loess_span, "\\n"))
cat(paste0("  rolling_window: ", opt\$rolling_window, "\\n"))
cat(paste0("  bin_size: ", opt\$bin_size, "\\n"))

################################################
## Create output directories                  ##
################################################

fig_dir <- "figures"
dir.create(file.path(fig_dir, "png"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(fig_dir, "pdf"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(fig_dir, "svg"), showWarnings = FALSE, recursive = TRUE)

################################################
## Load samplesheet and TSV files             ##
################################################

# Samplesheet must have columns: sample_id, condition
# Optionally: tsv_path (if not provided, TSVs are matched by filename)
ss <- fread(opt\$samplesheet)
stopifnot("sample_id" %in% names(ss))
stopifnot("condition" %in% names(ss))

cat(paste0("Samplesheet: ", nrow(ss), " samples, conditions: ",
           paste(unique(ss\$condition), collapse = ", "), "\\n"))

# Find all TSV files in current directory (staged by Nextflow)
tsv_files <- list.files(".", pattern = "\\\\.tsv\$", full.names = TRUE)
cat(paste0("Found ", length(tsv_files), " TSV files\\n"))

# Load and annotate each TSV
dat_list <- lapply(tsv_files, function(f) {
    d <- fread(f)
    fname <- basename(f)
    # Match to samplesheet by checking if sample_id is contained in filename
    matched <- ss[sapply(ss\$sample_id, function(sid) grepl(sid, fname, fixed = TRUE)), ]
    if (nrow(matched) == 0) {
        warning(paste0("No samplesheet match for: ", fname, " — skipping"))
        return(NULL)
    }
    d[, `:=`(
        sample_id = matched\$sample_id[1],
        condition = matched\$condition[1],
        filename  = fname
    )]
    return(d)
})

dat <- rbindlist(dat_list[!sapply(dat_list, is.null)])
cat(paste0("Combined data: ", nrow(dat), " rows x ", ncol(dat), " columns\\n"))
cat(paste0("Matched samples: ", dat[, uniqueN(sample_id)], "\\n"))

if (nrow(dat) == 0) {
    stop("No data loaded — check samplesheet sample_id values match TSV filenames")
}

################################################
## Per-condition summary                      ##
################################################

dat_summary <- dat[, .(
    mean_pct = mean(percent_modified, na.rm = TRUE),
    se_pct   = sd(percent_modified, na.rm = TRUE) / sqrt(.N),
    n_samples = .N
), by = .(condition, offset)]
setorder(dat_summary, condition, offset)

################################################
## Smoothing functions                        ##
################################################

apply_smoothing <- function(dt, method, roll_w, bin_sz, l_span) {
    if (method == "loess") {
        cat(paste0("Applying LOESS smoothing (span = ", l_span, ")\\n"))
        out <- dt[, {
            fit_mean <- loess(mean_pct ~ offset, data = .SD, span = l_span)
            # se_pct is NA when only 1 sample per condition; fall back to zero ribbon
            if (all(is.na(se_pct))) {
                smooth_se_vals <- rep(0, nrow(.SD))
            } else {
                sd_sub <- .SD[!is.na(se_pct)]
                fit_se <- loess(se_pct ~ offset, data = sd_sub, span = l_span)
                smooth_se_vals <- predict(fit_se, newdata = .SD)
            }
            .(offset     = offset,
              smooth_pct = predict(fit_mean),
              smooth_se  = smooth_se_vals)
        }, by = condition]

    } else if (method == "rolling_mean") {
        cat(paste0("Applying rolling mean (window = ", roll_w, " bp)\\n"))
        out <- dt[, {
            .(offset     = offset,
              smooth_pct = rollmean(mean_pct, k = roll_w, fill = NA, align = "center"),
              smooth_se  = rollmean(se_pct, k = roll_w, fill = NA, align = "center"))
        }, by = condition]

    } else if (method == "binned") {
        cat(paste0("Applying offset binning (bin size = ", bin_sz, " bp)\\n"))
        out <- dt[, {
            bin_mid <- floor(offset / bin_sz) * bin_sz + bin_sz / 2
            .SD[, .(smooth_pct = mean(mean_pct, na.rm = TRUE),
                     smooth_se  = mean(se_pct, na.rm = TRUE)),
                by = .(offset = bin_mid)]
        }, by = condition]

    } else {
        stop(paste0("Unknown smoothing method: ", method))
    }

    return(out[!is.na(smooth_pct)])
}

dat_smooth <- apply_smoothing(dat_summary, opt\$smooth, opt\$rolling_window, opt\$bin_size, opt\$loess_span)
cat(paste0("Smoothed data: ", nrow(dat_smooth), " rows\\n"))

################################################
## Build labels                               ##
################################################

method_label <- switch(opt\$smooth,
    loess        = paste0("LOESS (span = ", opt\$loess_span, ")"),
    rolling_mean = paste0("Rolling mean (", opt\$rolling_window, " bp window)"),
    binned       = paste0("Binned (", opt\$bin_size, " bp bins)")
)

method_suffix <- switch(opt\$smooth,
    loess        = "loess",
    rolling_mean = paste0("rollmean", opt\$rolling_window, "bp"),
    binned       = paste0("binned", opt\$bin_size, "bp")
)

# Detect unique conditions for dynamic color mapping (Okabe-Ito palette)
conditions <- sort(unique(dat\$condition))
okabe_ito <- c("#D55E00", "#0072B2", "#009E73", "#E69F00", "#56B4E9", "#CC79A7", "#F0E442", "#999999")
condition_colors <- setNames(okabe_ito[seq_along(conditions)], conditions)

n_per_cond <- dat[, .(n = uniqueN(sample_id)), by = condition]
cond_label <- paste(n_per_cond\$n, n_per_cond\$condition, collapse = " + ")

################################################
## Plot 1: Per-sample lines (smoothed)        ##
################################################

cat("Generating composite line plot\\n")

if (opt\$smooth == "loess") {
    p1 <- ggplot(dat, aes(x = offset, y = percent_modified, group = sample_id, color = condition)) +
        geom_smooth(method = "loess", span = opt\$loess_span, se = FALSE, linewidth = 0.8, alpha = 0.6) +
        scale_color_manual(values = condition_colors, name = "Condition") +
        labs(
            title = paste0("Methylation Around Features (", opt\$prefix, ")"),
            subtitle = paste0("Per-sample ", method_label, " | ", cond_label),
            x = "Offset from feature midpoint (bp)",
            y = "Percent modified (%)"
        ) +
        theme_bw(base_size = 25, base_family = "Arial") +
        theme(plot.title = element_text(face = "bold"),
              legend.position = "top",
              panel.grid.minor = element_blank())
} else {
    dat_per_sample <- dat[, {
        setorder(.SD, offset)
        if (opt\$smooth == "rolling_mean") {
            .(offset = offset,
              smooth_pct = rollmean(percent_modified, k = opt\$rolling_window, fill = NA, align = "center"))
        } else {
            bin_mid <- floor(offset / opt\$bin_size) * opt\$bin_size + opt\$bin_size / 2
            .SD[, .(smooth_pct = mean(percent_modified, na.rm = TRUE)), by = .(offset = bin_mid)]
        }
    }, by = .(sample_id, condition)]
    dat_per_sample <- dat_per_sample[!is.na(smooth_pct)]

    p1 <- ggplot(dat_per_sample, aes(x = offset, y = smooth_pct, group = sample_id, color = condition)) +
        geom_line(alpha = 0.6, linewidth = 0.8) +
        scale_color_manual(values = condition_colors, name = "Condition") +
        labs(
            title = paste0("Methylation Around Features (", opt\$prefix, ")"),
            subtitle = paste0("Per-sample ", method_label, " | ", cond_label),
            x = "Offset from feature midpoint (bp)",
            y = "Percent modified (%)"
        ) +
        theme_bw(base_size = 25, base_family = "Arial") +
        theme(plot.title = element_text(face = "bold"),
              legend.position = "top",
              panel.grid.minor = element_blank())
}

################################################
## Plot 2: Mean +/- SE ribbon (smoothed)      ##
################################################

cat("Generating mean +/- SE ribbon plot\\n")

p2 <- ggplot(dat_smooth, aes(x = offset, y = smooth_pct, color = condition, fill = condition)) +
    geom_ribbon(aes(ymin = smooth_pct - smooth_se, ymax = smooth_pct + smooth_se),
                alpha = 0.25, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = condition_colors, name = "Condition") +
    scale_fill_manual(values = condition_colors, name = "Condition") +
    labs(
        title = paste0("Methylation Around Features (", opt\$prefix, ")"),
        subtitle = paste0("Mean +/- SE | ", method_label, " | ", cond_label),
        x = "Offset from feature midpoint (bp)",
        y = "Percent modified (%)"
    ) +
    theme_bw(base_size = 25, base_family = "Arial") +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top",
          panel.grid.minor = element_blank())

################################################
## Save figures                               ##
################################################

cat("Saving figures\\n")

save_fig <- function(p, prefix, w = 12, h = 8) {
    ggsave(file.path(fig_dir, "png", paste0(prefix, ".png")), p, width = w, height = h, dpi = 300)
    ggsave(file.path(fig_dir, "pdf", paste0(prefix, ".pdf")), p, width = w, height = h, device = cairo_pdf)
    ggsave(file.path(fig_dir, "svg", paste0(prefix, ".svg")), p, width = w, height = h)
    cat(paste0("  Saved: ", prefix, " (png, pdf, svg)\\n"))
}

save_fig(p1, paste0(opt\$prefix, "_composite_lines_", method_suffix))
save_fig(p2, paste0(opt\$prefix, "_mean_ribbon_", method_suffix))

################################################
## Save combined data                         ##
################################################

fwrite(dat[, .(sample_id, condition, mod_code, offset, n_valid, n_mod, percent_modified)],
       paste0(opt\$prefix, "_combined.tsv"), sep = "\t")
cat(paste0("Saved combined TSV: ", nrow(dat), " rows\\n"))

fwrite(dat_summary, paste0(opt\$prefix, "_summary.tsv"), sep = "\t")
cat(paste0("Saved summary TSV: ", nrow(dat_summary), " rows\\n"))

################################################
## Versions                                   ##
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
dt.version <- as.character(packageVersion('data.table'))
gg.version <- as.character(packageVersion('ggplot2'))
zoo.version <- as.character(packageVersion('zoo'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    r-data.table:', dt.version),
        paste('    r-ggplot2:', gg.version),
        paste('    r-zoo:', zoo.version)
    ),
    'versions.yml'
)

cat("Done\\n")
