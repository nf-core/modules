#!/usr/bin/env Rscript

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x) {
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text = x, what = 'character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z) { length(z) <- 2; z })

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#'
#' @return Data frame

read_delim_flexible <- function(file) {
    ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1))

    if (ext == "tsv" || ext == "txt") {
        separator <- "\\t"
    } else if (ext == "csv") {
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }

    read.delim(file, sep = separator, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

#' Keep the result columns that can be written to a flat table
#'
#' @param x DEXSeqResults object
#'
#' @return DEXSeqResults object without the list columns

keep_valid_columns <- function(x) {
    keep <- vapply(x, function(col) is.numeric(col) || is.character(col), FUN.VALUE = logical(1))
    x[, keep, drop = FALSE]
}

#' Read the blocking variables of the contrast from the sample sheet
#'
#' A variable with a single level in the samples of the contrast has no effect to
#' adjust for, so it is dropped with a warning.
#'
#' @param sample_sheet Sample sheet holding only the samples of the contrast
#' @param blocking_variables Semicolon-separated sample sheet columns
#' @param contrast_variable Sample sheet column holding the contrast levels
#'
#' @return Named list of factors, named with make.names()

read_blocking <- function(sample_sheet, blocking_variables, contrast_variable) {
    columns <- trimws(unlist(strsplit(blocking_variables, ';', fixed = TRUE)))
    columns <- unique(columns[columns != ''])

    missing_columns <- setdiff(columns, colnames(sample_sheet))
    if (length(missing_columns) > 0) {
        stop(paste("Blocking variables missing from the sample sheet:", paste(missing_columns, collapse = ", ")))
    }

    names_model <- make.names(columns)
    clashes <- columns[columns == contrast_variable | names_model %in% c('sample', 'exon', 'condition')]
    if (length(clashes) > 0) {
        stop(paste0(
            "Blocking variables cannot be the contrast variable or be named 'sample', 'exon' or 'condition': ",
            paste(clashes, collapse = ", ")
        ))
    }
    if (anyDuplicated(names_model) > 0) {
        stop(paste("Blocking variables with the same name after make.names():", paste(columns, collapse = ", ")))
    }

    blocking <- list()

    for (i in seq_along(columns)) {
        values <- as.character(sample_sheet[[columns[i]]])

        if (any(is.na(values) | values == '')) {
            stop(paste0("Blocking variable '", columns[i], "' has missing values for samples in this contrast"))
        }

        values <- factor(values)

        if (nlevels(values) < 2) {
            warning(paste0("Blocking variable '", columns[i], "' has one level in this contrast and is left out of the model"))
            next
        }

        blocking[[names_model[i]]] <- values
    }

    blocking
}

#' Stop if the blocking variables cannot be estimated next to the contrast
#'
#' @param blocking Named list of factors from read_blocking()
#' @param condition Factor of the contrast levels
#' @param contrast_variable Sample sheet column holding the contrast levels

check_blocking <- function(blocking, condition, contrast_variable) {
    if (length(blocking) == 0) {
        return(invisible(NULL))
    }

    design_data <- as.data.frame(c(blocking, list(condition = condition)))
    rank_of <- function(terms) qr(model.matrix(reformulate(terms), data = design_data))[['rank']]
    ncol_of <- function(terms) ncol(model.matrix(reformulate(terms), data = design_data))

    for (b in names(blocking)) {
        if (rank_of(c(b, 'condition')) < ncol_of(c(b, 'condition'))) {
            stop(paste0(
                "Blocking variable '", b, "' is confounded with '", contrast_variable,
                "' in this contrast, so it cannot be adjusted for"
            ))
        }
    }

    if (rank_of(names(blocking)) < ncol_of(names(blocking))) {
        stop(paste(
            "The blocking variables are confounded with each other in this contrast:",
            paste(names(blocking), collapse = ", ")
        ))
    }

    terms <- c(names(blocking), 'condition')
    if (rank_of(terms) < ncol_of(terms)) {
        stop(paste0(
            "The blocking variables together are confounded with '", contrast_variable,
            "' in this contrast: ", paste(names(blocking), collapse = ", ")
        ))
    }
    if (rank_of(terms) >= nrow(design_data)) {
        stop(paste(
            "Too few samples in this contrast to adjust for the blocking variables:",
            paste(names(blocking), collapse = ", ")
        ))
    }

    invisible(NULL)
}

#' Build the full and reduced DEXSeq models
#'
#' The models already have a term per sample, so the blocking variables enter as
#' interactions with exon in both models, as the contrast does in the full model.
#'
#' @param blocking Named list of factors from read_blocking()
#'
#' @return List with the full and the reduced formula

dexseq_models <- function(blocking) {
    blocking_terms <- sprintf('%s:exon', names(blocking))
    reduced <- c('sample', 'exon', blocking_terms)

    list(
        full = reformulate(c(reduced, 'condition:exon')),
        reduced = reformulate(reduced)
    )
}

#' Plot the exon usage of the genes with the smallest p-values
#'
#' @param x DEXSeqResults object
#' @param file Path of the PDF file to write
#' @param n Number of genes to plot

save_plot_dexseq <- function(x, file, n) {
    x_valid <- x[is.finite(x[['pvalue']]) & is.finite(x[['padj']]), , drop = FALSE]

    if (nrow(x_valid) == 0) {
        warning("No valid results to plot.")
        return(invisible(NULL))
    }

    n <- min(n, nrow(x_valid))
    ids <- x_valid[['groupID']][order(x_valid[['pvalue']])[seq_len(n)]]

    pdf(file = file, width = 11.7, height = 8.3)
    on.exit(dev.off())

    for (i in ids) {
        tryCatch(
            plotDEXSeq(x, geneID = i, displayTranscripts = TRUE),
            error = function(e) warning(paste("Error plotting gene", i, ":", conditionMessage(e)))
        )
    }

    invisible(NULL)
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

opt <- list(
    output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    sample_file = '$samplesheet',
    count_dir = 'counts',
    annotation_file = '$gff',
    contrast_variable = '$contrast_variable',
    reference_level = '$reference',
    target_level = '$target',
    blocking_variables = '',
    sample_id_col = 'sample',
    count_suffix = '.txt',
    n_plot = 10,
    cores = as.integer('$task.cpus')
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for (ao in names(args_opt)) {
    if (! ao %in% names(opt)) {
        stop(paste("Invalid option:", ao))
    }
    opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
}

if (is.na(opt[['n_plot']]) || opt[['n_plot']] < 0) {
    stop("Option --n_plot must be a non-negative number")
}

if (opt[['reference_level']] == opt[['target_level']]) {
    stop("The reference and target levels must differ")
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(BiocParallel)
library(DEXSeq)

if (! is.na(opt[['cores']]) && opt[['cores']] > 1) {
    BPPARAM <- MulticoreParam(workers = opt[['cores']])
} else {
    BPPARAM <- SerialParam()
}

################################################
################################################
## READ IN SAMPLE SHEET AND COUNTS            ##
################################################
################################################

sample_sheet <- read_delim_flexible(opt[['sample_file']])

for (column in c(opt[['sample_id_col']], opt[['contrast_variable']])) {
    if (! column %in% colnames(sample_sheet)) {
        stop(paste0("Column '", column, "' is not in the sample sheet"))
    }
}

# The sample sheet can have several rows per sample, e.g. for sequencing runs

sample_sheet <- sample_sheet[! duplicated(sample_sheet[[opt[['sample_id_col']]]]), , drop = FALSE]

levels <- c(opt[['reference_level']], opt[['target_level']])
missing_levels <- setdiff(levels, sample_sheet[[opt[['contrast_variable']]]])

if (length(missing_levels) > 0) {
    stop(paste0(
        "Levels missing from the '", opt[['contrast_variable']], "' column of the sample sheet: ",
        paste(missing_levels, collapse = ", ")
    ))
}

# Keep the samples of the two levels being compared, with the reference first so
# that the fold changes are those of the target over the reference

sample_sheet <- sample_sheet[sample_sheet[[opt[['contrast_variable']]]] %in% levels, , drop = FALSE]

sample_data <- data.frame(
    sample = as.character(sample_sheet[[opt[['sample_id_col']]]]),
    condition = factor(sample_sheet[[opt[['contrast_variable']]]], levels = levels)
)

blocking <- read_blocking(sample_sheet, opt[['blocking_variables']], opt[['contrast_variable']])
check_blocking(blocking, sample_data[['condition']], opt[['contrast_variable']])

for (b in names(blocking)) {
    sample_data[[b]] <- blocking[[b]]
}

models <- dexseq_models(blocking)
message('Full model: ', paste(deparse(models[['full']], width.cutoff = 500L), collapse = ' '))
message('Reduced model: ', paste(deparse(models[['reduced']], width.cutoff = 500L), collapse = ' '))

count_files <- file.path(opt[['count_dir']], paste0(sample_data[['sample']], opt[['count_suffix']]))
missing_files <- count_files[! file.exists(count_files)]

if (length(missing_files) > 0) {
    stop(paste("Missing count files:", paste(basename(missing_files), collapse = ", ")))
}

################################################
################################################
## RUN DEXSEQ                                 ##
################################################
################################################

dxd <- DEXSeqDataSetFromHTSeq(
    countfiles = count_files,
    sampleData = sample_data,
    design = models[['full']],
    flattenedfile = opt[['annotation_file']]
)

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, BPPARAM = SerialParam())
dxd <- testForDEU(dxd, fullModel = models[['full']], reducedModel = models[['reduced']], BPPARAM = BPPARAM)
dxd <- estimateExonFoldChanges(dxd, BPPARAM = BPPARAM)

results <- DEXSeqResults(dxd)
gene_qvalues <- perGeneQValue(results)

################################################
################################################
## OUTPUTS                                    ##
################################################
################################################

prefix <- opt[['output_prefix']]

saveRDS(dxd, file = paste0(prefix, '.dexseq.dataset.rds'))
saveRDS(results, file = paste0(prefix, '.dexseq.results.rds'))
saveRDS(gene_qvalues, file = paste0(prefix, '.dexseq.gene_qvalues.rds'))

write.table(
    as.data.frame(keep_valid_columns(results)),
    file = paste0(prefix, '.dexseq.results.tsv'),
    sep = '\\t',
    quote = FALSE,
    row.names = FALSE
)

write.table(
    data.frame(groupID = names(gene_qvalues), padj = unname(gene_qvalues)),
    file = paste0(prefix, '.dexseq.gene_qvalues.tsv'),
    sep = '\\t',
    quote = FALSE,
    row.names = FALSE
)

if (opt[['n_plot']] > 0) {
    save_plot_dexseq(results, file = paste0(prefix, '.dexseq.plots.pdf'), n = opt[['n_plot']])
}

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste0(prefix, '.R_sessionInfo.log'))
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', paste(R.version[['major']], R.version[['minor']], sep = '.')),
        paste('    bioconductor-dexseq:', as.character(packageVersion('DEXSeq')))
    ),
    'versions.yml'
)
