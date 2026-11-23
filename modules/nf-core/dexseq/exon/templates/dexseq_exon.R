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
#' A variable with one level in the samples of the contrast would break the model
#' matrix, so it is left out with a warning.
#'
#' @return Named list of factors, named with make.names()

read_blocking <- function(sample_sheet, blocking_variables) {
    columns <- setdiff(strsplit(blocking_variables, ';', fixed = TRUE)[[1]], '')

    missing_columns <- setdiff(columns, colnames(sample_sheet))
    if (length(missing_columns) > 0) {
        stop("Blocking variables missing from the sample sheet: ", toString(missing_columns))
    }
    if (anyDuplicated(c('sample', 'exon', 'condition', make.names(columns))) > 0) {
        stop("Blocking variables cannot be named 'sample', 'exon' or 'condition', or share a name after make.names(): ", toString(columns))
    }

    incomplete <- columns[vapply(sample_sheet[columns], function(x) any(is.na(x) | x == ''), logical(1))]
    if (length(incomplete) > 0) {
        stop("Blocking variables with missing values in this contrast: ", toString(incomplete))
    }

    blocking <- lapply(setNames(sample_sheet[columns], make.names(columns)), function(x) factor(as.character(x)))
    single <- vapply(blocking, nlevels, integer(1)) < 2
    for (column in columns[single]) {
        warning("Blocking variable '", column, "' has one level in this contrast and is left out of the model")
    }

    blocking[! single]
}

#' Stop if the contrast cannot be told apart from the blocking variables
#'
#' DEXSeq does not check the rank of its models, so this fails with a clear message
#' instead of a degrees of freedom error from the model fit.

check_blocking <- function(blocking, condition, contrast_variable) {
    design_data <- data.frame(condition = condition)
    design_data[names(blocking)] <- blocking
    rank_of <- function(terms) qr(model.matrix(reformulate(c('1', terms)), design_data))[['rank']]

    if (rank_of(c(names(blocking), 'condition')) == rank_of(names(blocking))) {
        stop("Blocking variables ", toString(names(blocking)), " are confounded with '", contrast_variable, "' in this contrast")
    }
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

blocking <- read_blocking(sample_sheet, opt[['blocking_variables']])
check_blocking(blocking, sample_data[['condition']], opt[['contrast_variable']])
sample_data[names(blocking)] <- blocking

# The models have a term per sample, so blocking variables enter as interactions with exon

reduced_terms <- c('sample', 'exon', sprintf('%s:exon', names(blocking)))
full_model <- reformulate(c(reduced_terms, 'condition:exon'))
reduced_model <- reformulate(reduced_terms)
message("Full model: ", deparse1(full_model))

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
    design = full_model,
    flattenedfile = opt[['annotation_file']]
)

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, BPPARAM = SerialParam())
dxd <- testForDEU(dxd, reducedModel = reduced_model, BPPARAM = BPPARAM)
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
