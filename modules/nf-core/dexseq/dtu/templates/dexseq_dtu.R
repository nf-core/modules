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

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

opt <- list(
    output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    sample_file = '$samplesheet',
    count_file = '$counts',
    contrast_variable = '$contrast_variable',
    reference_level = '$reference',
    target_level = '$target',
    sample_id_col = 'sample',
    feature_id_col = 'feature_id',
    gene_id_col = 'gene_id',
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

count_table <- read_delim_flexible(opt[['count_file']])

for (column in c(opt[['feature_id_col']], opt[['gene_id_col']])) {
    if (! column %in% colnames(count_table)) {
        stop(paste0("Column '", column, "' is not in the count table"))
    }
}

missing_samples <- setdiff(sample_data[['sample']], colnames(count_table))

if (length(missing_samples) > 0) {
    stop(paste("Samples missing from the count table:", paste(missing_samples, collapse = ", ")))
}

# DEXSeq models integer counts, so estimated counts such as those of tximport are rounded

counts <- round(as.matrix(count_table[, sample_data[['sample']], drop = FALSE]))
rownames(counts) <- count_table[[opt[['feature_id_col']]]]

################################################
################################################
## RUN DEXSEQ                                 ##
################################################
################################################

dxd <- DEXSeqDataSet(
    countData = counts,
    sampleData = sample_data,
    design = ~ sample + exon + condition:exon,
    featureID = count_table[[opt[['feature_id_col']]]],
    groupID = count_table[[opt[['gene_id_col']]]]
)

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, BPPARAM = SerialParam())
dxd <- testForDEU(dxd, BPPARAM = BPPARAM)
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
