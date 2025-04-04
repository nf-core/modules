#!/usr/bin/env Rscript

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Check for Non-Empty, Non-Whitespace String
#'
#' This function checks if the input is non-NULL and contains more than just whitespace.
#' It returns TRUE if the input is a non-empty, non-whitespace string, and FALSE otherwise.
#'
#' @param input A variable to check.
#' @return A logical value: TRUE if the input is a valid, non-empty, non-whitespace string; FALSE otherwise.
#' @examples
#' is_valid_string("Hello World") # Returns TRUE
#' is_valid_string("   ")         # Returns FALSE
#' is_valid_string(NULL)          # Returns FALSE

is_valid_string <- function(input) {
    !is.null(input) && nzchar(trimws(input))
}

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = TRUE){

    ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1))

    if (ext == "tsv" || ext == "txt") {
        separator <- "\\t"
    } else if (ext == "csv") {
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }

    read.delim(
        file,
        sep = separator,
        header = header,
        row.names = row.names,
        check.names = check.names
    )
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

# I've defined these in a single array like this so that we could go back to an
# optparse-driven method in future with module bin/ directories, rather than
# the template

# Set defaults and classes

opt <- list(
    output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    count_file = '$counts',
    sample_file = '$samplesheet',
    contrast_variable = '$contrast_variable',
    reference_level = '$reference',
    target_level = '$target',
    blocking_variables = NULL,
    control_genes_file = '$control_genes_file',
    transcript_lengths_file = '$transcript_lengths_file',
    sizefactors_from_controls = FALSE,
    gene_id_col = "gene_id",
    sample_id_col = "experiment_accession",
    subset_to_contrast_samples = FALSE,
    exclude_samples_col = NULL,
    exclude_samples_values = NULL,
    test = "Wald",
    fit_type = "parametric",
    sf_type = 'ratio',
    min_replicates_for_replace = 7,
    use_t = FALSE,
    lfc_threshold = 0,
    alt_hypothesis = 'greaterAbs',
    independent_filtering = TRUE,
    p_adjust_method = 'BH',
    alpha = 0.1,
    minmu = 0.5,
    vs_method = 'vst', # 'rlog', 'vst', or 'rlog,vst'
    shrink_lfc = TRUE,
    cores = 1,
    vs_blind = TRUE,
    vst_nsub = 1000,
    round_digits = NULL
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }else{

        # Preserve classes from defaults where possible
        if (! is.null(opt[[ao]])){
            args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        }
        opt[[ao]] <- args_opt[[ao]]
    }
}
if ( ! is.null(opt\$round_digits)){
    opt\$round_digits <- as.numeric(opt\$round_digits)
}

# Check if required parameters have been provided

required_opts <- c('contrast_variable', 'reference_level', 'target_level', 'output_prefix')
missing <- required_opts[!unlist(lapply(opt[required_opts], is_valid_string)) | !required_opts %in% names(opt)]

if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid

for (file_input in c('count_file', 'sample_file')){
    if (! is_valid_string(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }

    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(DESeq2)
library(BiocParallel)

################################################
################################################
## READ IN COUNTS FILE AND SAMPLE METADATA    ##
################################################
################################################

count.table <-
    read_delim_flexible(
        file = opt\$count_file,
        header = TRUE,
        row.names = opt\$gene_id_col,
        check.names = FALSE
    )
sample.sheet <- read_delim_flexible(file = opt\$sample_file)

# Deal with spaces that may be in sample column
opt\$sample_id_col <- make.names(opt\$sample_id_col)

if (! opt\$sample_id_col %in% colnames(sample.sheet)){
    stop(paste0("Specified sample ID column '", opt\$sample_id_col, "' is not in the sample sheet"))
}

# Sample sheet can have duplicate rows for multiple sequencing runs, so uniqify
# before assigning row names

sample.sheet <- sample.sheet[! duplicated(sample.sheet[[opt\$sample_id_col]]), ]
rownames(sample.sheet) <- sample.sheet[[opt\$sample_id_col]]

# Check that all samples specified in the input sheet are present in the counts
# table. Assuming they are, subset and sort the count table to match the sample
# sheet

missing_samples <-
    sample.sheet[!rownames(sample.sheet) %in% colnames(count.table), opt\$sample_id_col]

if (length(missing_samples) > 0) {
    stop(paste(
        length(missing_samples),
        'specified samples missing from count table:',
        paste(missing_samples, collapse = ',')
    ))
} else{
    # Save any non-count data, will gene metadata etc we might need later
    noncount.table <-
        count.table[, !colnames(count.table) %in% rownames(sample.sheet), drop = FALSE]
    count.table <- count.table[, rownames(sample.sheet)]
}

################################################
################################################
## CHECK CONTRAST SPECIFICATION               ##
################################################
################################################

contrast_variable <- make.names(opt\$contrast_variable)
blocking.vars <- c()

if (!contrast_variable %in% colnames(sample.sheet)) {
    stop(
        paste0(
        'Chosen contrast variable \"',
        contrast_variable,
        '\" not in sample sheet'
        )
    )
} else if (any(!c(opt\$reflevel, opt\$treatlevel) %in% sample.sheet[[contrast_variable]])) {
    stop(
        paste(
        'Please choose reference and treatment levels that are present in the',
        contrast_variable,
        'column of the sample sheet'
        )
    )
} else if (is_valid_string(opt\$blocking_variables)) {
    blocking.vars = make.names(unlist(strsplit(opt\$blocking_variables, split = ';')))
    if (!all(blocking.vars %in% colnames(sample.sheet))) {
        missing_block <- paste(blocking.vars[! blocking.vars %in% colnames(sample.sheet)], collapse = ',')
        stop(
            paste(
                'Blocking variables', missing_block,
                'do not correspond to sample sheet columns.'
            )
        )
    }
}

# Optionally, subset to only the samples involved in the contrast

if (opt\$subset_to_contrast_samples){
    sample_selector <- sample.sheet[[contrast_variable]] %in% c(opt\$target_level, opt\$reference_level)
    selected_samples <- sample.sheet[sample_selector, opt\$sample_id_col]
    count.table <- count.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
}

# Optionally, remove samples with specified values in a given field (probably
# don't use this as well as the above)

if ((is_valid_string(opt\$exclude_samples_col)) && (is_valid_string(opt\$exclude_samples_values))){
    exclude_values = unlist(strsplit(opt\$exclude_samples_values, split = ';'))

    if (! opt\$exclude_samples_col %in% colnames(sample.sheet)){
        stop(paste(opt\$exclude_samples_col, ' specified to subset samples is not a valid sample sheet column'))
    }

    print(paste0('Excluding samples with values of ', opt\$exclude_samples_values, ' in ', opt\$exclude_samples_col))
    sample_selector <- ! sample.sheet[[opt\$exclude_samples_col]] %in% exclude_values

    selected_samples <- sample.sheet[sample_selector, opt\$sample_id_col]
    count.table <- count.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
}

# Now specify the model. Use cell-means style so we can be explicit with the
# contrasts

model <- '~ 0'

if (is_valid_string(opt\$blocking_variables)) {
    model <- paste(model, paste(blocking.vars, collapse = ' + '), sep=' + ')
}

# Make sure all the appropriate variables are factors

for (v in c(blocking.vars, contrast_variable)) {
    sample.sheet[[v]] <- as.factor(sample.sheet[[v]])
}

# Variable of interest goes last, see
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#multi-factor-designs

model <- paste(model, contrast_variable, sep = ' + ')

################################################
################################################
## Run DESeq2 processes                       ##
################################################
################################################

if (opt\$control_genes_file != ''){
    control_genes <- readLines(opt\$control_genes_file)
    if (! opt\$sizefactors_from_controls){
        count.table <- count.table[setdiff(rownames(count.table), control_genes),]
    }
}

dds <- DESeqDataSetFromMatrix(
    countData = round(count.table),
    colData = sample.sheet,
    design = as.formula(model)
)

# Build in transcript lengths. Copying what tximport does here:
# https://github.com/thelovelab/DESeq2/blob/6947d5bc629015fb8ffb2453a91b71665a164483/R/AllClasses.R#L409

if (opt\$transcript_lengths_file != ''){
    lengths <-
        read_delim_flexible(
            file = opt\$transcript_lengths_file,
            header = TRUE,
            row.names = opt\$gene_id_col,
            check.names = FALSE
        )
    lengths <- lengths[rownames(count.table), colnames(count.table)]
    dimnames(lengths) <- dimnames(dds)
    assays(dds)[["avgTxLength"]] <- lengths
}

if (opt\$control_genes_file != '' && opt\$sizefactors_from_controls){
    print(paste('Estimating size factors using', length(control_genes), 'control genes'))
    dds <- estimateSizeFactors(dds, controlGenes=rownames(count.table) %in% control_genes)
}

dds <- DESeq(
    dds,
    test = opt\$test,
    fitType = opt\$fit_type,
    minReplicatesForReplace = opt\$min_replicates_for_replace,
    useT = opt\$use_t,
    sfType = opt\$sf_type,
    parallel=TRUE, BPPARAM=MulticoreParam(opt\$cores)
)

comp.results <-
    results(
        dds,
        lfcThreshold = opt\$lfc_threshold,
        altHypothesis = opt\$alt_hypothesis,
        independentFiltering = opt\$independent_filtering,
        alpha = opt\$alpha,
        pAdjustMethod = opt\$p_adjust_method,
        minmu = opt\$minmu,
        contrast = c(
            contrast_variable,
            c(opt\$target_level, opt\$reference_level)
        )
    )

if (opt\$shrink_lfc){
    comp.results <- lfcShrink(dds,
        type = 'ashr',
        contrast = c(
            contrast_variable,
            c(opt\$target_level, opt\$reference_level)
        )
    )
}

# See https://support.bioconductor.org/p/97676/

if (opt\$transcript_lengths_file != ''){
    size_factors = estimateSizeFactorsForMatrix(counts(dds) / assays(dds)[["avgTxLength"]])
}else {
    size_factors = sizeFactors(dds)
}

################################################
################################################
## Generate outputs                           ##
################################################
################################################

contrast.name <-
    paste(opt\$target_level, opt\$reference_level, sep = "_vs_")
cat("Saving results for ", contrast.name, " ...\n", sep = "")

# Differential expression table- note very limited rounding for consistency of
# results

if (! is.null(opt\$round_digits)){
    comp.results <- apply(data.frame(comp.results), 2, function(x) round(x, opt\$round_digits))
}
comp.results <- `colnames<-`(
    data.frame(
        gene_id = rownames(comp.results),
        comp.results,
        check.names = FALSE
    ),
    c(opt\$gene_id_col, colnames(comp.results))  # Setting all column names
)

write.table(
    comp.results,
    file = paste(opt\$output_prefix, 'deseq2.results.tsv', sep = '.'),
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

# Dispersion plot

png(
    file = paste(opt\$output_prefix, 'deseq2.dispersion.png', sep = '.'),
    width = 600,
    height = 600
)
plotDispEsts(dds)
dev.off()

# R object for other processes to use

saveRDS(dds, file = paste(opt\$output_prefix, 'dds.rld.rds', sep = '.'))

# Size factors

sf_df = data.frame(
    sample = names(size_factors),
    data.frame(size_factors, check.names = FALSE),
    check.names = FALSE
)
colnames(sf_df) <- c('sample', 'sizeFactor')
write.table(
    sf_df,
    file = paste(opt\$output_prefix, 'deseq2.sizefactors.tsv', sep = '.'),
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

# Write specified matrices

normalised_matrix <- counts(dds, normalized = TRUE)
if (! is.null(opt\$round_digits)){
    normalised_matrix <- apply(normalised_matrix, 2, function(x) round(x, opt\$round_digits))
}
normalised_matrix <- `colnames<-`(
    data.frame(
        gene_id = rownames(counts(dds)),  # First column with row names from counts(dds)
        normalised_matrix,                # Other columns
        check.names = FALSE
    ),
    c(opt\$gene_id_col, colnames(normalised_matrix))  # Setting all column names
)

write.table(
    normalised_matrix,
    file = paste(opt\$output_prefix, 'normalised_counts.tsv', sep = '.'),
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

# Note very limited rounding for consistency of results

for (vs_method_name in strsplit(opt\$vs_method, ',')){
    if (vs_method_name == 'vst'){
        vs_mat <- assay(vst(dds, blind = opt\$vs_blind, nsub = opt\$vst_nsub))
    }else if (vs_method_name == 'rlog'){
        vs_mat <- assay(rlog(dds, blind = opt\$vs_blind, fitType = opt\$fit_type))
    }

    if (! is.null(opt\$round_digits)){
        vs_mat <- apply(vs_mat, 2, function(x) round(x, opt\$round_digits))
    }

    vs_mat <- `colnames<-`(
        data.frame(
            gene_id = rownames(counts(dds)),  # First column with row names from counts(dds)
            vs_mat,                           # Other columns from vs_mat
            check.names = FALSE
        ),
        c(opt\$gene_id_col, colnames(vs_mat))  # Setting all column names
    )

    write.table(
        vs_mat,
        file = paste(opt\$output_prefix, vs_method_name, 'tsv', sep = '.'),
        col.names = TRUE,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )
}

# Save model to file

write(model, file=paste(opt\$output_prefix, 'deseq2.model.txt', sep = '.'))

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(opt\$output_prefix, "R_sessionInfo.log", sep = '.'))
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
deseq2.version <- as.character(packageVersion('DESeq2'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-deseq2:', deseq2.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
