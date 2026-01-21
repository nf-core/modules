#!/usr/bin/env Rscript

# Filter rows based on the number of columns passing the abundance threshold. By
# default this will be any row with a value of 1 or more, which would be a
# permissive threshold for RNA-seq data.
#
# In RNA-seq studies it's often not enough to just remove genes not expressed in
# any sample. We also want to remove anything likely to be part of noise, or
# which has sufficiently low expression that differential analysis would not be
# useful. For that reason we might require a higher threshold than 1, and
# require that more than one sample passes.
#
# Often we want to know that a gene is expressed in a substantial enough number
# of sample that differential analysis worthwhile, so we may pick a threshold
# sample number related to group size. Note that we do not filter with an
# awareness of the groups themselves, since this adds bias towards discovery
# between those groups.
#
# This script supports the following options:
#
# - `--minimum_abundance`: Minimum value threshold for feature filtering.
#   Set to a numeric value (e.g., 0.5) to enable, or to 'false'/'null' to disable.
# - `--minimum_samples`: Minimum number of samples passing the abundance threshold.
# - `--minimum_samples_not_na`: Minimum number of non-NA values per feature.
# - `--grouping_variable`: Optional column in sample sheet for group-specific filtering.
# - `--minimum_proportion`: Proportion-based filtering threshold.
# - `--minimum_proportion_not_na`: Minimum proportion of non-NA values required per feature.
# - `--most_variant_features`: Optional integer specifying the number of most-variant features to retain.

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
#' @param nrows Passed to read.delim()
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, nrows = -1 ){

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
        check.names = FALSE
    )
}

#' Convert NULL-like inputs to R NULL.
#'
#' @param v Any input to be checked and parsed.
#'
#' @return NULL object if input is NULL-like, otherwise the original input.
parse_null <- function(v) {
    if (is.null(v)) return(NULL)
    if (length(v) == 0) return(NULL)
    vv <- as.character(v)
    if (vv %in% c("null", "NULL", "", "false", "FALSE")) return(NULL)
    v
}

#' Identify rows that are among the top n most variant
#'
#' @param matrix_data Matrix object
#'
#' @return output Boolean vector

most_variant_test <- function(matrix_data) {

    # Determine the indices of the top variant rows based on variance
    top_indices <- order(-apply(matrix_data, 1, var, na.rm = TRUE))[1:opt\$most_variant_features]

    # Return a boolean vector indicating if each row is among the top variant ones
    1:nrow(matrix_data) %in% top_indices
}

# Set up default options

opt <- list(
    abundance_matrix_file = '$abundance',   # Path to input abundance matrix (TSV/CSV)
    sample_file = '$samplesheet',           # Path to optional sample metadata file
    sample_id_col = NULL,                   # Column name in samplesheet matching matrix columns
    minimum_abundance = 1,                  # Abundance threshold (disabled if NULL/false)
    minimum_samples = 1,                    # Minimum number of samples passing abundance filter
    minimum_proportion = 0,                 # Proportion-based filtering threshold (optional)
    grouping_variable = NULL,               # Grouping variable for stratified filtering (optional)
    minimum_proportion_not_na = 0.5,        # Minimum proportion of non-NA values per feature
    minimum_samples_not_na = NULL,          # Minimum count of non-NA samples per feature (optional)
    most_variant_features = NULL            # Number of most variant rows to keep (optional)
)
opt_types <- lapply(opt, class)

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }else{

        # Allow explicit nulls from Nextflow CLI args
        args_opt[[ao]] <- parse_null(args_opt[[ao]])

        # Preserve classes from defaults where possible (only if not NULL now)
        if (!is.null(args_opt[[ao]]) && !is.null(opt[[ao]])) {
            args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        }
        opt[[ao]] <- args_opt[[ao]]
    }
}

abundance_matrix <- read_delim_flexible(opt\$abundance_matrix_file, row.names = 1)
feature_id_name <- colnames( read_delim_flexible(opt\$abundance_matrix_file, nrows = 1)[1])

# If a sample sheet was specified, validate the matrix against it

if (opt\$sample_file != ''){

    # Read the sample sheet and check against matrix

    samplesheet <- read_delim_flexible(
        opt\$sample_file,
        row.names = ifelse(is.null(opt\$sample_id_col), 1, opt\$sample_id_col)
    )
    missing_samples <- setdiff(rownames(samplesheet), colnames(abundance_matrix))

    if (length(missing_samples) > 0){
        stop(
            paste(
                paste(missing_samples, collapse = ', '),
                'not represented in supplied abundance matrix'
            )
        )
    }else{
        abundance_matrix <- abundance_matrix[,rownames(samplesheet)]
    }
}else{

    # If we're not using a sample sheet to select columns, then at least make
    # sure the ones we have are numeric (some upstream things like the RNA-seq
    # workflow have annotation columns as well)

    numeric_columns <- unlist(lapply(1:ncol(abundance_matrix), function(x) is.numeric(abundance_matrix[,x])))
    abundance_matrix <- abundance_matrix[,numeric_columns]
}

# If we want to define n based on the levels of a grouping variable...

if ((opt\$sample_file != '') && ( ! is.null(opt\$grouping_variable))){

    # Pick a minimum number of samples to pass threshold based on group size

    if (! opt\$grouping_variable %in% colnames(samplesheet)){
        stop(paste(opt\$grouping_variable, "not in supplied sample sheet"))
    }else{
        opt\$minimum_samples <- min(table(samplesheet[[opt\$grouping_variable]]))
        if ( opt\$minimum_proportion > 0){
            opt\$minimum_samples <- opt\$minimum_samples * opt\$minimum_proportion
        }
    }
}else if (opt\$minimum_proportion > 0){

    # Or if we want to define it based on a static proportion of the sample number

    opt\$minimum_samples <- ncol(abundance_matrix) * opt\$minimum_proportion
}

# Also set up filtering for NAs; use by default minimum_proportion_not_na; only
# use minimum_samples_not_na if it is provided (default NULL)
# -->NA test can always use minimum_samples_not_na as this will contain the correct
# value even if the proportion is to be used

if (is.null(opt\$minimum_samples_not_na)) {
    opt\$minimum_samples_not_na <- ncol(abundance_matrix) * opt\$minimum_proportion_not_na
}

# Define the tests

tests <- list(
    'na' = function(x) !any(is.na(x)) || sum(!is.na(x)) >= opt\$minimum_samples_not_na  # check if enough values in row are not NA
)

# Only apply abundance filter if a threshold was provided
if (!is.null(opt\$minimum_abundance)) {
    tests[["abundance"]] = function(x) sum(x >= opt\$minimum_abundance, na.rm = TRUE) >= opt\$minimum_samples
}

# Apply the functions row-wise on the abundance_matrix and store the result in a boolean matrix

boolean_matrix <- do.call(
    rbind,
    lapply(seq_len(nrow(abundance_matrix)), function(i) {
        vapply(tests, function(f) f(abundance_matrix[i, ]), logical(1))
    })
)

rownames(boolean_matrix) <- rownames(abundance_matrix)
colnames(boolean_matrix) <- names(tests)
# Apply the 'most_variant_test' function to identify the most variant rows and add
# the result to the boolean matrix

if (! is.null(opt\$most_variant_features)) {
    most_variant_vectors <- most_variant_test(abundance_matrix)
    boolean_matrix <- cbind(boolean_matrix, most_variant_vectors)
}

# We will retain features passing all tests

keep <- apply(boolean_matrix, 1, all)

# Write out the matrix retaining the specified rows and re-prepending the
# column with the feature identifiers

prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')

thresholds <- data.frame(
    rule = c(
        "minimum_samples_not_na",
        if (!is.null(opt\$minimum_abundance)) c("minimum_abundance", "minimum_samples"),
        if (!is.null(opt\$most_variant_features)) "most_variant_features"
    ),
    threshold = c(
        opt\$minimum_samples_not_na,
        if (!is.null(opt\$minimum_abundance)) c(opt\$minimum_abundance, opt\$minimum_samples),
        if (!is.null(opt\$most_variant_features)) opt\$most_variant_features
    ),
    check.names = FALSE
)

    write.table(
        thresholds,
        file = paste0(prefix, ".thresholds.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )


write.table(
    data.frame(rownames(abundance_matrix)[keep], abundance_matrix[keep,,drop = FALSE]),
    file = paste0(
        prefix,
        '.filtered.tsv'
    ),
    col.names = c(feature_id_name, colnames(abundance_matrix)),
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

# Write a boolean matrix returning specifying the status of each test

write.table(
    data.frame(rownames(abundance_matrix), boolean_matrix),
    file = paste0(prefix, '.tests.tsv'),
    col.names = c(feature_id_name, colnames(boolean_matrix)),
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink("R_sessionInfo.log")
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
