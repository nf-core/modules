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

string_to_logical <- function(input) {
  if (input == "FALSE") {
    FALSE
  } else if (input == "TRUE") {
    TRUE
  } else {
    stop(paste0(input, " is not a valid logical. Use 'FALSE' or 'TRUE'."))
  }
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

# Set defaults and classes
opt <- list(
    # File inputs
    seurat_object = '$seurat_object',
    assay = '$assay',

    # MULTIseqDemux Parameters
    quantile = 0.7,         # A numeric scalar for the quantile to use for demultiplexing
    autoThresh = TRUE,      # A logical scalar indicating whether to use automatic thresholding
    maxiter = 5L,           # An integer scalar specifying the maximum number of iterations
    qrangeFrom = 0.1,       # A numeric scalar specifying the start of the quantile range
    qrangeTo = 0.9,         # A numeric scalar specifying the end of the quantile range
    qrangeBy = 0.05,        # A numeric scalar specifying the step size for the quantile range
    verbose = TRUE,         # A logical scalar indicating whether to print verbose output

    # others
    prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix') # Prefix name for output files.
)
opt_types <- lapply(opt, class)

# Apply parameter overrides
args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }else{
        # Handle special cases for logicals
        if (opt_types[[ao]] == "logical") {
            opt[[ao]] <- string_to_logical(args_opt[[ao]])
        } else if (! is.null(opt[[ao]])){
            # Preserve classes from defaults where possible
            opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        } else {
            opt[[ao]] <- args_opt[[ao]]
        }
    }
}

# Set individual variables for backward compatibility and cleaner code
seuratObj <- opt\$seurat_object
assay <- opt\$assay
quantile <- opt\$quantile
autoThresh <- opt\$autoThresh
maxiter <- opt\$maxiter
qrangeFrom <- opt\$qrangeFrom
qrangeTo <- opt\$qrangeTo
qrangeBy <- opt\$qrangeBy
verbose <- opt\$verbose
prefix <- opt\$prefix

# Configure output precision
options(digits=5)

# Check if file exists
if (! file.exists(seuratObj)){
    stop(paste0(seuratObj, ' is not a valid file'))
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(Seurat)

################################################
################################################
## Main Process                               ##
################################################
################################################

# Loading Seurat object
hashtag <- readRDS(seuratObj)

# Demultiplex cells
if (autoThresh == TRUE) {
  hashtag <- MULTIseqDemux(hashtag, assay = assay, quantile = quantile, autoThresh = TRUE, maxiter = maxiter, qrange = seq(from = qrangeFrom, to = qrangeTo, by = qrangeBy), verbose = verbose)
} else {
  hashtag <- MULTIseqDemux(hashtag, assay = assay, quantile = quantile, verbose = verbose)
}

################################################
################################################
## SAVING RESULTS                             ##
################################################
################################################

# create a data frame to save the used parameters in a csv file
Argument <- c("seuratObjectPath", "quantile", "autoThresh", "maxiter", "qrangeFrom", "qrangeTo", "qrangeBy", "verbose", "assay")
Value <- c(seuratObj, quantile, autoThresh, maxiter, qrangeFrom, qrangeTo, qrangeBy, verbose, assay)
params <- data.frame(Argument, Value)
write.csv(params, paste0(prefix ,"_params_multiseqdemux.csv"))

# save the results from MULTIseqDemux()
write.csv(hashtag\$MULTI_ID, paste0(prefix , "_res_multiseqdemux.csv"))
saveRDS(hashtag, file = paste0(prefix ,"_multiseqdemux.rds"))

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
seurat.version <- as.character(packageVersion('Seurat'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    r-seurat:', seurat.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
