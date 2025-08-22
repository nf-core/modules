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

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

# Set defaults and classes

opt <- list(
    seurat_object = '$seurat_object',
    assay = '$assay',
    output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    quantile = 0.99,
    init = NULL,
    nstarts = 100L,
    kfunc = "clara",
    nsamples = 100L,
    seed = 42L,
    verbose = TRUE
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }else{
        # Handle special cases for logical and NULL values
        if (ao == "verbose") {
            opt[[ao]] <- as.logical(args_opt[[ao]])
        } else if (ao == "init" && args_opt[[ao]] == "NULL") {
            opt[[ao]] <- NULL
        } else if (! is.null(opt[[ao]])){
            # Preserve classes from defaults where possible
            opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        } else {
            opt[[ao]] <- args_opt[[ao]]
        }
    }
}

# Check file inputs are valid

if (! is_valid_string(opt\$seurat_object)) {
    stop("Please provide seurat_object", call. = FALSE)
}

if (! file.exists(opt\$seurat_object)){
    stop(paste0('Value of seurat_object: ', opt\$seurat_object, ' is not a valid file'))
}

if (! is_valid_string(opt\$assay)) {
    stop("Please provide assay", call. = FALSE)
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
hashtag <- readRDS(opt\$seurat_object)

# Demultiplex cells based on HTO enrichment
hashtag <- HTODemux(hashtag,
                   assay = opt\$assay,
                   positive.quantile = opt\$quantile,
                   init = opt\$init,
                   nstarts = opt\$nstarts,
                   kfunc = opt\$kfunc,
                   nsamples = opt\$nsamples,
                   seed = opt\$seed,
                   verbose = opt\$verbose)

################################################
################################################
## SAVING RESULTS                             ##
################################################
################################################

# create a data frame to save the used parameters in a csv file
init_value <- if (is.null(opt\$init)) "NULL" else opt\$init

Argument <- c("seuratObject", "quantile", "kfunc", "nstarts", "nsamples", "seed", "init", "assay", "verbose")
Value <- c(opt\$seurat_object, opt\$quantile, opt\$kfunc, opt\$nstarts, opt\$nsamples, opt\$seed, init_value, opt\$assay, opt\$verbose)
params <- data.frame(Argument, Value)

write.csv(params, paste0(opt\$output_prefix, "_params_htodemux.csv"))

# create csv files to save the results from HTODemux()
donors <- rownames(hashtag[[opt\$assay]])
assignment <- hashtag[[paste0(opt\$assay, "_classification")]]
assignment[[paste0(opt\$assay, "_classification")]][!assignment[[paste0(opt\$assay, "_classification")]] %in% c(donors, "Negative")] <- "Doublet"
write.csv(assignment, paste0(opt\$output_prefix, "_assignment_htodemux.csv"))
write.csv(hashtag[[paste0(opt\$assay, "_classification.global")]], paste0(opt\$output_prefix, "_classification_htodemux.csv"))
saveRDS(hashtag, file = paste0(opt\$output_prefix, "_htodemux.rds"))

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
        paste('    seurat:', seurat.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
