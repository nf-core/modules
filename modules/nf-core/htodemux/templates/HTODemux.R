#!/usr/bin/env Rscript

################################################
################################################
<<<<<<< HEAD
## USE PARAMETERS FROM NEXTFLOW               ##
################################################
################################################

# cast parameters from nextflow
seuratObj = '$seurat_object'
assay = '$assay'
options(digits=5)
quantile = as.double('$quantile')
init = NULL
if ('$init' != "NULL") {
    init = as.integer('$init')
}
nstarts = as.integer('$nstarts')
kfunc = '$kfunc'
nsamples = as.integer('$nsamples')
seed = as.integer('$seed')
verbose = as.logical('$verbose')
prefix = '$prefix'

# check if the file exists
if (! file.exists(seuratObj)){
    stop(paste0(seuratObj, ' is not a valid file'))
}

=======
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

# I've defined these in a single array like this so that we could go back to an
# optparse-driven method in future with module bin/ directories, rather than
# the template
# Set defaults and classes

opt <- list(
    seuratObj = '$seurat_object',
    assay = 'HTO',
    quantile = 0.99,
    init = NULL,
    nstarts = 100,
    kfunc = 'clara',
    nsamples = 100,
    seed = 42,
    prefix = '$prefix'
)
opt_types <- lapply(opt, class)

# apply parameter overrides
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
if ( ! is.null(opt\$init)){
    opt\$init <- as.numeric(opt\$init)
}

# Check if required parameters have been provided
required_opts <- c('seuratObj')
missing <- required_opts[!unlist(lapply(opt[required_opts], is_valid_string)) | !required_opts %in% names(opt)]
if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# check if the input path is a valid string
if (! is_valid_string(opt[['seuratObj']])) {
    stop(paste("Please provide", file_input), call. = FALSE)
}

# check if the file exists
if (! file.exists(opt[['seuratObj']])){
    stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
}


>>>>>>> df971f6e6 (add template)
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
<<<<<<< HEAD
hashtag <- readRDS(seuratObj)

# Demultiplex cells based on HTO enrichment
hashtag <- HTODemux(hashtag, assay = assay, positive.quantile = quantile, init = init, nstarts = nstarts, kfunc = kfunc, seed = seed, verbose = verbose)

=======
hashtag <- readRDS(opt\$seuratObj)

# Demultiplex cells based on HTO enrichment
if (opt\$kfunc == "clara") {
    hashtag <- HTODemux(hashtag, assay = opt\$assay, positive.quantile = opt\$quantile, init = opt\$init, nstarts = opt\$nstarts, kfunc = "clara", seed = opt\$seed)
} else {
    hashtag <- HTODemux(hashtag, assay = opt\$assay, positive.quantile = opt\$quantile, init = opt\$init, nstarts = opt\$nstarts, kfunc = "kmeans", seed = opt\$seed)
}
>>>>>>> df971f6e6 (add template)

################################################
################################################
## SAVING RESULTS                             ##
################################################
################################################

# create a data frame to save the used parameters in a csv file
<<<<<<< HEAD
if (is.null(init)) {
  init <- "NULL"
}

Argument <- c("seuratObject", "quantile", "kfunc", "nstarts", "nsamples", "seed", "init", "assay", "verbose")
Value <- c(seuratObj, quantile, kfunc, nstarts, nsamples, seed, init, assay, verbose)
params <- data.frame(Argument, Value)

write.csv(params, paste0(prefix ,"_params_htodemux.csv"))

# create csv files to save the results from HTODemux()
donors <- rownames(hashtag[[assay]])
assignment <- hashtag[[paste0(assay, "_classification")]]
assignment[[paste0(assay, "_classification")]][!assignment[[paste0(assay, "_classification")]] %in% c(donors, "Negative")] <- "Doublet"
write.csv(assignment, paste0(prefix ,"_assignment_htodemux.csv"))
write.csv(hashtag[[paste0(assay, "_classification.global")]], paste0(prefix ,"_classification_htodemux.csv"))
saveRDS(hashtag, file = paste0(prefix ,"_htodemux.rds"))
=======
if (is.null(opt\$init)) {
  init <- "NULL"
}

Argument <- c("seuratObject", "quantile", "kfunc", "nstarts", "nsamples", "seed", "init", "assay")
Value <- c(opt\$seuratObj, opt\$quantile, opt\$kfunc, opt\$nstarts, opt\$nsamples, opt\$seed, init, opt\$assay)
params <- data.frame(Argument, Value)

write.csv(params, paste0(opt\$prefix ,"_params_htodemux.csv"))

# create csv files to save the results from HTODemux()
donors <- rownames(hashtag[[opt\$assay]])
assignment <- hashtag[[paste0(opt\$assay, "_classification")]]
assignment[[paste0(opt\$assay, "_classification")]][!assignment[[paste0(opt\$assay, "_classification")]] %in% c(donors, "Negative")] <- "Doublet"
write.csv(assignment, paste0(opt\$prefix ,"_assignment_htodemux.csv"))
write.csv(hashtag[[paste0(opt\$assay, "_classification.global")]], paste0(opt\$prefix ,"_classification_htodemux.csv"))
<<<<<<< HEAD
saveRDS(hashtag, file = paste0(opt\$prefix ,"htodemux.rds"))
>>>>>>> df971f6e6 (add template)
=======
saveRDS(hashtag, file = paste0(opt\$prefix ,"_htodemux.rds"))
>>>>>>> 37633daa7 (add stub)

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
<<<<<<< HEAD
        paste('    seurat:', seurat.version)
=======
        paste('    bioconductor-deseq2:', seurat.version)
>>>>>>> df971f6e6 (add template)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
