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
#' @param header Boolean. TRUE if first row is header. False without header.
#' @param row.names The first column is used as row names by default.
#' Otherwise, give another number. Or use NULL when no row.names are present.
#'
#' @return output Data frame
read_delim_flexible <- function(file, header = TRUE, row.names = 1, check.names = TRUE){

    ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1))

    if (ext == "tsv" || ext == "txt") {
        separator <- "\\t"
    } else if (ext == "csv") {
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }

    mat <- read.delim(
        file,
        sep = separator,
        header = header,
        row.names = row.names,
        check.names = check.names
    )

    if ( (row.names == 'gene_id') & ('gene_name' %in% colnames(mat)) ){
        mat <- mat[, -which(colnames(mat) == 'gene_name')]
    } else if ( (row.names == 'gene_name') & ('gene_id' %in% colnames(mat)) ){
        mat <- mat[, -which(colnames(mat) == 'gene_id')]
    }

    return(mat)
}

#' Check if a variable can be numeric or not
#'
#' @param x Input variable
#' @return True if it can be numeric, False otherwise
can_be_numeric <- function(x) {
    stopifnot(is.atomic(x) || is.list(x)) # check if x is a vector
    numNAs <- sum(is.na(x))
    numNAs_new <- suppressWarnings(sum(is.na(as.numeric(x))))
    return(numNAs_new == numNAs)
}

#' Set the proper reference gene index.
#' This should be used for alr transformation only.
#'
#' @param ivar Reference variable given by user.
#' If it is 'null', then set the last column as reference (default).
#' Otherwise, it should refer to either gene name or gene index.
#' If the gene name is given, find its index.
#' @param mat Data matrix, with genes as columns
#'
#' @return The reference gene index
set_reference <- function(ivar, mat){
    if (is.na(ivar)){
        ivar <- ncol(mat)
    } else {
        isnumeric <- can_be_numeric(ivar)
        if (!isnumeric){
            genes <- colnames(mat)
            ivar <- which(genes == ivar)
        }
        ivar <- as.integer(ivar)
    }
    return(ivar)
}

#' Set the appropriate range for the sequence of cutoffs used in updateCutoffs.
#' Adjusts the interval to the different metrics.
#'
#' @param object propr object. Output from propr function.
#'
#' @return sequence of cutoff values.
seqCutoff <- function(object){
    matrix <- getMatrix(object)
    matrix[matrix <= 0] <- NA
    diag(matrix) <- NA
    min_cutoff <- round(min(matrix, na.rm = TRUE),3)
    max_cutoff <- round(max(matrix, na.rm = TRUE),3)
    step_cutoff <- (max_cutoff - min_cutoff)/ 20
    seq_cutoff <- seq(min_cutoff, max_cutoff, step_cutoff)
    return(seq_cutoff)
}

#' Extract the proportionality cutoff for a specified FDR value.
#' Gene pairs with a proportionality value higher than the extracted cutoff will be considered significantly proportional.
#'
#' @param object propr object. Output from propr function. updateCutoffs function should be applied to the object previous to valCutoff.
#' @param fdrVal FDR value to extract the cutoff for. Per default 0.05
#' @param metric Metric used to calculate the proportionality values. Options are 'cor', 'rho', 'phi', 'phs', 'vlr', 'pcor', 'pcor.shrink', 'pcor.bshrink'
#'
#' @return cutoff value. Proportionality values higher than this cutoff are considered significant.
valCutoff  <- function(object, metric, fdrVal = 0.05){
    fdr_df <- object@fdr
    print(fdr_df)
    # metric_up <- c("rho", "cor", "pcor", "pcor.shrink", "pcor.bshrink")

    if (prod(dim(fdr_df) == 0)){
        warning("Please run updateCutoff on propr first")
    }else{
        fdr_vals <- fdr_df\$FDR
        if(any(!is.na(fdr_vals))){ # if there is some defined value, continue, else out of range
            if(any(fdr_vals <= fdrVal)){ # if there is some value that is belowe the FDR threshold,
                fdr_threshold <- fdr_vals[which.max(fdr_vals <= fdrVal)] #choose the highest FDR that is lower than the threshold, else choose the lowest
            }else{
                warning("FDR is higher than the specified threshold for all proportionality values. Using the lowest fdr instead")
                fdr_threshold <- min(fdr_vals, na.rm = TRUE)
            }
            cutoff <- fdr_df\$cutoff[which(fdr_df\$FDR == fdr_threshold)] #select the corresponding cutoff value for the FDR
            print(cutoff)
        }else{
            stop("FDR not defined. This metric is not appropriate for the given dataset")
        }
    return(cutoff)
    }
}


#' Convert a proportionality matrix to an adjacency matrix based on a threshold.
#'
#' @param matrix proportionality matrix. Can be extracted from propr object with getMatrix().
#' @param cutoff Significant proportionality value extracted from valCutoff function.
#'
#' @return Adjacency matrix. Gene pairs with a proportionality value higher than the threshold will have 1, otherwise 0.
convert_to_adjacency <- function(matrix, cutoff, metric) {
    if (metric == 'cor' || metric == 'rho' || metric == 'pcor' || metric == 'pcor.shrink' || metric == 'pcor.bshrink'){
        adjacency <- ifelse(matrix > cutoff, 1, 0)
    } else {
        adjacency <- ifelse(matrix < cutoff, 1, 0)
    }
    return(adjacency)
}

################################################
################################################
## Parse arguments                            ##
################################################
################################################

opt <- list(
    count            = '$count',
    prefix           = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    transformation   = 'clr',
    reference        = NA,
    alpha            = NA,
    metric           = 'pcor.bshrink',
    permutation      = 0,
    cutoff_min       = NA,
    cutoff_max       = NA,
    cutoff_interval  = NA,
    ncores           = as.integer('$task.cpus'),
    features_id_col  = 'gene_id',
    fixseed          = FALSE,
    adjacency        = FALSE,
    fdrVal           = 0.05
)
opt_types <- list(
    count            = 'character',
    prefix           = 'character',
    transformation   = 'character',
    reference        = 'character',
    alpha            = 'numeric',
    metric           = 'character',
    permutation      = 'numeric',
    cutoff_min       = 'numeric',
    cutoff_max       = 'numeric',
    cutoff_interval  = 'numeric',
    ncores           = 'numeric',
    features_id_col  = 'character',
    fixseed          = 'logical',
    adjacency        = 'logical',
    fdrVal           = 'numeric'
)

# Apply parameter overrides
args_opt <- parse_args('$task.ext.args')

for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    } else {

        # Preserve classes from defaults where possible
        if (! is.null(opt[[ao]])){
            args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        }
        # set NA
        if (args_opt[[ao]] %in% c('NA', NA, 'null')){
            args_opt[[ao]] <- NA
        }
        opt[[ao]] <- args_opt[[ao]]
    }
}

# Check if required parameters have been provided
required_opts <- c('count')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]
if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid
for (file_input in c('count')){
    if (is.null(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }
    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

# check parameters
if (!opt\$metric %in% c('rho', 'phi', 'phs', 'cor', 'vlr', 'pcor', 'pcor.shrink', 'pcor.bshrink')) {
    stop('Please make sure you provided the correct metric')
}
if (opt\$metric == 'pcor.bshrink'){
    if (!is.na(opt\$alpha)) stop('Box-cox transformation is not implemented for pcor.bshrink yet.')
    if (!opt\$transformation %in% c('clr', 'alr')) stop('Please make sure you provided the correct transformation: clr or alr')
} else {
    if (!opt\$transformation %in% c('clr', 'alr', NA)) stop('Please make sure you provided the correct transformation: clr or alr. Or set NA if you dont want to transform the data.')
    if (is.na(opt\$transformation)) print('Warning: No transformation is required by user. We assume the input count data was already properly transformed.')
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(corpcor)
library(propr)

################################################
################################################
## Perform correlation analysis               ##
################################################
################################################

# read matrix
mat <- read_delim_flexible(
    opt\$count,
    header = TRUE,
    row.names = opt\$features_id_col,
    check.names = FALSE
)
mat <- t(mat)

# check zeros
# log transformation should be applied on non-zero data
# otherwise Inf values are generated
if (any(mat == 0)) print("Warning: Zeros will be replaced by minimum value before logratio analysis")

# set logratio transformation parameter -> ivar
# if alr, set the index of the reference gene as ivar
# if clr or NA, set same ivar
if (opt\$metric == 'pcor.bshrink'){
    opt\$ivar <- opt\$transformation
}else{
    if (opt\$transformation == 'alr'){
        opt\$ivar <- set_reference(opt\$reference, mat)
        gene_name <- colnames(mat)[opt\$ivar]
    } else {
        opt\$ivar <- opt\$transformation
    }
}

# Compute correlation coefficients
pro <- propr(
    mat,
    metric  = opt\$metric,
    ivar    = opt\$ivar,
    alpha   = opt\$alpha,
    p       = opt\$permutation,
    fixseed = opt\$fixseed
)

# update FDR by permutation, if required

if (opt\$permutation > 0) {
    cutoff <- seq(
        opt\$cutoff_min,
        opt\$cutoff_max,
        opt\$cutoff_interval
        )
    if (is.na(opt\$cutoff_min) ||  is.na(opt\$cutoff_max) || is.na(opt\$cutoff_interval)) {
        warning("cutoff values were not provided. Using the default cutoff values.")
        cutoff <- seqCutoff(pro)
    }
    m <- getMatrix(pro)
    diag(m) <- NA
    print((opt\$cutoff_max - opt\$cutoff_min)/2 + opt\$cutoff_min)
    print(max(m, na.rm = TRUE))
    if ((opt\$cutoff_max - opt\$cutoff_min)/2 + opt\$cutoff_min > max(m, na.rm = TRUE)) {
        warning("The provided cutoff values are out of range. Using the default cutoff values.")
        cutoff <- seqCutoff(pro)
    }
    pro <- updateCutoffs(pro, cutoff=cutoff, ncores=opt\$ncores)
}

# calculate cutoff and adjacency matrix, if required

if (opt\$adjacency == TRUE) {
    cutoff <- valCutoff(pro, opt\$metric, opt\$fdrVal)
    matrix <- getMatrix(pro)
    adj <- convert_to_adjacency(matrix, cutoff, opt\$metric)
}

################################################
################################################
## Generate outputs                           ##
################################################
################################################

saveRDS(
    pro,
    file = paste0(opt\$prefix, '.propr.rds')
)

write.table(
    round(pro@matrix, 8),  # round matrix decimals to avoid floating point inconsistencies
    file      = paste0(opt\$prefix, '.propr.tsv'),
    col.names = TRUE,
    row.names = TRUE,
    sep       = '\t',
    quote     = FALSE
)

if (opt\$permutation > 0) {
    write.table(
        pro@fdr,
        file      = paste0(opt\$prefix, '.fdr.tsv'),
        col.names = TRUE,
        row.names = FALSE,
        sep       = '\t',
        quote     = FALSE
    )
}

if (opt\$adjacency == TRUE) {
    write.table(
        adj,
        file      = paste0(opt\$prefix, '.adj.csv'),
        col.names = TRUE,
        row.names = TRUE,
        sep       = ',',
        quote     = FALSE
    )
}

################################################
################################################
## WARNINGS                                   ##
################################################
################################################

sink(paste0(opt\$prefix, ".warnings.log"))
print(warnings())
sink()

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste0(opt\$prefix, ".R_sessionInfo.log"))
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

propr.version <- as.character(packageVersion('propr'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-propr:', propr.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
