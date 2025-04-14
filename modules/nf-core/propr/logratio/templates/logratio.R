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


#' Set the proper reference gene index
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

#' Compute logratio
#'
#' @param mat Count data, with rows = samples and columns = genes
#' @param ivar 'clr' or gene index number for 'alr' transformation
#'
#' @return Logratio matrix
get_logratio <- function(mat, ivar){
    use <- propr:::ivar2index(mat, ivar)
    logX <- log(mat)
    logSet <- logX[, use, drop = FALSE]
    ref <- rowMeans(logSet)
    lr <- logX - ref
    return(lr)
}

#' Compute BoxCox-transformed ratio
#'
#' Since y = log(x) = [x^a-1]/a, ref = Mean[log(x)] = Mean[y]
#' Calculate log(x/ref) = log(x) - log(ref) = [x^a-1]/a - [ref^a-1]/a
#'
#' @param mat Count data, with rows = samples and columns = genes
#' @param ivar 'clr' or gene index number for 'alr' transformation
#' @param alpha parameter for BoxCox transformation
#'
#' @return BoxCox-transformed ratio matrix
get_boxcox <- function(mat, ivar, alpha){
    use <- propr:::ivar2index(mat, ivar)
    aX <- (mat^alpha - 1) / alpha
    aSet <- aX[, use, drop = FALSE]
    ref <- rowMeans(aSet)
    lr <- sweep(aX, 1, ref, "-")
    return(lr)
}


################################################
################################################
## Parse arguments                            ##
################################################
################################################

opt <- list(
    count          = '$count',
    prefix         = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    transformation = 'clr',       # it can be clr or alr
    reference      = NA,          # reference gene (name or index) for alr
    alpha          = NA,          # alpha for box-cox transformation
    feature_id_col = 'gene_id'    # column name of the feature ids
)
opt_types <- lapply(opt, class)
opt_types\$reference <- 'character'
opt_types\$alpha <- 'numeric'


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


if (!opt\$transformation %in% c('clr', 'alr')) stop('Please make sure you provided the correct lr_transformation')


################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(propr)
library(MASS)


################################################
################################################
## Compute ALR transformation                 ##
################################################
################################################


# read matrix
mat <- read_delim_flexible(
    opt\$count,
    header = TRUE,
    row.names = opt\$feature_id_col,
    check.names = FALSE
)
mat <- t(mat)


# check zeros
# log transformation should be applied on non-zero data
# otherwise Inf values are generated
if (any(mat == 0)) print("Zeros will be replaced by minimum value before logratio analysis")


# compute ALR/CLR
if (opt\$transformation == 'alr'){

    # get reference set
    opt\$ivar <- set_reference(opt\$reference, mat)
    gene_name <- colnames(mat)[opt\$ivar]
    opt\$reference <- gene_name

    # get alr
    if (is.na(opt\$alpha)){
        logratio <- get_logratio(mat, opt\$ivar)
    } else {
        logratio <- get_boxcox(mat, opt\$ivar, opt\$alpha)
    }

} else if (opt\$transformation == 'clr'){

    # get clr
    if (is.na(opt\$alpha)){
        logratio <- get_logratio(mat, 'clr')
    } else {
        logratio <- get_boxcox(mat, NA, opt\$alpha)
    }
}


################################################
################################################
## Generate outputs                           ##
################################################
################################################


write.table(
    t(logratio),
    file = paste0(opt\$prefix, '.logratio.tsv'),
    col.names = TRUE,
    row.names = TRUE,
    sep = '\t',
    quote = FALSE
)


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

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
propr.version <- as.character(packageVersion('propr'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    r-propr:', propr.version)
    ),
'versions.yml')


################################################
################################################
################################################
################################################
