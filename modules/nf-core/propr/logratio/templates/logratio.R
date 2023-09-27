#!/usr/bin/env Rscript


################################################
################################################
## Functions                                  ##
################################################
################################################

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

    read.delim(
        file,
        sep = separator,
        header = header,
        row.names = row.names,
        check.names = check.names
    )
}


#' Check if a variable can be numeric or not
#' 
#' @param x Input variable
#' @retur True if it can be numeric, False otherwise
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
#' 
#' @return The reference gene index
set_reference <- function(ivar, mat){
    if (ivar == 'null'){
        ivar <- ncol(mat)
    }else{
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
    aX <- (ct^alpha - 1) / alpha
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
    count     = '$count',
    prefix    = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    transform = '$transformation',
    reference = '$reference',
    alpha     = ifelse('$alpha' == 'null', NA, as.numeric('$alpha'))
)


if (!opt\$transform %in% c('clr', 'alr')) stop('Please make sure you provided the correct lr_transformation')


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
mat <- read_delim_flexible(opt\$count)
mat <- t(mat)


# check zeros
# log transformation should be applied on non-zero data
# otherwise Inf values are generated
if (any(mat == 0)) stop("There are missing values in the input matrix. Please handle the zeros before running this script")


# compute ALR/CLR
if (opt\$transform == 'alr'){

    # get reference set
    opt\$ivar <- set_reference(opt\$reference, mat)
    gene_name <- colnames(mat)[opt\$ivar]

    # get alr
    if (is.na(opt\$alpha)){
        logratio <- get_logratio(mat, opt\$ivar)
    }else{
        logratio <- get_boxcox(mat, opt\$ivar, opt\$alpha)
    }

    # set the reference column to NA
    # TODO decide if set reference column to NA, or remove column directly
    logratio[,opt\$ivar] = NA

}else if (opt\$transform == 'clr'){

    # get clr
    if (is.na(opt\$alpha)){
        logratio <- get_logratio(mat, 'clr')
    }else{
        logratio <- get_boxcox(mat, opt\$ivar, opt\$alpha)
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