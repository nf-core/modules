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

    if (!is.null(row.names)){
        if ( (row.names == 'gene_id') & ('gene_name' %in% colnames(mat)) ){
            mat <- mat[, -which(colnames(mat) == 'gene_name')]
        } else if ( (row.names == 'gene_name') & ('gene_id' %in% colnames(mat)) ){
            mat <- mat[, -which(colnames(mat) == 'gene_id')]
        }
    }

    return(mat)
}


################################################
################################################
## Parse arguments                            ##
################################################
################################################

opt <- list(
    prefix          = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    count           = '$count',
    samplesheet     = '$samplesheet',
    features_id_col = 'gene_id',
    obs_id_col      = 'sample',
    group_col       = 'treatment',
    metric          = 'theta_d',
    alpha           = NA,
    permutation     = 0,
    cutoff_min      = NA,
    cutoff_max      = NA,
    cutoff_interval = NA,
    ncores          = as.integer('$task.cpus')
)
opt_types <- list(
    prefix          = 'character',
    count           = 'character',
    samplesheet     = 'character',
    features_id_col = 'character',
    obs_id_col      = 'character',
    group_col       = 'character',
    metric          = 'character',
    alpha           = 'numeric',
    permutation     = 'numeric',
    cutoff_min      = 'numeric',
    cutoff_max      = 'numeric',
    cutoff_interval = 'numeric',
    ncores          = 'numeric'
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
required_opts <- c('count','samplesheet')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]
if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid
for (file_input in c('count','samplesheet')){
    if (is.null(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }
    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

# check parameters
if (!opt\$metric %in% c('theta_d', 'theta_e', 'theta_f')) stop('Please provide a valid differential proportionality metric')


################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(propr)


################################################
################################################
## Perform differential proportionality       ##
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
if (any(mat == 0)) print("Zeros will be replaced by minimun value before logratio analysis")

# parse group
# this creates a vector referring to the group id for each observation
samplesheet <- read_delim_flexible(
    opt\$samplesheet,
    header = TRUE,
    row.names = NULL,
    check.names = FALSE
)
tmp <- samplesheet[[opt\$group_col]]
names(tmp) <- samplesheet[[opt\$obs_id_col]]
group <- as.vector(tmp[rownames(mat)])
if (length(group) != nrow(mat)) stop('Error when parsing group')

# perform differential proportionality
pd <- propd(
    mat,
    group    = group,
    alpha    = opt\$alpha,
    weighted = FALSE,
    p        = opt\$permutation
)
if (opt\$metric == 'theta_d'){
    pd <- setDisjointed(pd)
}else if (opt\$metric == 'theta_e'){
    pd <- setEmergent(pd)
}else if (opt\$metric == 'theta_f'){
    pd <- setActive(pd, what = "theta_f")
}

# update FDR by permutation, if required
if (opt\$permutation > 0) {
    cutoff <- seq(
        opt\$cutoff_min,
        opt\$cutoff_max,
        opt\$cutoff_interval
    )
    pd <- updateCutoffs(pd, cutoff=cutoff, ncores=opt\$ncores)
}


################################################
################################################
## Generate outputs                           ##
################################################
################################################

write.table(
    getResults(pd),
    file      = paste0(opt\$prefix, '.propd.tsv'),
    col.names = TRUE,
    row.names = FALSE,
    sep       = '\t',
    quote     = FALSE
)

if (opt\$permutation > 0) {
    write.table(
        pd@fdr,
        file      = paste0(opt\$prefix, '.fdr.tsv'),
        col.names = TRUE,
        sep       = '\t',
        quote     = FALSE
    )
}


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
