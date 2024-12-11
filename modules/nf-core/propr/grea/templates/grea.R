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

#' Flexibly read CSV or TSV files (determined by file extension)
#'
#' @param file Input file
#' @param header Boolean. TRUE if first row is header. False without header.
#' @param row.names The first column is used as row names by default.
#' Otherwise, give another number. Or use NULL when no row.names are present.
#'
#' @return output Data frame
read_delim_flexible <- function(file, header = TRUE, row.names = 1, check.names = TRUE){

    ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1)) # Get the file extension

    if (ext == "tsv" || ext == "txt") { # If the file is a tsv or txt file
        separator <- "\\t" # Set the separator variable to tab
    } else if (ext == "csv") { # If the file is a csv file
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }

    mat <- read.delim( # Read the file
        file,
        sep = separator, # Set the separator defined above
        header = header,
        row.names = row.names,
        check.names = check.names
    )
}

#' Loads the .gmt file  and converts it into a knowledge database
#'
#' @param filename path of the .gmt file
#' @param nodes vector of node (eg. gene) names. Note that this set should be as
#' complete as possible. So it should not only contain the target genes but also
#' the background genes.
#' @return a list with:
#'     `db` A knowledge database (matrix) where each row is a graph node (eg. gene)
#'      and each column is a concept (eg. GO term, pathway, etc).
#'     `description` A list of descriptions for each concept.
load_gmt <- function(filename, nodes) {

    # read gmt file
    gmt <- readLines(filename)
    gmt <- strsplit(gmt, "\\t")

    # initialize database matrix
    db <- matrix(0, nrow = length(nodes), ncol = length(gmt))
    rownames(db) <- nodes
    colnames(db) <- sapply(gmt, function(entry) entry[[1]])

    # description of the concepts
    description <- list()

    # for concept in gmt
    for (i in 1:length(gmt)) {

        # get concept and description
        concept <- gmt[[i]][[1]]
        description[[concept]] <- gmt[[i]][[2]]

        # fill 1 if gene is in concept
        nodes_in_concept <- gmt[[i]][-c(1, 2)]
        nodes_in_concept <- nodes_in_concept[nodes_in_concept %in% nodes]
        db[nodes_in_concept, i] <- 1
    }

    return(list(db = db, description = description))
}

################################################
################################################
## Parse arguments                            ##
################################################
################################################

# Set defaults and classes

opt <- list(
    prefix           = ifelse('$task.ext.prefix' == 'null', '$meta.id',  '$task.ext.prefix'),

    # input data
    adj              = '$adj',          # adjacency matrix
    gmt              = '$gmt',          # knowledge database .gmt file

    # parameters for gene sets
    set_min          = 15,              # minimum number of genes in a set
    set_max          = 500,             # maximum number of genes in a set

    # parameters for permutation test
    permutation      = 100,             # number of permutations to perform

    # other parameters
    seed             = NA,              # seed for reproducibility
    round_digits     = NA,              # number of digits to round results
    ncores           = as.integer('$task.cpus')
)

opt_types <- list(
    prefix           = 'character',
    adj              = 'character',
    gmt              = 'character',
    set_min          = 'numeric',
    set_max          = 'numeric',
    permutation      = 'numeric',
    seed             = 'numeric',
    round_digits     = 'numeric',
    ncores           = 'numeric'
)

# Apply parameter overrides

args_ext <- ifelse('$task.ext.args' == 'null', '', '$task.ext.args')
args_opt <- parse_args(args_ext)
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    } else {

        # Preserve classes from defaults where possible
        args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])

        # handle NA, and avoid errors when NA is provided by user as character
        if (args_opt[[ao]] %in% c('NA', NA)) args_opt[[ao]] <- NA

        # replace values
        opt[[ao]] <- args_opt[[ao]]
    }
}

# Check if required parameters have been provided

required_opts <- c('adj', 'gmt') # defines a vector required_opts containing the names of the required parameters.
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]
if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid

for (file_input in c('adj', 'gmt')){
    if (is.null(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }
    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

# check parameters are valid

if (opt\$permutation < 0) {
    stop('permutation should be a positive integer')
}

print(opt)

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(propr)

################################################
################################################
## Enrichment analysis                        ##
################################################
################################################

# set seed when required

if (!is.na(opt\$seed)) {
    warning('Setting seed ', opt\$seed, ' for reproducibility')
    set.seed(opt\$seed)
}

# load adjacency matrix
# this matrix should have gene x gene dimensions

adj <- as.matrix(read_delim_flexible(
    opt\$adj,
    header = TRUE,
    row.names = 1,
    check.names = FALSE
))
if (nrow(adj) != ncol(adj)) {
    stop('Adjacency matrix should be a squared matrix that reflects the connections between all the nodes')
}
if (!all(rownames(adj) == colnames(adj))) {
    stop('Adjacency matrix row names are not equal to column names')
}

# load and process knowledge database

gmt <- load_gmt(
    opt\$gmt,
    rownames(adj)  # adj should contain all the nodes (target and background)
)

# filter gene sets
# gene sets with less than set_min or more than set_max genes are removed

idx <- which(colSums(gmt\$db) > opt\$set_min & colSums(gmt\$db) < opt\$set_max)
gmt\$db <- gmt\$db[, idx]
gmt\$description <- gmt\$description[idx]

# run GREA
# Basically, it calculates the odds ratio of the graph being enriched in each concept,
# and the FDR of the odds ratio through permutation tests

odds <- runGraflex(
    adj,
    gmt\$db,
    p=opt\$permutation,
    ncores=opt\$ncores
)
odds\$Description <- sapply(odds\$Concept, function(concept)
    gmt\$description[[concept]]
)

################################################
################################################
## Generate outputs                           ##
################################################
################################################

if (!is.na(opt\$round_digits)) {
    for (col in c('Odds', 'LogOR', 'FDR.under', 'FDR.over')){
        odds[,col] <- round(odds[,col], opt\$round_digits)
    }
}

write.table(
    odds,
    file      = paste0(opt\$prefix, '.grea.tsv'),
    col.names = TRUE,
    row.names = FALSE,
    sep       = '\\t',
    quote     = FALSE

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
