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

#' Converts the .gmt file into a df
#'
#' @param file_gmt_path path of the .gmt file provided by mygene module.
#' @return output dataframe a Dataframe: 1st column = GOterm, 2nd = Description, 3d to end = genes.
process_gmt_file <- function(file_gmt_path) {

    lines <- readLines(file_gmt_path)
    data_list <- list()

    for (line in lines) {
        fields <- strsplit(line, "\\t")[[1]] # Split the line based on the tab character
        go_term <- fields[1] # Extract the GO term

        # Create a data frame with the GO term in the first column
        # Fill in missing values with NA to ensure consistent column lengths
        data_list[[go_term]] <- data.frame(GOterm = go_term,
                                        Description = fields[2],
                                        GeneIDs = c(fields[3:length(fields)], rep(NA, max(0, 3 - length(fields)))))
    }

    gmt_df <- do.call(rbind, data_list) # Combine all data frames into a single data frame
    gmt_df\$GeneIDs <- as.character(gmt_df\$GeneIDs) # Convert gene IDs to character to avoid coercion

    return(gmt_df)
}

#' Converts the .gmt data frame into a knowledge matrix (contingency table)
#'
#' @param gmt_df .gmt df created by process_gmt_file
#' @return output dataframe. A knowledge database where each row is a graph node (gene)
#'  and each column is a concept (GO term).
gmt_to_K<- function(gmt_df){

    summ_df <- as.data.frame(gmt_df\$GeneIDs)
    summ_df <- cbind(summ_df, as.data.frame(gmt_df\$GOterm))
    colnames(summ_df)<- c("GeneIDs", "GOterm")
    summ_df<- unique(summ_df)

    summ_df\$value <- 1

    K <- table(summ_df\$GeneIDs, summ_df\$GOterm)
    K <- as.data.frame.matrix(K)

    return(K)
}

#' Expands knowledge matrix with missing genes to ensure same number of rows for A and K
#'
#' @param adjacency_matrix gene x gene correlation or proportionality adjacency matrix (output propr/propd)
#' @return output dataframe. A knowledge database where each row is a graph node (gene)
#'  and each column is a concept (GO term).
add_missing <- function(adjacency_matrix, knowledge_matrix){

    missing_genes <- setdiff(rownames(adjacency_matrix), rownames(knowledge_matrix))
    extra_rows <- data.frame(matrix(0, nrow = length(missing_genes), ncol = ncol(knowledge_matrix)))
    rownames(extra_rows) <- missing_genes
    colnames(extra_rows) <- colnames(knowledge_matrix)

    knowledge_matrix <- rbind(knowledge_matrix, extra_rows)
    return(knowledge_matrix)
}

################################################
################################################
## Parse arguments                            ##
################################################
################################################

opt <- list(
    adj              = '$adj',
    gmt              = '$gmt',
    prefix           = ifelse('$task.ext.prefix' == 'null', '$meta.id',  '$task.ext.prefix'),
    permutation      = 100,
    fixseed          = TRUE,
    ncores           = as.integer('$task.cpus')
)

opt_types <- list(
    adj              = 'character',
    gmt              = 'character',
    prefix           = 'character',
    permutation      = 'numeric',
    fixseed          = 'logical',
    ncores           = 'numeric'
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

# Read gene x gene adjacency matrix
A <- read_delim_flexible(opt\$adj, header = TRUE, row.names = 1, check.names = TRUE)

# Read and process gene x GO term matrix
gmt_df <- process_gmt_file(opt\$gmt)
K <- gmt_to_K(gmt_df)

# Ensure same number of rows (genes)
if (nrow(A) != nrow(K)){
    K <- add_missing(A, K)
}

# Run Graflex
G <- runGraflex(A, K, opt\$permutation, opt\$fixseed)

################################################
################################################
## Generate outputs                           ##
################################################
################################################

write.table(
    G,
    file      = paste0(opt\$prefix, '.go.tsv'),
    col.names = TRUE,
    row.names = TRUE,
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
