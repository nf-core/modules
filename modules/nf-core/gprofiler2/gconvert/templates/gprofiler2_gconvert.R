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

parse_args <- function(x) {
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z) { length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = FALSE) {

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

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

opt <- list(
    ids_tsv = '$ids_tsv',
    target = '$target',
    organism = 'hsapiens',
    numeric_ns = '',
    mthreshold = Inf,
    filter_na = TRUE,
    output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
)

opt_types <- lapply(opt, class)

args_opt <- parse_args('$task.ext.args')
allowed_args <- c('organism', 'numeric_ns', 'mthreshold', 'filter_na')
for (ao in names(args_opt)) {
    if (! ao %in% allowed_args) {
        stop(paste("Invalid option:", ao))
    } else {
        if (! is.null(opt[[ao]])) {
            args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        }
        opt[[ao]] <- args_opt[[ao]]
    }
}

required_opts <- c('ids_tsv', 'target', 'output_prefix')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]

if (length(missing) > 0) {
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

if (opt\$target == "") {
    stop("Please provide target.")
}

if (! file.exists(opt\$ids_tsv)) {
    stop(paste0('Value of ids_tsv: ', opt\$ids_tsv, ' is not a valid file'))
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(gprofiler2)

################################################
################################################
## READ IDS AND RUN GCONVERT                  ##
################################################
################################################

ids_table <- read_delim_flexible(
    file = opt\$ids_tsv
)

if (ncol(ids_table) < 1) {
    stop("The input ID table must contain at least one column.")
}

output_prefix <- paste0(opt\$output_prefix, ".gprofiler2")
output_file <- paste(output_prefix, 'gconvert', 'tsv', sep = '.')

query <- as.character(ids_table[[1]])
query <- query[!is.na(query) & query != ""]

if (length(query) > 0) {
    gconvert_results <- gconvert(
        query = query,
        organism = opt\$organism,
        target = opt\$target,
        numeric_ns = opt\$numeric_ns,
        mthreshold = opt\$mthreshold,
        filter_na = opt\$filter_na
    )

    write.table(
        gconvert_results,
        file = output_file,
        col.names = TRUE,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )
} else {
    print("No input IDs found, gProfiler2 g:Convert will be skipped.")
    file.create(output_file)
}

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
gprofiler2.version <- as.character(packageVersion('gprofiler2'))
gprofiler_data.version <- gprofiler2::get_version_info(opt[["organism"]])[["gprofiler_version"]]

writeLines(
    c(
        '"$task.process":',
        paste0('    r-base: ', r.version),
        paste0('    r-gprofiler2: ', gprofiler2.version),
        paste0('    gprofiler-data: ', gprofiler_data.version)
    ),
    con = 'versions.yml'
)

################################################
################################################
################################################
################################################
