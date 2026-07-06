#!/usr/bin/env Rscript

#written by Mo Tan (https://github.com/Schansiate) modified from original script by Oliver Wacker (https://github.com/WackerO)

# MIT License

# Copyright (c) QBiC

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

################################################
################################################
## Functions                                  ##
################################################
################################################

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
## Load libraries                             ##
################################################
################################################

library(gprofiler2)
library(nfcore.utils)

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
    output_prefix = '$task.ext.prefix'
)

opt <- process_inputs(
    opt,
    args = '$task.ext.args',
    keys_to_nullify = c('output_prefix'),
    expected_files = c('ids_tsv'),
    expected_double = c('mthreshold'),
    expected_boolean = c('filter_na'),
    required_opts = c('ids_tsv', 'target')
)

if (is.null(opt[['output_prefix']])) {
    opt[['output_prefix']] <- '$meta.id'
}

################################################
################################################
## READ IDS AND RUN GCONVERT                  ##
################################################
################################################

ids_table <- read_delim_flexible(
    file = opt[['ids_tsv']]
)

if (ncol(ids_table) < 1) {
    stop("The input ID table must contain at least one column.")
}

output_prefix <- paste0(opt[['output_prefix']], ".gprofiler2")
output_file <- paste(output_prefix, 'gconvert', 'tsv', sep = '.')

query <- as.character(ids_table[[1]])
query <- query[!is.na(query) & query != ""]

if (length(query) > 0) {
    gconvert_results <- gconvert(
        query = query,
        organism = opt[['organism']],
        target = opt[['target']],
        numeric_ns = opt[['numeric_ns']],
        mthreshold = opt[['mthreshold']],
        filter_na = opt[['filter_na']]
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
## R SESSION INFO AND VERSIONS FILE           ##
################################################
################################################

create_log_session_info()

create_versions_yml(
    packages = list("r-gprofiler2" = "gprofiler2"),
    task_name = '$task.process',
    versions_path = "versions.yml"
)

gprofiler_data.version <- gprofiler2::get_version_info(opt[["organism"]])[["gprofiler_version"]]

write(
    paste0('    gprofiler-data: ', gprofiler_data.version),
    file = "versions.yml",
    append = TRUE
)

################################################
################################################
################################################
################################################
