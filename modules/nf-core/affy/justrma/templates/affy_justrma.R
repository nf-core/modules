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
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = NULL){

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
        row.names = row.names
    )
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
    sample_file = '$samplesheet',
    celfiles_dir = '$celfiles_dir',
    description = '$description',
    file_name_col = 'file',
    sample_name_col = 'file',
    background = TRUE,
    normalize = TRUE,
    bg_version = 2,
    destructive=FALSE,
    cdfname = NULL,
    rm.mask = FALSE,
    rm.outliers = FALSE,
    rm.extra = FALSE
)
if (opt\$description == ''){
    opt\$description = NULL
}
opt_types <- lapply(opt, class)

# Apply parameter overrides

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

# Check file inputs are valid

for (file_input in c('sample_file')){
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

library(affy)

################################################
################################################
## READ IN SAMPLE SHEET WITH CELL FILE NAMES  ##
################################################
################################################

sample.sheet <- read_delim_flexible(file = opt\$sample_file)

if (! opt\$file_name_col %in% colnames(sample.sheet)){
    stop(paste0("Specified file name column '", opt\$file_name_col, "' is not in the sample sheet"))
}

################################################
################################################
## Run justRMA process                        ##
################################################
################################################

rownames(sample.sheet) <- sample.sheet[[opt\$file_name_col]]
eset <- justRMA(
  filenames = sample.sheet[[opt\$file_name_col]],
  celfile.path = opt\$celfiles_dir,,
  phenoData = sample.sheet,
  description = opt\$description, 
  rm.mask = opt\$rm.mask,
  rm.outliers = opt\$rm.outliers,
  rm.extra = opt\$rm.extra,
  sampleNames = sample.sheet[[opt\$sample_name_col]],
  normalize = opt\$normalize,
  background = opt\$background,
  bgversion = opt\$bgversion,
  destructive = opt\$destructive,
  cdfname = opt\$cdfname
)

################################################
################################################
## Generate outputs                           ##
################################################
################################################

# R object for other processes to use

saveRDS(eset, file = 'eset.rds')

# Write matrix

write.table(
    exprs(eset),
    file = 'matrix.tsv',
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

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
affy.version <- as.character(packageVersion('affy'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-affy:', affy.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
