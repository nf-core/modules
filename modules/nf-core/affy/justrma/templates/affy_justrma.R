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

#' Install the right CDF for a given cel file
#'
#' @param celfile A valid path to a CEL file
#' @param annotation Boolean indication wheter to install the annotation
#'        package
#'
#' @return output The CDF environment or a list detailing the failed locations.

install_cdf_db <- function(celfile, annotation = FALSE){
    library(affyio)
    headdetails <- read.celfile.header(celfile)
    ref.cdfName <- headdetails[[1]]
    cleaned.cdfName <- cleancdfname(ref.cdfName, addcdf = FALSE)

    exts = 'cdf'
    if (annotation){
        exts <- c(exts, '.db')
    }
    options(timeout=600)
    for (package in paste0(cleaned.cdfName, exts)){
        install.packages(
            package,
            lib = 'libs',
            repos = BiocManager::repositories(),
            dependencies = c("Depends", "Imports")
        )
    }
    cleaned.cdfName
}

#' Round numeric dataframe columns to fixed decimal places by applying
#' formatting and converting back to numerics
#'
#' @param dataframe A data frame
#' @param columns Which columns to round (assumes all of them by default)
#' @param digits How many decimal places to round to?
#'
#' @return output Data frame

round_dataframe_columns <- function(df, columns = NULL, digits = 8){
    if (is.null(columns)){
        columns <- colnames(df)
    }

    df[,columns] <- format(data.frame(df[, columns]), nsmall = digits)

    # Convert columns back to numeric

    for (c in columns) {
        df[[c]][grep("^ *NA\$", df[[c]])] <- NA
        df[[c]] <- as.numeric(df[[c]])
    }
    df
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
    bgversion = 2,
    destructive=FALSE,
    cdfname = NULL,
    rm.mask = FALSE,
    rm.outliers = FALSE,
    rm.extra = FALSE,
    build_annotation = FALSE,
    keep.log2 = TRUE
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

# Install the CDF in a path we can definitely write to (for some reason affy
# ignores .libPaths() left to its own devices)

dir.create('libs')
.libPaths('libs')
first_cel <- file.path(opt\$celfiles_dir, sample.sheet[[opt\$file_name_col]][1])
cdf_name <- install_cdf_db(first_cel, annotation = opt\$build_annotation)

# Run the main function

rownames(sample.sheet) <- sample.sheet[[opt\$file_name_col]]
eset <- justRMA(
    filenames = sample.sheet[[opt\$file_name_col]],
    celfile.path = opt\$celfiles_dir,
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

# This should happen as part of the above, not sure why it doesn't
sampleNames(eset) <- sample.sheet[[opt\$sample_name_col]]

# unlog2 the data, if keep.log2 is set to FALSE
if (!opt\$keep.log2) exprs(eset) <- 2**exprs(eset)

################################################
################################################
## Generate outputs                           ##
################################################
################################################

if (opt\$build_annotation){

    # Make some annotation

    dbname <- paste0(cdf_name, '.db')
    library(dbname, character.only = TRUE)
    anno <- select(
        get(dbname),
        keys=keys(get(dbname), keytype="PROBEID"),
        columns=c('ENSEMBL', 'ENTREZID', 'SYMBOL', 'GENENAME', 'GENETYPE'),
        keytype="PROBEID"
    )

    # Remove duplicates by probe
    anno <- do.call(
        rbind,
        lapply(
            split(
                anno,
                anno\$PROBEID
            ),
            function(x) apply(x, 2, function(y) paste(unique(y), collapse=', '))
        )
    )

    write.table(
        anno,
        file = paste0(cdf_name, '.annotation.tsv'),
        col.names = TRUE,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )
}

# R object for other processes to use

output_prefix <- '$prefix'
saveRDS(eset, file = paste0(output_prefix, '_eset.rds'))

# Write matrix

write.table(
    data.frame(
        probe_id = rownames(eset),
        round_dataframe_columns(as.data.frame(exprs(eset))),
        check.names = FALSE
    ),
    file = paste0(output_prefix, '_matrix.tsv'),
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
