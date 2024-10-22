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

    df[,columns] <- format(
        data.frame(df[, columns], check.names = FALSE),
        nsmall = digits
    )

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

opt <- list(
    querygse = '$querygse',
    metacols = NULL
)
args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }else{
        opt[[ao]] <- args_opt[[ao]]
    }
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(GEOquery)

################################################
################################################
## Do the GEO query retrieval                 ##
################################################
################################################

# Fetch data for GSE number

eset <- getGEO(
    GEO = opt\$querygse,
    destdir = getwd()
)[[1]]

# Write probeset annotation. If supplied, Parse metadata columns from nextflow
# parameters to subset on the feature metadata file

probeset_annotation = fData(eset)
if (! is.null(opt\$metacols)){
    feature_cols = strsplit(opt\$metacols,',')[[1]]
    probeset_annotation <- probeset_annotation[,feature_cols]
}

################################################
################################################
## Generate outputs                           ##
################################################
################################################

output_prefix <- ifelse('$task.ext.prefix' == 'null', '', '$task.ext.prefix')

write.table(
    probeset_annotation,
    paste0(output_prefix,'annotation.tsv'),
    col.names=TRUE,
    row.names=FALSE,
    sep="\t",
    quote=FALSE
)

# If data is not log scale, transform it as needed for limma downstream

if(max(exprs(eset),na.rm=T) > 20) { # a bit dirty, needs proper solution later...
    exprs(eset) <- log2(exprs(eset) + 1)
}

saveRDS(eset, file = paste0(output_prefix, 'eset.rds'))

# Write intensity matrix (normalised)

write.table(
    data.frame(
        probe_id = rownames(eset),
        round_dataframe_columns(as.data.frame(exprs(eset))),
        check.names = FALSE
    ),
    file = paste0(output_prefix, 'matrix.tsv'),
    col.names = TRUE, row.names = FALSE,
    sep = '\t', quote = FALSE
)

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(output_prefix, "R_sessionInfo.log", sep = '.'))
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
geoquery.version <- as.character(packageVersion("GEOquery"))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-geoquery:', geoquery.version)
    ),
    'versions.yml')

################################################
################################################
################################################
################################################
