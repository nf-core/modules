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

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = F) {

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

#' Round numeric dataframe columns to fixed decimal places by applying
#' formatting and converting back to numerics
#'
#' @param dataframe A data frame
#' @param columns Which columns to round (assumes all of them by default)
#' @param digits How many decimal places to round to?
#'
#' @return output Data frame
# TODO check if this is necessary
round_dataframe_columns <- function(df, columns = NULL, digits = 8) {
    if (is.null(columns)) {
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

# I've defined these in a single array like this so that we could go back to an
# optparse-driven method in future with module bin/ directories, rather than
# the template

# Set defaults and classes

opt <- list(
    intensities_file = '$intensities',
    sample_file = '$samplesheet',
    contrast_variable = '$contrast_variable',
    protein_id_col = 'Majority protein IDs',
    sample_id_col = 'sample',
    measure_col_prefix = 'intensities',
    normfuns = 'normalizeMedian',
    plotSampleDistributions_method = 'violin',
    plotMV_loess = T,
    palette_name = 'Set1'
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)) {
    if (! ao %in% names(opt)) {
        stop(paste("Invalid option:", ao))
    } else {

        # Preserve classes from defaults where possible
        if (! is.null(opt[[ao]])) {
            args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        }
        opt[[ao]] <- args_opt[[ao]]
    }
}

# Check if required parameters have been provided

required_opts <- c('intensities_file', 'sample_file', 'contrast_variable')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]

if (length(missing) > 0) {
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid

for (file_input in c('intensities_file', 'sample_file')) {
    if (is.null(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }

    if (! file.exists(opt[[file_input]])) {
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(limma)
library(plotly)
library(proteus)

################################################
################################################
# READ IN INTENSITIES FILE AND SAMPLE METADATA #
################################################
################################################

intensities.table <-
    read_delim_flexible(
        file = opt\$intensities_file,
        check.names = FALSE
    )

sample.sheet <-
    read_delim_flexible(
        file = opt\$sample_file,
        check.names=FALSE
    )

if (! opt\$protein_id_col %in% colnames(intensities.table)) {
    stop(paste0("Specified protein ID column '", opt\$protein_id_col, "' is not in the intensities table"))
}

if (! opt\$sample_id_col %in% colnames(sample.sheet)) {
    stop(paste0("Specified sample ID column '", opt\$sample_id_col, "' is not in the sample sheet"))
}

# Add metadata columns that are necessary for proteus

sample.sheet\$sample <- sample.sheet[[opt\$sample_id_col]]
sample.sheet\$condition <- sample.sheet[[opt\$contrast_variable]]

# Add prefix for proteinGroups measurement columns to the sample IDs from the sampesheet
measure.cols <- setNames(paste0(opt\$measure_col_prefix, sample.sheet[[opt\$sample_id_col]]), sample.sheet[[opt\$sample_id_col]])

# Check that all samples specified in the input sheet are present in the intensities table

missing_columns <- paste0(opt\$measure_col_prefix, sample.sheet[[opt\$sample_id_col]])
missing_columns <- missing_columns[!missing_columns %in% colnames(intensities.table)]
if (length(missing_columns) > 0) {
    stop(paste(
        length(missing_columns),
        'specified samples do not have a(n)',
        opt\$measure_col_prefix,
        'column in intensities table. The following columns are missing:',
        paste(missing_columns, collapse = ', ')
    ))
}

################################################
################################################
## CHECK AND FORMAT NORMFUN AND FILTERFUN     ##
################################################
################################################

valid_normfuns <- c("normalizeMedian", "normalizeQuantiles")
normfuns <- opt\$normfuns

# Check validity of normfun(s)
invalid_normfuns <- normfuns[!(normfuns %in% valid_normfuns)]
if (length(invalid_normfuns)>0) {
    stop(paste0("Invalid normfuns argument(s): ",
        paste(invalid_normfuns, collapse=", "),
        ". Valid normfuns are: ",
        paste(valid_normfuns, collapse=", "),
        "."))
}

################################################
################################################
## Run Proteus processes and generate outputs ##
################################################
################################################

output_prefix <- opt\$contrast_variable

# Replace proteus default ID column with user param and re-set the names of the resulting object (gsub sets the names to NULL)

proteinColumns <- setNames(gsub("Majority protein IDs", opt\$protein_id_col, proteus::proteinColumns), names(proteus::proteinColumns))
proteinGroups <- readProteinGroups(
    file=opt\$intensities_file,
    meta=sample.sheet,
    measure.cols=measure.cols,
    data.cols=proteinColumns
)

# Generate plots for all requested normalizations; also, save normalized protein groups for limma

for (normfun in normfuns) {
    proteinGroups.normalized <- normalizeData(proteinGroups, norm.fun = eval(parse(text=normfun))) # Proteus also accepts other norm.funs, e.g. from limma

    # Apply log2 and remove NAs as these will otherwise mess with some of the following modules

    proteinGroups.normalized\$tab <- na.omit(log2(proteinGroups.normalized\$tab))
    
    png(paste0(output_prefix, '.proteus.', normfun, '_normalised_distributions.png'), width = 5*300, height = 5*300, res = 300, pointsize = 8) 
    print(
        plotSampleDistributions(proteinGroups.normalized, title=paste0("Sample distributions after applying\n", normfun), fill="condition", method=opt\$plotSampleDistributions_method)
         + scale_fill_brewer(palette=opt\$palette_name, name=opt\$contrast_variable)
         + theme(plot.title = element_text(size = 12)) 
        )
    dev.off()
    
    png(paste0(output_prefix, '.proteus.', normfun, '_normalised_mean_variance_relationship.png'), width = 5*300, height = 5*300, res = 300, pointsize = 8) 
    print(
        plotMV(proteinGroups.normalized, with.loess=opt\$plotMV_loess) 
         + ggtitle(paste0("Sample mean variance relationship after applying\n", normfun))
         + scale_fill_distiller(palette=opt\$palette_name)
         + theme(plot.title = element_text(size = 12)) 
        )
    dev.off()

    png(paste0(output_prefix, '.proteus.', normfun, '_normalised_dendrogram.png'), width = 5*300, height = 5*300, res = 300, pointsize = 8)
    print(
        plotClustering(proteinGroups.normalized)
         + ggtitle(paste0("Sample clustering after applying\n", normfun))
         + theme(plot.title = element_text(size = 12))
        )
    dev.off()
    
    # R object for other processes to use
    
    saveRDS(proteinGroups.normalized, file = paste0(output_prefix, '.proteus.', normfun, 'normalised_proteingroups.rds'))

    # Write normalized intensities matrix
    
    out_df <- data.frame(
        proteinGroups.normalized\$tab,
        check.names = FALSE
    )
    out_df[[opt\$protein_id_col]] <- rownames(proteinGroups.normalized\$tab) # proteus saves the IDs as rownames; make column from those
    out_df <- out_df[c(opt\$protein_id_col, colnames(out_df)[colnames(out_df) != opt\$protein_id_col])] # move ID column to first position
    write.table(
        out_df,
        file = paste(output_prefix, 'proteus', normfun, 'normalised_proteingroups_tab', 'tsv', sep = '.'),
        col.names = TRUE,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )
}

# Process and save raw table

proteinGroups\$tab <- na.omit(log2(proteinGroups\$tab))

# Generate raw distribution plot

png(paste0(output_prefix, '.proteus.raw_distributions.png'), width = 5*300, height = 5*300, res = 300, pointsize = 8) 
print(
    plotSampleDistributions(proteinGroups, title="Raw sample distributions", fill="condition", method=opt\$plotSampleDistributions_method)
        + scale_fill_brewer(palette=opt\$palette_name, name=opt\$contrast_variable)
        + theme(plot.title = element_text(size = 12)) 
    )
dev.off()

# R object for other processes to use

saveRDS(proteinGroups, file = paste0(output_prefix, '.proteus.raw_proteingroups.rds'))

# Write raw intensities matrix

out_df <- data.frame(
        proteinGroups\$tab,
        check.names = FALSE
    )
out_df[[opt\$protein_id_col]] <- rownames(proteinGroups\$tab) # proteus saves the IDs as rownames; make column from those
out_df <- out_df[c(opt\$protein_id_col, colnames(out_df)[colnames(out_df) != opt\$protein_id_col])] # move ID column to first position


write.table(
    out_df,
    file = paste(output_prefix, 'proteus', 'raw_proteingroups_tab', 'tsv', sep = '.'),
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
limma.version <- as.character(packageVersion('limma'))
plotly.version <- as.character(packageVersion('plotly'))
proteus.version <- as.character(packageVersion('proteus'))
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-limma:', limma.version),
        paste('    r-plotly:', plotly.version),
        paste('    r-proteus-bartongroup:', proteus.version)
    ),
'versions.yml')
################################################
################################################
################################################
################################################