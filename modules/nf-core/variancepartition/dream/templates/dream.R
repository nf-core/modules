#!/usr/bin/env Rscript

# Load libraries
library("variancePartition")
library("edgeR")
library("BiocParallel")
library("limma")

# Load auxiliary helping functions

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse
parse_args <- function(x) {
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text = x, what = 'character', quiet = TRUE))
    args_vals <- lapply(args_vals, function(z) { length(z) <- 2; z })
    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[ !is.na(parsed_args) ]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#'
#' @return output Data frame
#'
read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = TRUE, stringsAsFactors = FALSE) {
    ext <- tools::file_ext(basename(file))
    if (ext == "tsv" || ext == "txt") {
        separator <- "\\t"
    } else if (ext == "csv") {
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }
    cat("Reading file", basename(file), "with", ext, "separator\n")
    df <- read.delim( file, sep = separator, header = header, row.names = row.names,
                    check.names = check.names, stringsAsFactors = stringsAsFactors)
    return(df)
}


# Options list
opt <- list(
    output_prefix      = "$meta.contrast_id",       # Prefix for output files
    count_file         = "$counts",                 # File containing raw counts
    sample_file        = "$samplesheet",            # File containing sample information
    contrast_variable  = "$meta.contrast_variable", # Variable for contrast (e.g., "treatment")
    contrast_reference    = "$meta.contrast_reference",# Reference level for the contrast
    contrast_target       = "$meta.contrast_target",   # Target level for the contrast (e.g., "mCherry")
    blocking_factors   = "$meta.blocking_factors",  # Blocking variables (e.g., "sample_number")
    sample_id_col      = "sample",    # Column name for sample IDs
    threads            = "$task.cpus",              # Number of threads for multithreading
    subset_to_contrast_samples = FALSE,            # Whether to subset to contrast samples
    exclude_samples_col = NULL,                    # Column for excluding samples
    exclude_samples_values = NULL,                 # Values for excluding samples
    number             = Inf,                     # Maximum number of results
    adjust.method      = "BH",                    # Adjustment method for topTable
    p.value            = 1,                       # P-value threshold for topTable
    lfc                = 0,                       # Log fold-change threshold for topTable
    confint            = FALSE,                   # Whether to compute confidence intervals in topTable
    ndups              = NULL,                    # Number of duplicates for lmFit
    spacing            = NULL,                    # Spacing for lmFit
    block              = NULL,                    # Block design for lmFit
    correlation        = NULL,                    # Correlation for lmFit
    method             = "ls",                    # Method for lmFit
    proportion         = 0.01,                    # Proportion for eBayes
    stdev_coef_lim     = "0.1,4",                 # Standard deviation coefficient limits for eBayes
    trend              = FALSE,                   # Whether to use trend in eBayes
    robust             = FALSE,                   # Whether to use robust method in eBayes
    winsor_tail_p      = "0.05,0.1",              # Winsor tail probabilities for eBayes
    ddf                = "adaptive",              # 'Satterthwaite', 'Kenward-Roger', or 'adaptive'
    reml               = FALSE
)

# Load external arguments to opt list
args_opt <- parse_args("$task.ext.args")
for (ao in names(args_opt)) {
    if (!ao %in% names(opt)) {
        stop(paste("Invalid option:", ao))
    }
    opt[[ao]] <- args_opt[[ao]]
}

# If there are no blocking factors, convert string "null" to NULL
if (!is.null(opt\$blocking_factors) && tolower(opt\$blocking_factors) == "null") {
    opt\$blocking_factors <- NULL
}

# Load metadata
metadata <- read_delim_flexible(opt\$sample_file, header = TRUE, stringsAsFactors = TRUE)
rownames(metadata) <- metadata[[opt\$sample_id_col]]

# Ensure contrast variable is factor, then relevel
if (!is.null(opt\$contrast_reference) && opt\$contrast_reference != "") {
    metadata[[opt\$contrast_variable]] <- factor(metadata[[opt\$contrast_variable]])
    metadata[[opt\$contrast_variable]] <- relevel(metadata[[opt\$contrast_variable]], ref = opt\$contrast_reference)
}

# Exclude samples in metadata if specified
if (!is.null(opt\$exclude_samples_col) && !is.null(opt\$exclude_samples_values)) {
  metadata <- metadata[!(metadata[[opt\$exclude_samples_col]] %in% opt\$exclude_samples_values), ]
}

# Load count data
countMatrix <- read_delim_flexible(opt\$count_file, header = TRUE, stringsAsFactors = FALSE)
rownames(countMatrix) <- countMatrix\$gene_id # count_file/matrix must have a gene_id column.
countMatrix <- countMatrix[, rownames(metadata), drop = FALSE]

# Standard usage of limma/voom
dge <- DGEList(countMatrix)
dge <- calcNormFactors(dge)

# Print sample names
cat("Response sample names (dge counts):\n")
print(colnames(dge\$counts))
cat("Metadata sample names (rownames(metadata)):\n")
print(rownames(metadata))

# Construct model formula using contrast and, if available, blocking variables
if (!is.null(opt\$blocking_factors) && opt\$blocking_factors != "") {
  form <- as.formula(paste0("~ ", opt\$contrast_variable, " + (1 | ", opt\$blocking_factors, ")"))
} else {
  form <- as.formula(paste0("~ ", opt\$contrast_variable))
}

# Parallel processing setup
threads <- as.numeric(opt\$threads)
param <- SnowParam(threads, "SOCK", progressbar = TRUE)

# Estimate weights using the linear mixed model of DREAM
vobjDream <- voomWithDreamWeights(dge, form, metadata, BPPARAM = param)

# Fit the DREAM model with ddf and reml options
fitmm <- dream(vobjDream, form, metadata, ddf = opt\$ddf, reml = opt\$reml)

# Parse stdev_coef_lim and winsor_tail_p into numeric vectors
stdev_coef_lim_vals <- as.numeric(unlist(strsplit(opt\$stdev_coef_lim, ",")))
winsor_tail_p_vals  <- as.numeric(unlist(strsplit(opt\$winsor_tail_p, ",")))

# Fit the empirical Bayes model with additional parameters
fitmm <- eBayes(fitmm, proportion = opt\$proportion,
                stdev.coef.lim = stdev_coef_lim_vals,
                trend = opt\$trend, robust = opt\$robust,
                winsor.tail.p = winsor_tail_p_vals)

# Display design matrix
head(fitmm\$design, 3)
print(colnames(fitmm\$design))

# Write coefficient name from contrast variable and target level (e.g., "treatmentmCherry")
coef_name <- paste0(opt\$contrast_variable, opt\$contrast_target)

# Get topTable results
results <- topTable(fitmm, coef = coef_name,
                    adjust.method = opt\$adjust.method, p.value = opt\$p.value,
                    lfc = opt\$lfc, confint = opt\$confint)

# Export topTable results
write.table(results, file = paste(opt\$output_prefix, 'dream.results.tsv', sep = '.'),
            col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE )

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
edger.version <- as.character(packageVersion('edgeR'))
variancePartition.version <- as.character(packageVersion('variancePartition'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-edger:', edger.version),
        paste('    variancePartition:', variancePartition.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
