#!/usr/bin/env Rscript

# Load libraries
library("variancePartition")
library("edgeR")
library("BiocParallel")
library("limma")

# Load auxiliary helping functions

#' Check for Non-Empty, Non-Whitespace String
#'
#' This function checks if the input is non-NULL and contains more than just whitespace.
#' It returns TRUE if the input is a non-empty, non-whitespace string, and FALSE otherwise.
#'
#' @param input A variable to check.
#' @return A logical value: TRUE if the input is a valid, non-empty, non-whitespace string; FALSE otherwise.
#' @examples
#' is_valid_string("Hello World") # Returns TRUE
#' is_valid_string("   ")         # Returns FALSE
#' is_valid_string(NULL)          # Returns FALSE
is_valid_string <- function(input) {
  !is.null(input) && nzchar(trimws(input))
}

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

#
#' Turn “null” or empty strings into actual NULL
#'
#' @param x Input option
#'
#' @return NULL or x
#'
nullify <- function(x) {
  if (is.character(x) && (tolower(x) == "null" || x == "")) NULL else x
}

# Options list
opt <- list(
    output_prefix              = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    count_file                 = "$counts",               # File containing raw counts
    sample_file                = "$samplesheet",          # File containing sample information
    blocking_variables         = NULL,
    contrast_variable          = "$contrast_variable",    # Variable for contrast (e.g., "treatment")
    contrast_reference         = "$reference",            # Reference level for the contrast
    contrast_target            = "$target",               # Target level for the contrast (e.g., "mCherry")
    contrast_string            = "$comparison",           # Full (complex) contrast expression comparison needed if using formula
    sample_id_col              = "sample",                # Column name for sample IDs
    threads                    = "$task.cpus",            # Number of threads for multithreading
    subset_to_contrast_samples = FALSE,                   # Whether to subset to contrast samples
    exclude_samples_col        = NULL,                    # Column for excluding samples
    exclude_samples_values     = NULL,                    # Values for excluding samples
    adjust.method              = "BH",                    # Adjustment method for topTable
    p.value                    = 1,                       # P-value threshold for topTable
    lfc                        = 0,                       # Log fold-change threshold for topTable
    confint                    = FALSE,                   # Whether to compute confidence intervals in topTable
    proportion                 = 0.01,                    # Proportion for eBayes
    stdev_coef_lim             = "0.1,4",                 # Standard deviation coefficient limits for eBayes
    trend                      = FALSE,                   # Whether to use trend in eBayes
    robust                     = FALSE,                   # Whether to use robust method in eBayes
    winsor_tail_p              = "0.05,0.1",              # Winsor tail probabilities for eBayes
    ddf                        = "adaptive",              # 'Satterthwaite', 'Kenward-Roger', or 'adaptive'
    reml                       = FALSE,
    round_digits               = NULL,
    formula                    = "$formula",              # User-specified formula (e.g. "~ + (1 | sample_number)")
    apply_voom                 = FALSE                    # Whether to apply `voomWithDreamWeights`
)

# Load external arguments to opt list
args_opt <- parse_args("$task.ext.args")
for (ao in names(args_opt)) {
    if (!ao %in% names(opt)) {
        stop(paste("Invalid option:", ao))
    }
    opt[[ao]] <- args_opt[[ao]]
}

# If there is no option supplied, convert string "null" to NULL
keys <- c("formula", "contrast_string", "contrast_target", "contrast_variable", "blocking_variables", "contrast_reference")
opt[keys] <- lapply(opt[keys], nullify)

opt\$threads      <- as.numeric(opt\$threads)
opt\$apply_voom   <- as.logical(opt\$apply_voom)
opt\$proportion   <- as.numeric(opt\$proportion)
opt\$trend        <- as.logical(opt\$trend)
opt\$robust       <- as.logical(opt\$robust)
opt\$reml         <- as.logical(opt\$reml)
opt\$p.value      <- as.numeric(opt\$p.value)
opt\$lfc          <- as.numeric(opt\$lfc)
opt\$confint      <- as.logical(opt\$confint)

if (!is.null(opt\$round_digits)){
  opt\$round_digits <- as.numeric(opt\$round_digits)
}

# Load metadata
metadata <- read_delim_flexible(opt\$sample_file, header = TRUE, stringsAsFactors = TRUE)
rownames(metadata) <- metadata[[opt\$sample_id_col]]

# Check if required parameters have been provided
if (is_valid_string(opt\$formula)) {
  contrast_tuple <- c('contrast_variable', 'contrast_reference', 'contrast_target', 'blocking_variables')
  offending <- vapply(
    opt[contrast_tuple],
    is_valid_string,
    logical(1)
  )
  offending_opts <- contrast_tuple[offending]
  if (length(offending_opts) > 0) {
    stop(paste("When 'formula' is provided, contrasts must be specified only via 'contrast_string'.\n",
        "The following options should not be set:",
        paste(offending_opts, collapse = ', ')))
  }
  required_opts <- c('output_prefix', 'contrast_string')
} else {
  required_opts <- c('contrast_variable', 'contrast_reference', 'contrast_target', 'output_prefix')
}
missing <- required_opts[!unlist(lapply(opt[required_opts], is_valid_string)) | !required_opts %in% names(opt)]

if (length(missing) > 0) {
  stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

if (!is_valid_string(opt\$formula)) {
  contrast_variable <- make.names(opt\$contrast_variable)
  blocking.vars <- c()

  if (!contrast_variable %in% colnames(metadata)) {
    stop(
      paste0(
        'Chosen contrast variable \"',
        contrast_variable,
        '\" not in sample sheet'
      )
    )
  } else if (any(!c(opt\$contrast_reference, opt\$contrast_target) %in% metadata[[contrast_variable]])) {
    stop(
      paste(
        'Please choose reference and treatment levels that are present in the',
        contrast_variable,
        'column of the sample sheet'
      )
    )
  } else if (is_valid_string(opt\$blocking_variables)) {
    blocking.vars = make.names(unlist(strsplit(opt\$blocking_variables, split = ';')))
    if (!all(blocking.vars %in% colnames(metadata))) {
      missing_block <- paste(blocking.vars[! blocking.vars %in% colnames(metadata)], collapse = ',')
      stop(
        paste(
          'Blocking variables', missing_block,
          'do not correspond to sample sheet columns.'
        )
      )
    }
  }
}

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
countMatrix <- read_delim_flexible(opt\$count_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
rownames(countMatrix) <- countMatrix\$gene_id # count_file/matrix must have a gene_id column.
countMatrix <- countMatrix[, rownames(metadata), drop = FALSE]

if (opt\$subset_to_contrast_samples) {
  sample_selector <- metadata[[contrast_variable]] %in%
    c(opt\$contrast_target, opt\$contrast_reference)

  selected_samples <- rownames(metadata)[sample_selector]

  metadata <- metadata[selected_samples, , drop = FALSE]
  countMatrix <- countMatrix[, selected_samples, drop = FALSE]
}

# Construct model formula using user-provided formula if available; if not, generate formula from variables
if (is_valid_string(opt\$formula)) {
    form <- as.formula(opt\$formula)
} else {
    form <- '~ 0'

    if (is_valid_string(opt\$blocking_variables)) {
        form <- paste(form, paste(blocking.vars, collapse = ' + '), sep=' + ')
    }

    # Make sure all the appropriate variables are factors
    for (v in c(blocking.vars, contrast_variable)) {
        metadata[[v]] <- as.factor(make.names(metadata[[v]]))
    }

    # Variable of interest goes last
    form <- as.formula(paste(form, contrast_variable, sep = ' + '))
}
print(form)

# Parallel processing setup
threads <- as.numeric(opt\$threads)

# Optionally apply voom
if (as.logical(opt\$apply_voom)) {
    # Standard usage of limma/voom
    dge <- DGEList(countMatrix)
    dge <- calcNormFactors(dge)
    vobjDream <- voomWithDreamWeights(dge, form, metadata, BPPARAM = MulticoreParam(threads))

    # Write normalized counts matrix to a TSV file
     normalized_counts <- vobjDream\$E
    if (!is.null(opt\$round_digits)) {
        normalized_counts <- apply(normalized_counts, 2, function(x) round(x, opt\$round_digits))
    }
    normalized_counts_with_genes <- data.frame(gene_id = rownames(normalized_counts), normalized_counts, check.names = FALSE, row.names = NULL)
    write.table(normalized_counts_with_genes,
        file = paste(opt\$output_prefix, "normalised_counts.tsv", sep = '.'),
        sep = "	",
        quote = FALSE,
        row.names = FALSE)
} else {
    # Assume countMatrix roughly follows a normal distribution
    vobjDream <- countMatrix
}

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

# If contrast_string is provided, use that for makeContrast
if (!is.null(opt\$contrast_string)) {
    contrast_string <- opt\$contrast_string
} else {
    # Construct the contrast_string
    treatment_target <- make.names(paste0(opt\$contrast_variable, opt\$contrast_target))
    treatment_reference <- make.names(paste0(opt\$contrast_variable, opt\$contrast_reference))
    contrast_string <- paste0(treatment_target, "-", treatment_reference)
}

# Use makeContrasts if contrast_string exists
if (is_valid_string(contrast_string)) {
    cat("Using contrast string:", contrast_string, "\n")

    colnames(fitmm\$design) <- make.names(colnames(fitmm\$design))
    contrast_matrix <- makeContrasts(contrast = contrast_string, levels = colnames(fitmm\$design))
    fit2 <- contrasts.fit(fitmm, contrast_matrix)
    fit2 <- eBayes(fit2, proportion = opt\$proportion,
                  stdev.coef.lim = stdev_coef_lim_vals,
                  trend = opt\$trend, robust = opt\$robust,
                  winsor.tail.p = winsor_tail_p_vals)
    results <- topTable(fit2, number = Inf,
                        adjust.method = opt\$adjust.method,
                        p.value = opt\$p.value, lfc = opt\$lfc, confint = opt\$confint)
}

results\$gene_id <- rownames(results)
results <- results[, c("gene_id", setdiff(names(results), "gene_id"))]

# Round results if required
if (!is.null(opt\$round_digits)) {
    numeric_columns <- vapply(results, is.numeric, logical(1))
    results[numeric_columns] <- lapply(results[numeric_columns], round, digits = opt\$round_digits)
}

# Export topTable results
write.table(results, file = paste(opt\$output_prefix, 'dream.results.tsv', sep = '.'),
            col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE )

# Save model to file
write(deparse(form), file=paste(opt\$output_prefix, 'dream.model.txt', sep = '.'))

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
