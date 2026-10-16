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
    parsed_args[ ( ! parsed_args %in%  c('', 'null')) & ! is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = FALSE){

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

#' Turn "null" or empty strings into actual NULL
#'
#' @param x Input option
#'
#' @return NULL or x
#'
nullify <- function(x) {
  if (is.character(x) && (tolower(x) == "null" || x == "")) NULL else x
}

#' Stop if an option is not one of the allowed values
#'
#' @param name Option name in opt
#' @param choices Allowed values

check_choice <- function(name, choices) {
    if (! opt[[name]] %in% choices) {
        stop(paste0("Invalid ", name, " '", opt[[name]], "': must be one of ", paste(choices, collapse = ', ')))
    }
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

# Options settable through ext.args (--option value):
#   probe_id_col                Feature ID column in the count matrix (gene_id)
#   sample_id_col               Sample ID column in the sample sheet (experiment_accession)
#   blocking_variables          ';'-separated sample sheet columns added to the
#                               contrast-variable design; ignored with a formula
#   subset_to_contrast_samples  Keep only reference and target samples (FALSE)
#   exclude_samples_col         Sample sheet column used to exclude samples
#   exclude_samples_values      ';'-separated values of exclude_samples_col to drop
#   filter_by_expr              Filter features with filterByExpr() (FALSE)
#   norm_method                 normLibSizes() method: TMM, TMMwsp, RLE,
#                               upperquartile or none (TMM)
#   robust                      Robust estimateDisp() and glmQLFit() (TRUE)
#   test                        QLF (glmQLFit/glmQLFTest) or LRT (glmFit/glmLRT)
#   adjust.method               p.adjust() method used by topTags() (BH)
#   prior_count                 Prior count for log2-CPM normalised counts (2)
#   round_digits                Decimal places for output tables (no rounding)
opt <- list(
    output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    count_file = '$counts',
    sample_file = '$samplesheet',
    contrast_variable = '$contrast_variable',
    reference_level = '$reference',
    target_level = '$target',
    formula = '$formula',
    contrast_string = '$comparison',
    blocking_variables = NULL,
    probe_id_col = "gene_id",
    sample_id_col = "experiment_accession",
    subset_to_contrast_samples = FALSE,
    exclude_samples_col = NULL,
    exclude_samples_values = NULL,
    filter_by_expr = FALSE,
    norm_method = "TMM",
    robust = TRUE,
    test = "QLF",
    adjust.method = "BH",
    prior_count = 2,
    round_digits = NULL
)
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
if ( ! is.null(opt\$round_digits)){
    opt\$round_digits <- as.numeric(opt\$round_digits)
}

# If there is no option supplied, convert string "null" to NULL
keys <- c("formula", "contrast_string", "contrast_variable", "reference_level", "target_level")
opt[keys] <- lapply(opt[keys], nullify)

check_choice('test', c('QLF', 'LRT'))
check_choice('norm_method', c('TMM', 'TMMwsp', 'RLE', 'upperquartile', 'none'))
check_choice('adjust.method', p.adjust.methods)

# Check file inputs are valid

for (file_input in c('count_file', 'sample_file')){
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

library(limma)
library(edgeR)

################################################
################################################
## READ IN COUNTS FILE AND SAMPLE METADATA    ##
################################################
################################################

count.table <-
    read_delim_flexible(
        file = opt\$count_file,
        header = TRUE,
        row.names = opt\$probe_id_col
    )
# Syntactic column names let sample sheet columns be used in model formulas;
# user-supplied column names are passed through make.names() to match
sample.sheet <- read_delim_flexible(file = opt\$sample_file, check.names = TRUE)

opt\$sample_id_col <- make.names(opt\$sample_id_col)

if (! opt\$sample_id_col %in% colnames(sample.sheet)){
    stop(paste0("Specified sample ID column '", opt\$sample_id_col, "' is not in the sample sheet"))
}

# Sample sheet can have duplicate rows for multiple sequencing runs, so uniqify
# before assigning row names

sample.sheet <- sample.sheet[! duplicated(sample.sheet[[opt\$sample_id_col]]), ]
rownames(sample.sheet) <- sample.sheet[[opt\$sample_id_col]]

# Check that all samples specified in the input sheet are present in the count
# table. Assuming they are, subset and sort the count table to match the
# sample sheet, dropping any non-count annotation columns (e.g. gene_name)

missing_samples <-
    sample.sheet[!rownames(sample.sheet) %in% colnames(count.table), opt\$sample_id_col]

if (length(missing_samples) > 0) {
    stop(paste(
        length(missing_samples),
        'specified samples missing from count table:',
        paste(missing_samples, collapse = ',')
    ))
} else{
    count.table <- count.table[, rownames(sample.sheet), drop = FALSE]
}

################################################
################################################
## CHECK CONTRAST SPECIFICATION               ##
################################################
################################################

contrast_variable <- NULL
blocking.vars <- c()

if (!is.null(opt\$contrast_variable)) {
    contrast_variable <- make.names(opt\$contrast_variable)
    if (!contrast_variable %in% colnames(sample.sheet)) {
        stop(
            paste0(
            'Chosen contrast variable \"',
            contrast_variable,
            '\" not in sample sheet'
            )
        )
    } else if (any(!c(opt\$reference_level, opt\$target_level) %in% sample.sheet[[contrast_variable]])) {
        stop(
            paste(
            'Please choose reference and target levels that are present in the',
            contrast_variable,
            'column of the sample sheet'
            )
        )
    } else if (!is.null(opt\$blocking_variables)) {
        blocking.vars = make.names(unlist(strsplit(opt\$blocking_variables, split = ';')))
        if (!all(blocking.vars %in% colnames(sample.sheet))) {
            missing_block <- paste(blocking.vars[! blocking.vars %in% colnames(sample.sheet)], collapse = ',')
            stop(
                paste(
                    'Blocking variables', missing_block,
                    'do not correspond to sample sheet columns.'
                )
            )
        }
    }
} else if (!is.null(opt\$blocking_variables)) {
    warning("blocking_variables is ignored when a formula is supplied; include blocking terms in the formula instead")
}

# Optionally, subset to only the samples involved in the contrast

if (opt\$subset_to_contrast_samples){
    if (is.null(contrast_variable)){
        stop("subset_to_contrast_samples requires a contrast_variable")
    }
    sample_selector <- sample.sheet[[contrast_variable]] %in% c(opt\$target_level, opt\$reference_level)
    selected_samples <- sample.sheet[sample_selector, opt\$sample_id_col]
    count.table <- count.table[, selected_samples, drop = FALSE]
    sample.sheet <- sample.sheet[selected_samples, ]
}

# Optionally, remove samples with specified values in a given field (probably
# don't use this as well as the above)

if ((! is.null(opt\$exclude_samples_col)) && (! is.null(opt\$exclude_samples_values))){
    exclude_values = unlist(strsplit(opt\$exclude_samples_values, split = ';'))
    opt\$exclude_samples_col <- make.names(opt\$exclude_samples_col)

    if (! opt\$exclude_samples_col %in% colnames(sample.sheet)){
        stop(paste(opt\$exclude_samples_col, ' specified to subset samples is not a valid sample sheet column'))
    }

    print(paste0('Excluding samples with values of ', opt\$exclude_samples_values, ' in ', opt\$exclude_samples_col))
    sample_selector <- ! sample.sheet[[opt\$exclude_samples_col]] %in% exclude_values

    selected_samples <- sample.sheet[sample_selector, opt\$sample_id_col]
    count.table <- count.table[, selected_samples, drop = FALSE]
    sample.sheet <- sample.sheet[selected_samples, ]
}

################################################
################################################
## Build the model                            ##
################################################
################################################

if (!is.null(opt\$formula)) {
    model <- opt\$formula
    model_formula <- as.formula(model)
    cat("Using user-specified formula:\n   ", deparse(model_formula), "\n")
    design <- model.matrix(model_formula, data = sample.sheet)
    colnames(design) <- make.names(colnames(design))
    cat("Column names after make.names():\n   ", paste(colnames(design), collapse = ", "), "\n")

} else {
    if (is.null(contrast_variable)){
        stop("Please provide either a formula and comparison, or a contrast variable with reference and target levels")
    }

    # Put the contrast variable first in zero-intercept designs so each
    # contrast level is represented directly in the design matrix.
    model_vars <- c(contrast_variable, blocking.vars)
    model <- paste('~ 0 +', paste(model_vars, collapse = '+'))

    for (v in model_vars) {
        sample.sheet[[v]] <- as.factor(sample.sheet[[v]])
    }

    design <- model.matrix(as.formula(model), data = sample.sheet)

    # The contrast variable comes first, so its levels are the leading columns
    contrast_levels <- levels(sample.sheet[[contrast_variable]])
    colnames(design)[seq_along(contrast_levels)] <- paste0(contrast_variable, '.', contrast_levels)
    colnames(design) <- make.names(colnames(design))
    cat("Final column names after make.names():\n")
    print(colnames(design))
}

################################################
################################################
## Contrast                                   ##
################################################
################################################

if (!is.null(opt\$contrast_string)) {
    cat("Using contrast string:", opt\$contrast_string, "\n")
    contrast_string <- as.character(opt\$contrast_string)
} else {
    treatment_target <- make.names(paste0(contrast_variable, ".", opt\$target_level))
    treatment_reference <- make.names(paste0(contrast_variable, ".", opt\$reference_level))

    if ((treatment_target %in% colnames(design)) && (treatment_reference %in% colnames(design))) {
        contrast_string <- paste0(treatment_target, "-", treatment_reference)
    } else if (treatment_target %in% colnames(design)) {
        contrast_string <- paste0(treatment_target, "- 0")
    } else if (treatment_reference %in% colnames(design)) {
        contrast_string <- paste0("0 - ", treatment_reference)
    } else {
        stop(paste0(treatment_target, " and ", treatment_reference, " not found in design matrix"))
    }
    cat("Using contrast string:", contrast_string, "\n")
}
contrast.matrix <- makeContrasts(contrasts = contrast_string, levels = colnames(design))

################################################
################################################
## Run edgeR                                  ##
################################################
################################################

count.matrix <- as.matrix(count.table)

# A 'group' column in the sample sheet becomes the DGEList group; it is unused
# because every edgeR call below is given the design
y <- DGEList(counts = count.matrix, samples = sample.sheet)

if (as.logical(opt\$filter_by_expr)) {
    keep <- filterByExpr(y, design)
    cat("Keeping", sum(keep), "of", length(keep), "features after filterByExpr()\n")
    y <- y[keep, , keep.lib.sizes = FALSE]
}

y <- normLibSizes(y, method = opt\$norm_method)

y <- estimateDisp(y, design, robust = as.logical(opt\$robust))

if (opt\$test == 'QLF') {
    fit <- glmQLFit(y, design, robust = as.logical(opt\$robust))
    res <- glmQLFTest(fit, contrast = contrast.matrix)
} else {
    fit <- glmFit(y, design)
    res <- glmLRT(fit, contrast = contrast.matrix)
}

comp.results <- topTags(
    res,
    n = Inf,
    adjust.method = opt\$adjust.method,
    sort.by = 'none'
)\$table

################################################
################################################
## Generate outputs                           ##
################################################
################################################

cat("Saving results for ", contrast_string, " ...\n", sep = "")

# Differential expression table

if (! is.null(opt\$round_digits)){
    comp.results <- round(comp.results, opt\$round_digits)
}
out_df <- cbind(
    setNames(data.frame(rownames(comp.results)), opt\$probe_id_col),
    data.frame(comp.results, check.names = FALSE, row.names = NULL)
)
write.table(
    out_df,
    file = paste(opt\$output_prefix, 'edger.results.tsv', sep = '.'),
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

# Normalised counts (log2 CPM) for all features, including any removed by
# filterByExpr(), using the effective library sizes of the fitted data

norm_counts <- cpm(
    count.matrix,
    lib.size = getNormLibSizes(y),
    log = TRUE,
    prior.count = as.numeric(opt\$prior_count)
)
if (! is.null(opt\$round_digits)){
    norm_counts <- round(norm_counts, opt\$round_digits)
}
norm_counts_df <- data.frame(
    setNames(data.frame(rownames(norm_counts)), opt\$probe_id_col),
    norm_counts,
    check.names = FALSE,
    row.names = NULL
)
write.table(
    norm_counts_df,
    file = paste(opt\$output_prefix, 'normalised_counts.tsv', sep = '.'),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
)

# Library sizes and normalisation factors

write.table(
    data.frame(
        sample = rownames(y\$samples),
        lib.size = y\$samples\$lib.size,
        norm.factors = y\$samples\$norm.factors
    ),
    file = paste(opt\$output_prefix, 'edger.norm_factors.tsv', sep = '.'),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
)

# Plots

png(file = paste(opt\$output_prefix, 'edger.mean_difference.png', sep = '.'), width = 600, height = 600)
plotMD(res)
dev.off()

png(file = paste(opt\$output_prefix, 'edger.bcv.png', sep = '.'), width = 600, height = 600)
plotBCV(y)
dev.off()

# R object for other processes to use

saveRDS(fit, file = paste(opt\$output_prefix, 'DGEGLM.edger.rds', sep = '.'))

# Save model and contrast to file

write(
    c(model, paste('contrast:', contrast_string), paste('test:', opt\$test)),
    file = paste(opt\$output_prefix, 'edger.model.txt', sep = '.')
)

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(opt\$output_prefix, "R_sessionInfo.log", sep = '.'))
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
edger.version <- as.character(packageVersion('edgeR'))
limma.version <- as.character(packageVersion('limma'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-edger:', edger.version),
        paste('    bioconductor-limma:', limma.version)
    ),
'versions.yml')
