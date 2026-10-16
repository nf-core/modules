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

nullify <- function(x) {
    if (is.character(x) && (tolower(x) == "null" || x == "")) NULL else x
}

#' Write a table with the feature identifiers as the first column
#'
#' @param x Matrix or data frame with feature identifiers as row names
#' @param id_col Name of the identifier column
#' @param file Output file

write_feature_table <- function(x, id_col, file){
    out_df <- cbind(
        setNames(data.frame(rownames(x)), id_col),
        data.frame(x, check.names = FALSE)
    )
    write.table(
        out_df,
        file = file,
        col.names = TRUE,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

opt <- list(
    output_prefix              = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    count_file                 = '$counts',
    sample_file                = '$samplesheet',
    contrast_variable          = '$contrast_variable',
    reference_level            = '$reference',
    target_level               = '$target',
    formula                    = '$formula',
    contrast_string            = '$comparison',
    blocking_variables         = NULL,
    gene_id_col                = "gene_id",
    sample_id_col              = "experiment_accession",
    subset_to_contrast_samples = FALSE,
    exclude_samples_col        = NULL,
    exclude_samples_values     = NULL,
    filter_by_expr             = FALSE,     # filterByExpr
    min_count                  = 10,        # filterByExpr
    norm_method                = "TMM",     # normLibSizes
    test_method                = "qlf",     # 'qlf' (glmQLFit/glmQLFTest) or 'lrt' (glmFit/glmLRT)
    robust                     = TRUE,      # estimateDisp, glmQLFit
    prior_count                = 2,         # cpm
    adjust_method              = "BH",      # topTags
    round_digits               = NULL
)
opt_types <- lapply(opt, class)

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
        opt[[ao]] <- args_opt[[ao]]
    }
}
if ( ! is.null(opt\$round_digits)){
    opt\$round_digits <- as.numeric(opt\$round_digits)
}

keys <- c("formula", "contrast_string", "contrast_variable", "reference_level", "target_level")
opt[keys] <- lapply(opt[keys], nullify)

if (! opt\$test_method %in% c('qlf', 'lrt')){
    stop(paste0("Invalid test_method '", opt\$test_method, "': must be 'qlf' or 'lrt'"))
}

if (is.null(opt\$contrast_string) && any(sapply(opt[c('contrast_variable', 'reference_level', 'target_level')], is.null))){
    stop("Please provide either a comparison string, or a contrast variable with reference and target levels")
}

for (file_input in c('count_file', 'sample_file')){
    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

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
        row.names = opt\$gene_id_col
    )
sample.sheet <- read_delim_flexible(file = opt\$sample_file)

# Deal with spaces that may be in sample column
opt\$sample_id_col <- make.names(opt\$sample_id_col)

if (! opt\$sample_id_col %in% colnames(sample.sheet)){
    stop(paste0("Specified sample ID column '", opt\$sample_id_col, "' is not in the sample sheet"))
}

# Sample sheet can have duplicate rows for multiple sequencing runs, so uniqify
# before assigning row names

sample.sheet <- sample.sheet[! duplicated(sample.sheet[[opt\$sample_id_col]]), ]
rownames(sample.sheet) <- sample.sheet[[opt\$sample_id_col]]

missing_samples <-
    sample.sheet[!rownames(sample.sheet) %in% colnames(count.table), opt\$sample_id_col]

if (length(missing_samples) > 0) {
    stop(paste(
        length(missing_samples),
        'specified samples missing from count table:',
        paste(missing_samples, collapse = ',')
    ))
}

# Drop non-count columns (e.g. gene names) and order samples as in the sample sheet
count.table <- count.table[, rownames(sample.sheet)]

################################################
################################################
## CHECK CONTRAST SPECIFICATION               ##
################################################
################################################

blocking.vars <- c()

if (!is.null(opt\$contrast_variable)) {
    contrast_variable <- make.names(opt\$contrast_variable)

    if (!contrast_variable %in% colnames(sample.sheet)) {
        stop(paste0('Chosen contrast variable "', contrast_variable, '" not in sample sheet'))
    } else if (any(!c(opt\$reference_level, opt\$target_level) %in% sample.sheet[[contrast_variable]])) {
        stop(
            paste(
                'Please choose reference and target levels that are present in the',
                contrast_variable,
                'column of the sample sheet'
            )
        )
    } else if (!is.null(opt\$blocking_variables)) {
        blocking.vars <- make.names(unlist(strsplit(opt\$blocking_variables, split = ';')))
        if (!all(blocking.vars %in% colnames(sample.sheet))) {
            missing_block <- paste(blocking.vars[! blocking.vars %in% colnames(sample.sheet)], collapse = ',')
            stop(paste('Blocking variables', missing_block, 'do not correspond to sample sheet columns.'))
        }
    }

    if (opt\$subset_to_contrast_samples){
        selected_samples <- rownames(sample.sheet)[sample.sheet[[contrast_variable]] %in% c(opt\$target_level, opt\$reference_level)]
        count.table <- count.table[, selected_samples]
        sample.sheet <- sample.sheet[selected_samples, ]
    }
}

# Optionally, remove samples with specified values in a given field

if ((! is.null(opt\$exclude_samples_col)) && (! is.null(opt\$exclude_samples_values))){
    exclude_values <- unlist(strsplit(opt\$exclude_samples_values, split = ';'))

    if (! opt\$exclude_samples_col %in% colnames(sample.sheet)){
        stop(paste(opt\$exclude_samples_col, ' specified to subset samples is not a valid sample sheet column'))
    }

    print(paste0('Excluding samples with values of ', opt\$exclude_samples_values, ' in ', opt\$exclude_samples_col))
    selected_samples <- rownames(sample.sheet)[! sample.sheet[[opt\$exclude_samples_col]] %in% exclude_values]
    count.table <- count.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
}

################################################
################################################
## Build the model design                     ##
################################################
################################################

if (!is.null(opt\$formula)) {
    model <- opt\$formula
    cat("Using user-specified formula:\n   ", model, "\n")
    design <- model.matrix(as.formula(model), data = sample.sheet)
} else {
    # A zero-intercept design represents each contrast level directly
    model_vars <- c(contrast_variable, blocking.vars)
    model <- paste('~ 0 +', paste(model_vars, collapse = '+'))

    for (v in model_vars) {
        sample.sheet[[v]] <- as.factor(sample.sheet[[v]])
    }

    design <- model.matrix(as.formula(model), data = sample.sheet)
    colnames(design) <- sub(
        paste0('^', contrast_variable),
        paste0(contrast_variable, '.'),
        colnames(design)
    )
}
colnames(design) <- make.names(colnames(design))
cat("Design matrix columns:\n   ", paste(colnames(design), collapse = ", "), "\n")

if (!is.null(opt\$contrast_string)) {
    contrast_string <- as.character(opt\$contrast_string)
} else {
    treatment_target <- make.names(paste0(contrast_variable, ".", opt\$target_level))
    treatment_reference <- make.names(paste0(contrast_variable, ".", opt\$reference_level))

    # A level absent from the design (e.g. dropped as collinear) is compared against zero
    if ((treatment_target %in% colnames(design)) && (treatment_reference %in% colnames(design))) {
        contrast_string <- paste0(treatment_target, "-", treatment_reference)
    } else if (treatment_target %in% colnames(design)) {
        contrast_string <- paste0(treatment_target, "- 0")
    } else if (treatment_reference %in% colnames(design)) {
        contrast_string <- paste0("0 - ", treatment_reference)
    } else {
        stop(paste0(treatment_target, " and ", treatment_reference, " not found in design matrix"))
    }
}
cat("Using contrast string:", contrast_string, "\n")
contrast.matrix <- makeContrasts(contrasts = contrast_string, levels = colnames(design))

################################################
################################################
## Run edgeR                                  ##
################################################
################################################

dge <- DGEList(counts = as.matrix(count.table), samples = sample.sheet)

if (opt\$filter_by_expr) {
    keep <- filterByExpr(dge, design = design, min.count = opt\$min_count)
    cat("Keeping", sum(keep), "of", length(keep), "features after filterByExpr\n")
    dge <- dge[keep, , keep.lib.sizes = FALSE]
}

dge <- normLibSizes(dge, method = opt\$norm_method)
dge <- estimateDisp(dge, design = design, robust = opt\$robust)

if (opt\$test_method == 'qlf') {
    fit <- glmQLFit(dge, design = design, robust = opt\$robust)
    test <- glmQLFTest(fit, contrast = contrast.matrix)
} else {
    fit <- glmFit(dge, design = design)
    test <- glmLRT(fit, contrast = contrast.matrix)
}

comp.results <- topTags(
    test,
    n = Inf,
    adjust.method = opt\$adjust_method,
    sort.by = 'none'
)\$table

normalised_counts <- cpm(dge, log = TRUE, prior.count = opt\$prior_count)

################################################
################################################
## Generate outputs                           ##
################################################
################################################

if (! is.null(opt\$round_digits)){
    comp.results <- round(comp.results, opt\$round_digits)
    normalised_counts <- round(normalised_counts, opt\$round_digits)
}

write_feature_table(comp.results, opt\$gene_id_col, paste(opt\$output_prefix, 'edger.results.tsv', sep = '.'))
write_feature_table(normalised_counts, opt\$gene_id_col, paste(opt\$output_prefix, 'normalised_counts.tsv', sep = '.'))

png(
    file = paste(opt\$output_prefix, 'edger.mean_difference.png', sep = '.'),
    width = 600,
    height = 600
)
plotMD(test)
invisible(dev.off())

png(
    file = paste(opt\$output_prefix, 'edger.bcv.png', sep = '.'),
    width = 600,
    height = 600
)
plotBCV(dge)
invisible(dev.off())

saveRDS(fit, file = paste(opt\$output_prefix, 'DGEGLM.edger.rds', sep = '.'))

write(model, file = paste(opt\$output_prefix, 'edger.model.txt', sep = '.'))

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

writeLines(
    c(
        '"${task.process}":',
        paste('    bioconductor-edger:', as.character(packageVersion('edgeR')))
    ),
    'versions.yml'
)
