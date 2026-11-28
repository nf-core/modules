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
#' @param check.names Passed to read.delim()
#' @param character_cols Columns to read as character, e.g. identifiers with leading zeros
#' @param comment_char Passed to read.delim() as comment.char
#' @param na_strings Passed to read.delim() as na.strings
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, check.names = FALSE, character_cols = NULL, comment_char = "", na_strings = "NA"){

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
        check.names = check.names,
        comment.char = comment_char,
        na.strings = na_strings,
        colClasses = if (is.null(character_cols)) NA else setNames(rep("character", length(character_cols)), character_cols)
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

#' Round the numeric columns of a data frame, leaving annotation columns as they are
#'
#' @param x Data frame
#' @param digits Number of decimal places, or NULL for no rounding
#'
#' @return Data frame

round_numeric <- function(x, digits){
    if (is.null(digits)) return(x)
    numeric_cols <- vapply(x, is.double, logical(1))
    x[numeric_cols] <- lapply(x[numeric_cols], round, digits = digits)
    x
}

#' Write a results table without row names
#'
#' @param x Data frame
#' @param file Output file

write_results_table <- function(x, file){
    write.table(
        x,
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
    feature_id_col             = "feature_id",
    annotation_cols            = NULL,
    count_col_suffix           = NULL,
    sample_id_col              = "experiment_accession",
    subset_to_contrast_samples = FALSE,
    exclude_samples_col        = NULL,
    exclude_samples_values     = NULL,
    filter_by_expr             = FALSE,     # filterByExpr
    min_count                  = 10,        # filterByExpr
    norm_method                = "TMM",     # normLibSizes
    test_method                = "qlf",     # 'qlf' (glmQLFit) or 'lrt' (glmFit)
    robust                     = TRUE,      # estimateDisp, glmQLFit
    adjust_method              = "BH",      # diffSpliceDGE, topSpliceDGE
    plot_n                     = 10,        # plotSpliceDGE
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
            converted <- suppressWarnings(as(args_opt[[ao]], opt_types[[ao]]))
            if (is.na(converted)){
                stop(paste0("Invalid value for option ", ao, ": '", args_opt[[ao]], "' is not a valid ", opt_types[[ao]]))
            }
            args_opt[[ao]] <- converted
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

if (!is.null(opt\$formula) && is.null(opt\$contrast_string)){
    stop("A formula requires a comparison string")
}
if (!is.null(opt\$contrast_string) && is.null(opt\$formula) && is.null(opt\$contrast_variable)){
    stop("A comparison string requires either a formula or a contrast variable")
}
if (is.null(opt\$contrast_string) && any(sapply(opt[c('contrast_variable', 'reference_level', 'target_level')], is.null))){
    stop("Please provide either a formula with a comparison string, or a contrast variable with reference and target levels")
}
if (opt\$subset_to_contrast_samples && any(sapply(opt[c('contrast_variable', 'reference_level', 'target_level')], is.null))){
    stop("subset_to_contrast_samples requires a contrast variable with reference and target levels")
}
if (!is.null(opt\$formula) && !is.null(opt\$blocking_variables)){
    stop("blocking_variables cannot be used with a formula: add the blocking variables to the formula instead")
}
if (!is.null(opt\$reference_level) && identical(opt\$reference_level, opt\$target_level)){
    stop(paste0("Reference and target levels must differ (both are '", opt\$reference_level, "')"))
}
if (identical(opt\$gene_id_col, opt\$feature_id_col)){
    stop("gene_id_col and feature_id_col must be different columns")
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

annotation_cols <- c(opt\$gene_id_col, opt\$feature_id_col)
if (! is.null(opt\$annotation_cols)){
    annotation_cols <- unique(c(annotation_cols, unlist(strsplit(opt\$annotation_cols, split = ';'))))
}

# Skipping '#' lines lets featureCounts tables be read as they are
count.table <- read_delim_flexible(file = opt\$count_file, character_cols = annotation_cols, comment_char = '#')

if (! is.null(opt\$count_col_suffix)){
    count_cols <- colnames(count.table)
    has_suffix <- endsWith(count_cols, opt\$count_col_suffix)
    count_cols[has_suffix] <- substr(count_cols[has_suffix], 1, nchar(count_cols[has_suffix]) - nchar(opt\$count_col_suffix))
    colnames(count.table) <- count_cols
}

missing_annotation <- annotation_cols[! annotation_cols %in% colnames(count.table)]
if (length(missing_annotation) > 0){
    stop(paste('Annotation columns missing from count table:', paste(missing_annotation, collapse = ',')))
}

# Headers and column options both go through make.names() so they stay aligned. Sample IDs
# are read as character so that IDs such as '001' still match the count table columns, and
# empty cells are read as missing rather than as an empty factor level.
opt\$sample_id_col <- make.names(opt\$sample_id_col)
sample.sheet <- read_delim_flexible(
    file = opt\$sample_file,
    check.names = TRUE,
    character_cols = opt\$sample_id_col,
    na_strings = c("", "NA")
)

if (! opt\$sample_id_col %in% colnames(sample.sheet)){
    stop(paste0("Specified sample ID column '", opt\$sample_id_col, "' is not in the sample sheet"))
}

# Sample sheet can have duplicate rows for multiple sequencing runs, so uniqify
# before assigning row names

sample.sheet <- sample.sheet[! duplicated(sample.sheet[[opt\$sample_id_col]]), ]
rownames(sample.sheet) <- sample.sheet[[opt\$sample_id_col]]

feature.annotation <- count.table[, annotation_cols, drop = FALSE]
feature.ids <- paste(feature.annotation[[opt\$gene_id_col]], feature.annotation[[opt\$feature_id_col]], sep = ':')
if (any(duplicated(feature.ids))) {
    stop(paste0(
        "Combinations of '", opt\$gene_id_col, "' and '", opt\$feature_id_col, "' must be unique, e.g. ",
        paste(head(unique(feature.ids[duplicated(feature.ids)]), 5), collapse = ', '),
        " occur more than once. Count overlapping exons from a flattened annotation"
    ))
}
rownames(feature.annotation) <- feature.ids

################################################
################################################
## CHECK CONTRAST SPECIFICATION               ##
################################################
################################################

contrast_variable <- if (is.null(opt\$contrast_variable)) NULL else make.names(opt\$contrast_variable)
blocking.vars <- c()

if (!is.null(contrast_variable)) {
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
        sample.sheet <- sample.sheet[selected_samples, ]
    }
}

# Optionally, remove samples with specified values in a given field

if ((! is.null(opt\$exclude_samples_col)) && (! is.null(opt\$exclude_samples_values))){
    opt\$exclude_samples_col <- make.names(opt\$exclude_samples_col)
    exclude_values <- unlist(strsplit(opt\$exclude_samples_values, split = ';'))

    if (! opt\$exclude_samples_col %in% colnames(sample.sheet)){
        stop(paste(opt\$exclude_samples_col, ' specified to subset samples is not a valid sample sheet column'))
    }

    print(paste0('Excluding samples with values of ', opt\$exclude_samples_values, ' in ', opt\$exclude_samples_col))
    selected_samples <- rownames(sample.sheet)[! sample.sheet[[opt\$exclude_samples_col]] %in% exclude_values]
    sample.sheet <- sample.sheet[selected_samples, ]

    if (is.null(opt\$contrast_string) && any(!c(opt\$reference_level, opt\$target_level) %in% sample.sheet[[contrast_variable]])) {
        stop('Sample exclusion removed all samples of the reference or target level')
    }
}

# Samples are checked against the count table only once excluded samples are gone
missing_samples <- rownames(sample.sheet)[! rownames(sample.sheet) %in% colnames(count.table)]

if (length(missing_samples) > 0) {
    stop(paste(
        length(missing_samples),
        'specified samples missing from count table:',
        paste(missing_samples, collapse = ',')
    ))
}

# Drop non-count columns and order samples as in the sample sheet
count.table <- count.table[, rownames(sample.sheet), drop = FALSE]
rownames(count.table) <- feature.ids

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
    # Blocking variables are checked on the samples that enter the model, i.e. after any subsetting
    single_level <- vapply(blocking.vars, function(v) length(unique(na.omit(sample.sheet[[v]]))) < 2, logical(1))
    for (v in blocking.vars[single_level]) {
        warning("Blocking variable '", v, "' has one level in this contrast and is left out of the model", call. = FALSE)
    }
    blocking.vars <- blocking.vars[! single_level]

    # A zero-intercept design represents each contrast level directly
    model_vars <- c(contrast_variable, blocking.vars)
    model <- paste('~ 0 +', paste(model_vars, collapse = '+'))

    for (v in model_vars) {
        sample.sheet[[v]] <- as.factor(sample.sheet[[v]])
    }

    design <- model.matrix(as.formula(model), data = sample.sheet)

    contrast_levels <- levels(sample.sheet[[contrast_variable]])
    contrast_cols <- match(paste0(contrast_variable, contrast_levels), colnames(design))
    colnames(design)[contrast_cols] <- paste0(contrast_variable, '.', contrast_levels)
}

# model.matrix() silently drops samples with missing values in the model variables
if (nrow(design) != nrow(sample.sheet)) {
    model_cols <- intersect(all.vars(as.formula(model)), colnames(sample.sheet))
    na_samples <- rownames(sample.sheet)[! complete.cases(sample.sheet[, model_cols, drop = FALSE])]
    stop(paste0(
        "Samples with missing values in the model variables (", paste(model_cols, collapse = ', '), "): ",
        paste(na_samples, collapse = ', ')
    ))
}

# Checked before fitting so that rank deficiency gives a clear message
if (is.null(opt\$formula) && length(blocking.vars) > 0) {
    design_rank <- qr(design)\$rank
    if (design_rank < ncol(design)) {
        stop(paste0(
            "Blocking variables ", paste(blocking.vars, collapse = ', '), " are confounded with '", contrast_variable, "' in this contrast"
        ))
    }
    if (design_rank >= nrow(design)) {
        stop(paste0(
            "Blocking variables ", paste(blocking.vars, collapse = ', '), " leave no residual degrees of freedom in this contrast"
        ))
    }
}

colnames(design) <- make.names(colnames(design))
cat("Design matrix columns:\n   ", paste(colnames(design), collapse = ", "), "\n")

if (!is.null(opt\$contrast_string)) {
    contrast_string <- as.character(opt\$contrast_string)
} else {
    treatment_target <- make.names(paste0(contrast_variable, ".", opt\$target_level))
    treatment_reference <- make.names(paste0(contrast_variable, ".", opt\$reference_level))

    if (!all(c(treatment_target, treatment_reference) %in% colnames(design))) {
        stop(paste0(treatment_target, " and/or ", treatment_reference, " not found in design matrix"))
    }
    contrast_string <- paste0(treatment_target, "-", treatment_reference)
}
cat("Using contrast string:", contrast_string, "\n")
contrast.matrix <- makeContrasts(contrasts = contrast_string, levels = colnames(design))

################################################
################################################
## Run edgeR                                  ##
################################################
################################################

dge <- DGEList(counts = as.matrix(count.table), samples = sample.sheet, genes = feature.annotation)

if (opt\$filter_by_expr) {
    keep <- filterByExpr(dge, design = design, min.count = opt\$min_count)
    cat("Keeping", sum(keep), "of", length(keep), "features after filterByExpr\n")
    dge <- dge[keep, , keep.lib.sizes = FALSE]
}

dge <- normLibSizes(dge, method = opt\$norm_method)
dge <- estimateDisp(dge, design = design, robust = opt\$robust)

if (opt\$test_method == 'qlf') {
    fit <- glmQLFit(dge, design = design, robust = opt\$robust)
} else {
    fit <- glmFit(dge, design = design)
}

splice <- diffSpliceDGE(
    fit,
    contrast = contrast.matrix,
    geneid = opt\$gene_id_col,
    exonid = opt\$feature_id_col
)

exon.results  <- topSpliceDGE(splice, test = "exon",  number = Inf, FDR = 1)
simes.results <- topSpliceDGE(splice, test = "Simes", number = Inf, FDR = 1)
gene.results  <- topSpliceDGE(splice, test = "gene",  number = Inf, FDR = 1)

# topSpliceDGE adjusts with BH, so other methods are applied to its unadjusted p-values
if (opt\$adjust_method != "BH") {
    exon.results\$FDR  <- p.adjust(exon.results\$P.Value,  method = opt\$adjust_method)
    simes.results\$FDR <- p.adjust(simes.results\$P.Value, method = opt\$adjust_method)
    gene.results\$FDR  <- p.adjust(gene.results\$P.Value,  method = opt\$adjust_method)
}

################################################
################################################
## Generate outputs                           ##
################################################
################################################

write_results_table(round_numeric(exon.results, opt\$round_digits),  paste(opt\$output_prefix, 'edger.diffsplice.results.tsv', sep = '.'))
write_results_table(round_numeric(simes.results, opt\$round_digits), paste(opt\$output_prefix, 'edger.diffsplice.simes.tsv', sep = '.'))
write_results_table(round_numeric(gene.results, opt\$round_digits),  paste(opt\$output_prefix, 'edger.diffsplice.gene.tsv',  sep = '.'))

pdf(file = paste(opt\$output_prefix, 'edger.diffsplice.pdf', sep = '.'))
plot_genes <- head(as.character(simes.results[[opt\$gene_id_col]]), opt\$plot_n)
for (gene in plot_genes) {
    plotSpliceDGE(splice, geneid = gene, genecolname = opt\$gene_id_col)
}
if (nrow(simes.results) == 0) {
    plot.new()
    title(main = "No genes with more than one feature were tested")
} else if (length(plot_genes) == 0) {
    plot.new()
    title(main = "No genes plotted (plot_n is 0)")
}
invisible(dev.off())

png(
    file = paste(opt\$output_prefix, 'edger.bcv.png', sep = '.'),
    width = 600,
    height = 600
)
plotBCV(dge)
invisible(dev.off())

saveRDS(fit, file = paste(opt\$output_prefix, 'DGEGLM.edger.rds', sep = '.'))
saveRDS(splice, file = paste(opt\$output_prefix, 'DGELRT.edger.rds', sep = '.'))

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
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    bioconductor-edger:', as.character(packageVersion('edgeR'))),
        paste('    bioconductor-limma:', as.character(packageVersion('limma')))
    ),
    'versions.yml'
)
