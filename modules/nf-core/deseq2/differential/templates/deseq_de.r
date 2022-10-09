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
  args_vals <- unlist(lapply(args_list, function(y) strsplit(y, ' +')))
  
  as.list(structure(args_vals[c(FALSE, TRUE)], names = args_vals[c(TRUE, FALSE)]))
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
    count_file = '$counts',
    sample_file = '$samplesheet',
    contrast_variable = '$contrast_meta.variable',
    reference_level = '$contrast_meta.reference',
    treatment_level = '$contrast_meta.target',
    blocking_variables = '$contrast_meta.blocking',
    gene_id_col = "gene_id",
    sample_id_col = "experiment_accession",
    test = "Wald",
    fit_type = "parametric",
    sf_type = 'ratio',
    min_replicates_for_replace = 7,
    use_t = FALSE,
    lfc_threshold = 0,
    alt_hypothesis = 'greaterAbs',
    independent_filtering = 'TRUE',
    p_adjust_method = 'BH',
    alpha = 0.1,
    minmu = 0.5,
    write_normalised = TRUE,
    write_variance_stabilised = TRUE,
    vs_method = 'vst',
    random_seed = 0,
    round_results = FALSE
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }else{
        opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
    }
}

# Check file inputs are valid

for (file_input in c('count_file', 'sample_file')){
    if (is.null(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }

    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

if (opt\$random_seed > 0){
    set.seed(opt\$random_seed)
}


################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(DESeq2)

################################################
################################################
## READ IN COUNTS FILE AND SAMPLE METADATA    ##
################################################
################################################

count.table <-
    read.delim(
        file = opt\$count_file,
        header = TRUE,
        row.names = opt\$gene_id_col
    )
sample.sheet <- read.csv(file = opt\$sample_file)

if (! opt\$sample_id_col %in% colnames(sample.sheet)){
    stop(paste0("Specified sample ID column '", opt\$sample_id_col, "' is not in the sample sheet"))
}

# Sample sheet can have duplicate rows for multiple sequencing runs, so uniqify
# before assigning row names

sample.sheet <- sample.sheet[! duplicated(sample.sheet[[opt\$sample_id_col]]), ]
rownames(sample.sheet) <- sample.sheet[[opt\$sample_id_col]]

# Check that all samples specified in the input sheet are present in the counts
# table. Assuming they are, subset and sort the count table to match the sample
# sheet

missing_samples <-
    sample.sheet[!rownames(sample.sheet) %in% colnames(count.table), opt\$sample_id_col]

if (length(missing_samples) > 0) {
    stop(paste(
        len(missing_samples),
        'specified samples missing from count table:',
        paste(missing_samples, collapse = ',')
    ))
} else{
    # Save any non-count data, will gene metadata etc we might need later
    noncount.table <-
        count.table[, !colnames(count.table) %in% rownames(sample.sheet), drop = FALSE]
    count.table <- count.table[, rownames(sample.sheet)]
}

################################################
################################################
## CHECK CONTRAST SPECIFICATION               ##
################################################
################################################

blocking.vars <- c()

if (!opt\$contrast_variable %in% colnames(sample.sheet)) {
    stop(
        paste0(
        'Chosen contrast variable \"',
        opt\$contrast_variable,
        '\" not in sample sheet'
        )
    )
} else if (any(!c(opt\$reflevel, opt\$treatlevel) %in% sample.sheet[[opt\$contrast_variable]])) {
    stop(
        paste(
        'Please choose reference and treatment levels that are present in the',
        opt\$contrast_variable,
        'column of the sample sheet'
        )
    )
} else if (!is.null(opt\$blocking)) {
    blocking.vars = unlist(strsplit(opt\$blocking, split = ','))
    if (!all(blocking.vars %in% colnames(sample.sheet))) {
        stop(
            paste0(
                'One or more of the blocking variables specified (',
                opt\$blocking,
                ') do not correspond to sample sheet columns.'
            )
        )
    }
}

# Now specify the model. Use cell-means style so we can be explicit with the
# contrasts

model <- '~ 0 +'

if (!is.null(opt\$blocking_variables)) {
    model <- paste(model, paste(blocking.vars, collapse = '+'))
}

# Make sure all the appropriate variables are factors

for (v in c(blocking.vars, opt\$contrast_variable)) {
    sample.sheet[[v]] <- as.factor(sample.sheet[[v]])
}

# Variable of interest goes last, see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#multi-factor-designs

model <- paste(model, opt\$contrast_variable, sep = ' + ')

################################################
################################################
## Run DESeq2 processes                       ##
################################################
################################################

dds <- DESeqDataSetFromMatrix(
    countData = round(count.table),
    colData = sample.sheet,
    design = as.formula(model)
)

dds <- DESeq(
    dds,
    test = opt\$test,
    fitType = opt\$fit_type,
    minReplicatesForReplace = opt\$min_replicates_for_replace,
    useT = opt\$use_t,
    sfType = opt\$sf_type
)

comp.results <-
    results(
        dds,
        lfcThreshold = opt\$lfc_threshold,
        altHypothesis = opt\$alt_hypothesis,
        independentFiltering = opt\$independent_filtering,
        alpha = opt\$alpha,
        pAdjustMethod = opt\$p_adjust_method,
        minmu = opt\$minmu,
        contrast = c(
            opt\$contrast_variable,
            c(opt\$treatment_level, opt\$reference_level)
        )
    )

################################################
################################################
## Generate outputs                           ##
################################################
################################################

# Common function to round numerics in data frames - mostly useful for testing,
# since it doesn't seem to be possible to make outputs entirely reproducible
# across machienes, even with set.seed(), but the differences are tiny.

prepare_results <- function(x){
    if (opt\$round_results){
        format(x, nsmall=8)
    }else{
        x
    }
}


contrast.name <-
    paste(opt\$treatment_level, opt\$reference_level, sep = "_vs_")
cat("Saving results for ", contrast.name, " ...\n", sep = "")

# Differential expression table

write.table(
    prepare_results(data.frame(gene_id = rownames(comp.results), comp.results)),
    file = 'deseq2.results.tsv',
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

# Dispersion plot

png(
    'deseq2.dispersion.png',
    width = 600,
    height = 600
)
plotDispEsts(dds)
dev.off()

# R object for other processes to use

save(dds, file = 'dds.rld.RData')

# Size factors

sf_df = data.frame(sample = names(sizeFactors(dds)), data.frame(sizeFactors(dds)))
colnames(sf_df) <- c('sample', 'sizeFactor')
write.table(
    sf_df,
    file = 'deseq2.sizefactors.tsv',
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

# Write specified matrices

if (opt\$write_normalised){
    write.table(
        data.frame(gene_id=rownames(counts(dds)), counts(dds, normalized = TRUE)),
        file = 'normalised_counts.tsv',
        col.names = TRUE,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )
}

if (opt\$write_variance_stabilised){
    if (opt\$vs_method == 'vst'){
        vs_func = vst
    }else{
        vs_func = rlog
    }
    write.table(
        prepare_results(data.frame(gene_id=rownames(counts(dds)), assay(vs_func(dds)))),
        file = 'variance_stabilised_counts.tsv',
        col.names = TRUE,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )
}

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
deseq2.version <- as.character(packageVersion('DESeq2'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-deseq2:', deseq2.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
