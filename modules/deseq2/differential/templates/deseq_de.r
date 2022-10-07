#!/usr/bin/env Rscript

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

# I've defined these in a single array like this so that we could go back to an
# optparse-driven method in future with module bin/ directories, rather than
# the template

opt <- list(
    'count_file' = '$counts',
    'sample_file' = '$samplesheet',
    'contrast_variable' = '$meta.variable',
    'reference_level' = '$meta.reference',
    'treatment_level' = '$meta.target',
    'blocking_variables' = '$meta.blocking',
    'gene_id_col' = '$params.deseq_gene_id_col',
    'sample_id_col' = '$params.deseq_sample_id_col',
    'test' = '$params.deseq_test',
    'fit_type' = '$params.deseq_fit_type',
    'sf_type' = '$params.deseq_sf_type',
    'min_replicates_for_replace' = $params.deseq_min_replicates_for_replace,
    'use_t' = $params.deseq_use_t,
    'lfc_threshold' = $params.deseq_lfc_threshold,
    'alt_hypothesis' = '$params.deseq_alt_hypothesis',
    'independent_filtering' = $params.deseq_independent_filtering,
    'p_adjust_method' = '$params.deseq_p_adjust_method',
    'alpha' = $params.deseq_alpha,
    'minmu' = $params.deseq_minmu,
    'write_normalised' = $params.deseq_write_normalised,
    'write_variance_stabilised' = $params.deseq_write_variance_stabilised,
    'vs_method' = '$params.deseq_vs_method',
    'outdir' = '.'
)

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

contrast.name <-
    paste(opt\$treatment_level, opt\$reference_level, sep = "_vs_")
cat("Saving results for ", contrast.name, " ...\n", sep = "")

if (!file.exists(opt\$outdir)) {
    dir.create(opt\$outdir, recursive = TRUE)
}

# Differential expression table

write.table(
    data.frame(gene_id = rownames(comp.results), comp.results),
    file = file.path(opt\$outdir, 'deseq2.results.tsv'),
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

# Dispersion plot

png(
    file.path(opt\$outdir, 'deseq2.dispersion.png'),
    width = 600,
    height = 600
)
plotDispEsts(dds)
dev.off()

# R object for other processes to use

save(dds, file = file.path(opt\$outdir, 'dds.rld.RData'))

# Size factors

sf_df = data.frame(sample = names(sizeFactors(dds)), data.frame(sizeFactors(dds)))
colnames(sf_df) <- c('sample', 'sizeFactor')
write.table(
    sf_df,
    file = file.path(opt\$outdir, 'deseq2.sizefactors.tsv'),
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

# Write specified matrices

if (opt\$write_normalised){
    write.table(
        data.frame(gene_id=rownames(counts(dds)), counts(dds, normalized = TRUE)),
        file = file.path(opt\$outdir, 'normalised_counts.tsv'),
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
        data.frame(gene_id=rownames(counts(dds)), assay(vs_func(dds))),
        file = file.path(opt\$outdir, 'variance_stabilised_counts.tsv'),
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

RLogFile <- file.path(opt\$outdir, "R_sessionInfo.log")

sink(RLogFile)
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
