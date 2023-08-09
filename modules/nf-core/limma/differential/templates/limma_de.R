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

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = TRUE){

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
    output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    count_file = '$intensities',
    sample_file = '$samplesheet',
    contrast_variable = '$contrast_variable',
    reference_level = '$reference',
    target_level = '$target',
    blocking_variables = NULL,
    probe_id_col = "probe_id",
    sample_id_col = "experiment_accession",
    subset_to_contrast_samples = FALSE,
    exclude_samples_col = NULL,
    exclude_samples_values = NULL,
    ndups = NULL,                # lmFit
    spacing = NULL,              # lmFit
    block = NULL,                # lmFit
    correlation = NULL,          # lmFit
    method = 'ls',               # lmFit
    proportion = 0.01,           # eBayes
    stdev.coef.lim = '0.1,4',    # eBayes
    trend = FALSE,               # eBayes
    robust = FALSE,              # eBayes
    winsor.tail.p = '0.05,0.1',  # eBayes
    adjust.method = "BH",        # topTable
    p.value = 1,                 # topTable
    lfc = 0,                     # topTable
    confint = FALSE              # topTable
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

# Check if required parameters have been provided

required_opts <- c('contrast_variable', 'reference_level', 'target_level', 'output_prefix')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]

if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
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

# Convert things to vectors where needed
vector_opt <- c('stdev.coef.lim', 'winsor.tail.p')
opt[vector_opt] = lapply(strsplit(unlist(opt[vector_opt]), ','), as.numeric)

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(limma)

################################################
################################################
## READ IN COUNTS FILE AND SAMPLE METADATA    ##
################################################
################################################

intensities.table <-
    read_delim_flexible(
        file = opt\$count_file,
        header = TRUE,
        row.names = opt\$probe_id_col,
        check.names = FALSE
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

# Check that all samples specified in the input sheet are present in the
# intensities table. Assuming they are, subset and sort the count table to
# match the sample sheet

missing_samples <-
    sample.sheet[!rownames(sample.sheet) %in% colnames(intensities.table), opt\$sample_id_col]

if (length(missing_samples) > 0) {
    stop(paste(
        length(missing_samples),
        'specified samples missing from count table:',
        paste(missing_samples, collapse = ',')
    ))
} else{
    # Save any non-count data, will gene metadata etc we might need later
    nonintensities.table <-
        intensities.table[, !colnames(intensities.table) %in% rownames(sample.sheet), drop = FALSE]
    intensities.table <- intensities.table[, rownames(sample.sheet)]
}

################################################
################################################
## CHECK CONTRAST SPECIFICATION               ##
################################################
################################################

contrast_variable <- make.names(opt\$contrast_variable)
blocking.vars <- c()

if (!contrast_variable %in% colnames(sample.sheet)) {
    stop(
        paste0(
        'Chosen contrast variable \"',
        contrast_variable,
        '\" not in sample sheet'
        )
    )
} else if (any(!c(opt\$reflevel, opt\$treatlevel) %in% sample.sheet[[contrast_variable]])) {
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
# Optionally, subset to only the samples involved in the contrast

if (opt\$subset_to_contrast_samples){
    sample_selector <- sample.sheet[[contrast_variable]] %in% c(opt\$target_level, opt\$reference_level)
    selected_samples <- sample.sheet[sample_selector, opt\$sample_id_col]
    intensities.table <- intensities.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
}

# Optionally, remove samples with specified values in a given field (probably
# don't use this as well as the above)

if ((! is.null(opt\$exclude_samples_col)) && (! is.null(opt\$exclude_samples_values))){
    exclude_values = unlist(strsplit(opt\$exclude_samples_values, split = ';'))

    if (! opt\$exclude_samples_col %in% colnames(sample.sheet)){
        stop(paste(opt\$exclude_samples_col, ' specified to subset samples is not a valid sample sheet column'))
    }

    print(paste0('Excluding samples with values of ', opt\$exclude_samples_values, ' in ', opt\$exclude_samples_col))
    sample_selector <- ! sample.sheet[[opt\$exclude_samples_col]] %in% exclude_values

    selected_samples <- sample.sheet[sample_selector, opt\$sample_id_col]
    intensities.table <- intensities.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
}

# Now specify the model. Use cell-means style so we can be explicit with the
# contrasts

model <- paste('~ 0 +', contrast_variable)

if (!is.null(opt\$blocking_variables)) {
    model <- paste(model, paste(blocking.vars, collapse = '+'), sep = '+')
}

# Make sure all the appropriate variables are factors

for (v in c(blocking.vars, contrast_variable)) {
    sample.sheet[[v]] <- as.factor(sample.sheet[[v]])
}

################################################
################################################
## Run Limma processes                       ##
################################################
################################################

# Generate the design

design <- model.matrix(
    as.formula(model),
    data=sample.sheet
)
colnames(design) <- sub(
    contrast_variable,
    paste0(contrast_variable, '.'), colnames(design)
)

# Prepare for and run lmFit()

lmfit_args = c(
    list(
        object = as.matrix(intensities.table),
        design = design
    ),
    opt[c('ndups', 'spacing', 'block', 'method')]
)

if (! is.null(opt\$block)){
    lmfit_args[['block']] = sample.sheet[[opt\$block]]
}
if (! is.null(opt\$correlation)){
    lmfit_args[['correlation']] = opt\$correlation
}

fit <- do.call(lmFit, lmfit_args)

# Contrasts bit
contrast <- paste(paste(contrast_variable, c(opt\$target_level, opt\$reference_level), sep='.'), collapse='-')
contrast.matrix <- makeContrasts(contrasts=contrast, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)

# Prepare for and run ebayes

ebayes_args = c(
    list(
        fit = fit2
    ),
    opt[c('proportion', 'stdev.coef.lim', 'trend', 'robust', 'winsor.tail.p')]
)

fit2 <- do.call(eBayes, ebayes_args)

# Run topTable() to generate a results data frame

toptable_args = c(
    list(
        fit = fit2,
        sort.by = 'none',
        number = nrow(intensities.table)
    ),
    opt[c('adjust.method', 'p.value', 'lfc', 'confint')]
)

comp.results <- do.call(topTable, toptable_args)[rownames(intensities.table),]

################################################
################################################
## Generate outputs                           ##
################################################
################################################

contrast.name <-
    paste(opt\$target_level, opt\$reference_level, sep = "_vs_")
cat("Saving results for ", contrast.name, " ...\n", sep = "")

# Differential expression table- note very limited rounding for consistency of
# results

write.table(
    data.frame(
        probe_id = rownames(comp.results),
        comp.results
    ),
    file = paste(opt\$output_prefix, 'limma.results.tsv', sep = '.'),
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

# Dispersion plot

png(
    file = paste(opt\$output_prefix, 'limma.mean_difference.png', sep = '.'),
    width = 600,
    height = 600
)
plotMD(intensities.table)
dev.off()

# R object for other processes to use

saveRDS(fit2, file = paste(opt\$output_prefix, 'MArrayLM.limma.rds', sep = '.'))

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
limma.version <- as.character(packageVersion('limma'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-limma:', limma.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
