#!/usr/bin/env Rscript

################################################
################################################
## Functions                                  ##
################################################
################################################
cat("Initializing functions\n")

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = TRUE){

    ## Get file extension
    ext <- tools::file_ext(basename(file))

    ## Define separator
    if (ext == "tsv" || ext == "txt") {
        separator <- '\t'
    } else if (ext == "csv") {
        separator <- ','
    } else {
        stop(paste("Unknown separator for", ext))
    }

    ## Read file
    cat("Reading file", basename(file), "with", ext, "separator\n")

    df <- read.delim(
        file,
        sep = separator,
        header = header,
        row.names = row.names,
        check.names = check.names
    )

    return(df)
}

#' Check if a formula contains mixed model components (random effects)
#'
#' This function checks if a given formula contains random effects (e.g., terms like `(1|group)`)
#' by internally using the `findbars()` function to identify random effect terms.
#' Duplicated from `variancePartition:::.isMixedModelFormula()`, to avoid relying on an internal function
#'
#' @param formula A formula object. For example, `~ condition + (1|group)`.
#'
#' @return A logical value: `TRUE` if the formula contains random effects, otherwise `FALSE`.
#'
#' @details The function checks for the presence of random effects in the formula by looking for
#'          terms that include a `|` symbol, which indicates random effects in a mixed model.
#'          The `findbars()` function is defined internally to identify these terms.
#'
#' @examples
#' # Example formula with random effects
#' form <- ~ 0 + condition + (1|group)
#' is_mixed <- my_isMixedModelFormula(form)
#' print(is_mixed)  # Output: TRUE
#'
#' # Example formula without random effects
#' form_no_random <- ~ 0 + condition
#' is_mixed_no_random <- my_isMixedModelFormula(form_no_random)
#' print(is_mixed_no_random)  # Output: FALSE

isMixedModelFormula <- function(formula) {

    !is.null(lme4::findbars(as.formula(formula)))

}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################
library(optparse)
cat("Parsing arguments\n")

# Define the full list of options
option_list <- list(
    make_option(c("-o", "--output_prefix"), type = "character", default = "dream_analysis",
        help = "Prefix for output files [default:  %default]"),
    make_option(c("-c", "--count_file"), type = "character", default = NULL,
        help = "File containing raw counts [default:  %default]"),
    make_option(c("-s", "--sample_file"), type = "character", default = NULL,
        help = "File containing sample information [default:  %default]"),
    make_option(c("--contrast_variable"), type = "character", default = NULL,
        help = "Variable for contrast [default:  %default]"),
    make_option(c("--reference_level"), type = "character", default = NULL,
        help = "Reference level for the contrast [default:  %default]"),
    make_option(c("--target_level"), type = "character", default = NULL,
        help = "Target level for the contrast [default:  %default]"),
    make_option(c("--blocking_variables"), type = "character", default = NULL,
        help = "Blocking variables for the analysis [default: %default]"),
    make_option(c("--sample_id_col"), type = "character", default = "sample",
        help = "Column name for sample IDs [default: %default]"),
    make_option(c("--subset_to_contrast_samples"), type = "logical", default = FALSE,
        help = "Whether to subset to contrast samples [default: %default]"),
    make_option(c("--exclude_samples_col"), type = "character", default = NULL,
        help = "Column for excluding samples [default: %default]"),
    make_option(c("--exclude_samples_values"), type = "character", default = NULL,
        help = "Values for excluding samples [default: %default]"),
    make_option(c("--threads"), type = "numeric", default = 1,
        help = "Number of threads for multithreading [default: %default]"),
    make_option(c("--ddf"), type = "character", default = "adaptive",
        help = "Specifiy 'Satterthwaite', 'Kenward-Roger', or 'adaptative' method for dream() [default: %default]"),
    make_option(c("--reml"), type = "logical", default = TRUE,
        help = "Use restricted maximum likelihood to fit linear mixed model with dream() [default: %default]"),
    make_option(c("--proportion"), type = "numeric", default = 0.01,
        help = "Proportion for eBayes [default: %default]"),
    make_option(c("--stdev_coef_lim"), type = "character", default = "0.1,4",
        help = "Standard deviation coefficient limits for eBayes [default: %default]"),
    make_option(c("--trend"), type = "logical", default = FALSE,
        help = "Whether to use trend in eBayes [default: %default]"),
    make_option(c("--robust"), type = "logical", default = FALSE,
        help = "Whether to use robust method in eBayes [default: %default]"),
    make_option(c("--winsor_tail_p"), type = "character", default = "0.05,0.1",
        help = "Winsor tail probabilities for eBayes [default: %default]"),
    make_option(c("--adjust_method"), type = "character", default = "BH",
        help = "Adjustment method for topTable [default: %default]"),
    make_option(c("--p_value"), type = "numeric", default = 1,
        help = "P-value threshold for topTable [default: %default]"),
    make_option(c("--lfc"), type = "numeric", default = 0,
        help = "Log fold-change threshold for topTable [default: %default]"),
    make_option(c("--confint"), type = "logical", default = FALSE,
        help = "Whether to compute confidence intervals in topTable [default: %default]")
)

# Parse options
description <-
"This script performs differential gene expression analysis using the R 'dream' package from 'variancePartition'.\n\n
It requires the following inputs:
    1. A counts table (gene expression data)
    2. A variable, reference and contrast levels (for specifying the contrasts in the analysis)
    3. A sample sheet (sample metadata and experimental setup)."

opt <- parse_args(OptionParser(option_list = option_list, description = description))

# Check if required parameters have been provided
cat("Validating required arguments\n")

required_opts <- c("contrast_variable", "reference_level", "target_level", "output_prefix", "count_file", "sample_file")
missing <- required_opts[sapply(required_opts, function(o) is.null(opt[[o]]))]

if (length(missing) > 0) {
    stop(paste("Missing required arguments:", paste(missing, collapse = ", ")), call. = FALSE)
}

# Check if required file inputs are valid
for (file_input in c("count_file", "sample_file")) {
    if (!file.exists(opt[[file_input]])) {
        stop(paste0("Value of ", file_input, ": ", opt[[file_input]], " is not a valid file"), call. = FALSE)
    }
}

## Check default values for ddf options
ddf_valid <- c("Satterthwaite", "Kenward-Roger", "adaptive")
if ( !opt$ddf %in% ddf_valid ) {
    stop(paste0("'--ddf '", opt$ddf, "' is not a valid option from '", paste(ddf_valid, collapse = "', '"), "'"), call. = FALSE)
}

## Check adjust method options
adjust_valid <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
if (!opt$adjust_method %in% adjust_valid && !is.null(opt$adjust_method)) {
        stop(paste0("'--adjust_method '", opt$adjust_method, "' is not a valid option from '", paste(adjust_valid, collapse = "', '"), "', or NULL"), call. = FALSE)
}

# Convert specific options to numeric vectors
vector_opt <- c("stdev_coef_lim", "winsor_tail_p")
opt[vector_opt] <- lapply(strsplit(unlist(opt[vector_opt]), ","), as.numeric)

# Save first version of RData, useful for debuging with all original parameters already set
cat("Exporting preliminary RData\n")

work_dir <- getwd()                         ## for dev purposes
save.image("dream_de.RData")
#setwd(work_dir); load("dream_de.RData")    ## for dev purposes

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################
cat("Importing libraries\n")

library(edgeR)
library(variancePartition)

################################################
################################################
## READ IN COUNTS FILE AND SAMPLE METADATA    ##
################################################
################################################

intensities.table <- read_delim_flexible(file = opt$count_file)
sample.sheet <- read_delim_flexible(file = opt$sample_file)

# Deal with spaces that may be in sample column
opt$sample_id_col <- make.names(opt$sample_id_col)

if (! opt$sample_id_col %in% colnames(sample.sheet)){
    stop(paste0("Specified sample ID column '", opt$sample_id_col, "' is not in the sample sheet"))
}

# Sample sheet can have duplicate rows for multiple sequencing runs, so uniqify
# before assigning row names
sample.sheet <- sample.sheet[! duplicated(sample.sheet[[opt$sample_id_col]]), ]
rownames(sample.sheet) <- sample.sheet[[opt$sample_id_col]]

# Check that all samples specified in the input sheet are present in the
# intensities table. Assuming they are, subset and sort the count table to
# match the sample sheet

missing_samples <-
    sample.sheet[!rownames(sample.sheet) %in% colnames(intensities.table), opt$sample_id_col]

if (length(missing_samples) > 0) {
    stop(paste(
        length(missing_samples),
        'specified samples missing from count table:',
        paste(missing_samples, collapse = ',')
    ))
} else {
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
cat("Validating contrasts\n")

contrast_variable <- make.names(opt$contrast_variable)
blocking.vars <- c()

if (!contrast_variable %in% colnames(sample.sheet)) {
    stop(
        paste0(
        'Chosen contrast variable "',
        contrast_variable,
        '" not in sample sheet'
        )
    )
} else if (any(!c(opt$reference_level, opt$target_level) %in% sample.sheet[[contrast_variable]])) {
    stop(
        paste0(
        'Please choose reference and target levels that are present in the ',
        contrast_variable,
        ' column of the sample sheet'
        )
    )
} else if (!is.null(opt$blocking_variables)) {
    blocking.vars = make.names(unlist(strsplit(opt$blocking_variables, split = ';')))
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

# Handle conflicts between blocking variables and block
if (!is.null(opt$block) && !is.null(opt$blocking_variables)) {
    if (opt$block %in% blocking.vars) {
        warning(paste("Variable", opt$block, "is specified both as a fixed effect and a random effect. It will be treated as a random effect only."))
        blocking.vars <- setdiff(blocking.vars, opt$block)
        if (length(blocking.vars) == 0) {
            opt$blocking_variables <- NULL
        } else {
            opt$blocking_variables <- paste(blocking.vars, collapse = ';')
        }
    }
}

# Optionally, subset to only the samples involved in the contrast

if (opt$subset_to_contrast_samples){
    sample_selector <- sample.sheet[[contrast_variable]] %in% c(opt$target_level, opt$reference_level)
    selected_samples <- sample.sheet[sample_selector, opt$sample_id_col]
    intensities.table <- intensities.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
}

# Optionally, remove samples with specified values in a given field (probably
# don't use this as well as the above)

if ((! is.null(opt$exclude_samples_col)) && (! is.null(opt$exclude_samples_values))){
    exclude_values = unlist(strsplit(opt$exclude_samples_values, split = ';'))

    if (! opt$exclude_samples_col %in% colnames(sample.sheet)){
        stop(paste(opt$exclude_samples_col, ' specified to subset samples is not a valid sample sheet column'))
    }

    print(paste0('Excluding samples with values of ', opt$exclude_samples_values, ' in ', opt$exclude_samples_col))
    sample_selector <- ! sample.sheet[[opt$exclude_samples_col]] %in% exclude_values

    selected_samples <- sample.sheet[sample_selector, opt$sample_id_col]
    intensities.table <- intensities.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
}

################################################
################################################
## Build the Model Formula                    ##
################################################
################################################

## TODO: START ADAPTING TO DREAM IN HERE!
## NEW STARTS HERE!
cat("Creating formula\n")

# Build the model formula with blocking variables first
model_vars <- c()

if (!is.null(opt$blocking_variables)) {
    cat("opt$blocking_variables:", paste(opt$blocking_variables, collapse = ' '), "\n")
    print(opt$blocking_variables)

    # Include blocking variables (including pairing variables if any)
    for (VARIABLE in opt$blocking_variables) {
        cat("Adding", VARIABLE, "factor to formula\n")
        ## Create a string to reconstruct Wilkinson formula " (1 | variable )"
        model_vars <- c(model_vars, paste0("(1 | ", VARIABLE, ")"))

        ## Convert the variable into factor
        if (!is.numeric(sample.sheet[[ VARIABLE ]])) {
            cat("Converting into factor\n")
            sample.sheet[[ VARIABLE ]] <- as.factor(sample.sheet[[ VARIABLE ]])
        }
    }
}

# Construct the model formula
## Expected structure:
## "~ 0 + fixed_effect + (1 | random_variable_1) + (1 | random_variable_N)"
model <- paste(
    '~ 0 +',
    contrast_variable,
    if (length(model_vars) > 0) { paste0(" +", paste(model_vars, collapse = ' + ')) } else { NULL },
    sep = "")  ## TODO: This is limited to additive models! not possible for interaction relations

cat("Model vars:", model_vars, "\n")
cat("Model:", model, "\n")

# Construct the formula
form <- as.formula(model)
cat("Formula:", deparse(form), "\n")

## TODO: Check if the model is mixed or not, could be useful to report it later
cat("Checking for mixed formula\n")
mixed_form <- isMixedModelFormula(form)

################################################
################################################
## Run Dream processes                        ##
################################################
################################################

# Generate the design matrix
cat("Creating design matrix\n")

design <- model.matrix(
    form,
    sample.sheet
)

# Specify parallel processing
param <- SnowParam(as.numeric(opt$threads), "SOCK", progressbar = TRUE)

# Create a DGEList object for RNA-seq data
cat("Creating DGEList\n")
dge <- DGEList(counts = intensities.table)

## Calculate normalization factors
cat("Calculating normalization factors\n")
dge <- calcNormFactors(dge)

# estimate weights using linear mixed model of dream
cat("Normalizing data\n")
vobjDream <- voomWithDreamWeights(dge, form, sample.sheet, BPPARAM = param)

# Create and export variance plot
cat("Analyzing variance\n")
vp <- fitExtractVarPartModel(exprObj = intensities.table, formula = update(form, ~ . - 0), data = sample.sheet, REML = opt$reml)

cat("Creating variance plot\n")
var_plot <- plotVarPart(sortCols(vp))

cat("Exporting variance plot\n")
png(
    file = paste(opt$output_prefix, 'dream.var_plot.png', sep = '.'),
    width = 600,
    height = 300
)
plot(var_plot)
dev.off()

## Set contrast (this can be scaled for more than one comparison)
cat("Building contrasts\n")
L <- variancePartition::makeContrastsDream(
    form,
    sample.sheet,
    contrasts = c(                                                      ## This is a named vector for all the contrast we'd like to run.
        setNames(                                                       ## For now, it's just one contrast but we should be able to scale it with the yml somehow
            paste(                                                      ## For a variable named `treatment`, with levels A and B
                paste0(opt$contrast_variable, opt$target_level),            ## Create target value: treatmentB
                paste0(opt$contrast_variable, opt$reference_level),         ## Create reference value: treatmentA
                sep = " - "),                                               ## Set a difference between them
            opt$output_prefix)                                          ## Assign the name to the comparison (`contrast id`` from the yml)
        )
    )

# Visualize contrast matrix
cat("Creating and exporting contrast plot\n")
contrasts_plot <- plotContrasts(L)

# Export contrast plot
png(
    file = paste(opt$output_prefix, 'dream.contrasts_plot.png', sep = '.'),
    width = 600,
    height = 300
)
plot(contrasts_plot)
dev.off()


# Fit the dream model on each gene
# For the hypothesis testing, by default,
# dream() uses the KR method for <= 20 samples,
# otherwise it uses the Satterthwaite approximation (this is the behavior indicated with `ddf = adaptive`)
# We can force it to use the others methods with is ddf argument
cat("Fitting model with dream()\n")
fitmm <-
    dream(
        exprObj = vobjDream,
        formula = form,
        data    = sample.sheet,
        L       = L,
        ddf     = opt$ddf,
        reml    = opt$reml
    )

# Adjust results with Empirical Bayes
## Create list of arguments for eBayes
cat("Adjusting results with eBayes\n")
ebayes_args <- list(
    fit = fitmm
)

if (! is.null(opt$proportion)){
    ebayes_args[['proportion']] <- as.numeric(opt$proportion)
}
if (! is.null(opt$stdev.coef.lim)){
    ebayes_args[['stdev.coef.lim']] <- as.numeric(opt$stdev.coef.lim)
}
if (! is.null(opt$trend)){
    ebayes_args[['trend']] <- as.logical(opt$trend)
}
if (! is.null(opt$robust)){
    ebayes_args[['robust']] <- as.logical(opt$robust)
}
if (! is.null(opt$winsor.tail.p)){
    ebayes_args[['winsor.tail.p']] <- as.numeric(opt$winsor.tail.p)
}

## Run variancePartition::eBayes
fitmm <- do.call(variancePartition::eBayes, ebayes_args)

# get names of available coefficients and contrasts for testing
colnames(fitmm)

# Get results of hypothesis test on coefficients of interest (only one coeff for now)
cat("Exporting results with topTable()\n")
for (COEFFICIENT in opt$output_prefix) {

    ## Initialize topTable() arguments
    toptable_args <- list(
        fit = fitmm,
        coef = COEFFICIENT,
        sort.by = 'none',
        number = nrow(intensities.table)
    )

    ## Complete list with extra arguments, if they were provided
    if (! is.null(opt$adjust.method)){
        toptable_args[['adjust.method']] <- opt$adjust_method
    }
    if (! is.null(opt$p.value)){
        toptable_args[['p.value']] <- as.numeric(opt$p_value)
    }
    if (! is.null(opt$lfc)){
        toptable_args[['lfc']] <- as.numeric(opt$lfc)
    }
    if (! is.null(opt$confint)){
        toptable_args[['confint']] <- as.logical(opt$confint)
    }

    ## generate topTable
    comp.results <- do.call(variancePartition::topTable, toptable_args)[rownames(intensities.table),]

    ## Export topTable
    write.table(
        comp.results,
        file = paste(opt$output_prefix, 'dream.results.tsv', sep = '.'),
        col.names = TRUE,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )
}

## END OF NEW

################################################
################################################
## Generate other outputs                     ##
################################################
################################################

# Dispersion plot

png(
    file = paste(opt$output_prefix, 'dream.mean_difference.png', sep = '.'),
    width = 600,
    height = 600
)
plotMD(fitmm)
dev.off()

# R object for other processes to use
saveRDS(fitmm, file = paste(opt$output_prefix, 'MArrayMM.dream.rds', sep = '.'))

# Save model to file
write(model, file=paste(opt$output_prefix, 'dream.model.txt', sep = '.'))

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(opt$output_prefix, "R_sessionInfo.log", sep = '.'))
print(sessionInfo())
sink()

################################################
################################################
################################################
################################################

save.image("dream_de.RData")
