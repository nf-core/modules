#!/usr/bin/env Rscript

################################################
################################################
## Functions                                  ##
################################################
################################################

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

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

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

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = TRUE){

    ext <- tolower(tail(strsplit(basename(file), split = "\\.")[[1]], 1))

    if (ext == "tsv" || ext == "txt") {
        separator <- "\t"
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

round_dataframe_columns <- function(df, columns = NULL, digits = 8){
    if (is.null(columns)){
        columns <- colnames(df)
    }

    df[,columns] <- format(
        data.frame(df[, columns], check.names = FALSE),
        nsmall = digits
    )

    # Convert columns back to numeric

    for (c in columns) {
        df[[c]][grep("^ *NA\\$", df[[c]])] <- NA
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
    output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    count_file = '$counts',
    sample_file = '$samplesheet',
    sample_treatment_col = '$sample_treatment_col',
    reference_level = '$reference',
    target_level = '$target',
    sample_id_col = "sample",
    samples_pairing_col = NULL,
    samples_batch_col = NULL,
    subset_to_contrast_samples = FALSE,
    exclude_samples_col = NULL,
    exclude_samples_values = NULL,
    dataType = "RNAseq",             # for anota2seqDataSetFromMatrix
    normalize = TRUE,                # for anota2seqDataSetFromMatrix
    transformation = "TMM-log2",     # for anota2seqDataSetFromMatrix
    filterZeroGenes = TRUE,          # for anota2seqDataSetFromMatrix
    varCutOff = NULL,                # for anota2seqDataSetFromMatrix
    performQC = TRUE,                # for anota2seqRun
    onlyGroup = FALSE,               # for anota2seqRun
    performROT = TRUE,               # for anota2seqRun
    generateSingleGenePlots = FALSE, # for anota2seqRun
    analyzeBuffering = TRUE,         # for anota2seqRun
    analyzemRNA = TRUE,              # for anota2seqRun
    useRVM = TRUE,                   # for anota2seqRun
    correctionMethod = "BH",         # for anota2seqRun
    useProgBar = FALSE,              # for anota2seqRun
    getRVM = TRUE,                   # for anota2seqGetOutput
    output = 'full'                  # for anota2seqGetOutput
)
opt_types <- lapply(opt, class)

# Apply parameter overrides
# Assuming task.ext.args is passed as a command line argument
args <- commandArgs(trailingOnly = TRUE)
task_ext_args <- paste(args, collapse = " ")

# Parse the arguments
args_opt <- parse_args(task_ext_args)
#args_opt <- parse_args('$task.ext.args')
print(args_opt)
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

required_opts <- c('sample_treatment_col', 'reference_level', 'target_level', 'output_prefix')
missing <- required_opts[!unlist(lapply(opt[required_opts], is_valid_string)) | !required_opts %in% names(opt)]

if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid

for (file_input in c('count_file', 'sample_file')){
    if (! is_valid_string(opt[[file_input]])) {
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

library(anota2seq)

################################################
################################################
## READ IN COUNTS FILE AND SAMPLE METADATA    ##
################################################
################################################

count.table <-
    read_delim_flexible(
        file = opt$count_file,
        header = TRUE,
        row.names = opt$gene_id_col,
        check.names = FALSE
    )
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

# Check that all samples specified in the input sheet are present in the counts
# table. Assuming they are, subset and sort the count table to match the sample
# sheet

missing_samples <-
    sample.sheet[!rownames(sample.sheet) %in% colnames(count.table), opt$sample_id_col]

if (length(missing_samples) > 0) {
    stop(paste(
        length(missing_samples),
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

sample_treatment_col <- make.names(opt$sample_treatment_col)

if (!sample_treatment_col %in% colnames(sample.sheet)) {
    stop(
        paste0(
        'Chosen contrast variable \"',
        sample_treatment_col,
        '\" not in sample sheet'
        )
    )
} else if (any(!c(opt$reflevel, opt$treatlevel) %in% sample.sheet[[sample_treatment_col]])) {
    stop(
        paste(
        'Please choose reference and treatment levels that are present in the',
        sample_treatment_col,
        'column of the sample sheet'
        )
    )
}

# Optionally, subset to only the samples involved in the contrast

if (opt$subset_to_contrast_samples){
    sample_selector <- sample.sheet[[sample_treatment_col]] %in% c(opt$target_level, opt$reference_level)
    selected_samples <- sample.sheet[sample_selector, opt$sample_id_col]
    count.table <- count.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
}

# Optionally, remove samples with specified values in a given field (probably
# don't use this as well as the above)

if ((is_valid_string(opt$exclude_samples_col)) && (is_valid_string(opt$exclude_samples_values))){
    exclude_values = unlist(strsplit(opt$exclude_samples_values, split = ';'))

    if (! opt$exclude_samples_col %in% colnames(sample.sheet)){
        stop(paste(opt$exclude_samples_col, ' specified to subset samples is not a valid sample sheet column'))
    }

    print(paste0('Excluding samples with values of ', opt$exclude_samples_values, ' in ', opt$exclude_samples_col))
    sample_selector <- ! sample.sheet[[opt$exclude_samples_col]] %in% exclude_values

    selected_samples <- sample.sheet[sample_selector, opt$sample_id_col]
    count.table <- count.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
}

################################################
################################################
## CHECK CONDITIONS AND REPLICATES            ##
################################################
################################################

# Count the number of conditions and replicates
condition_counts <- table(sample.sheet[[sample_treatment_col]][sample.sheet$type=="riboseq"]) # only count riboseq samples?
num_conditions <- length(condition_counts)
min_replicates <- min(condition_counts)


# Check conditions and decide on analysis
if (num_conditions >= 2 && min_replicates >= 3) {
    # Proceed with analysis as normal
    cat("Proceeding with full anota2seq analysis.\n")
    proceed_with_analysis <- TRUE
    multiple_contrasts <- FALSE # indicate whether to create a contrast matrix with multiple contrasts
} else if (num_conditions == 2 && min_replicates < 3) {
    # Skip analysis with warning
    cat("WARNING: anota2seq needs at least three replicates within each condition when there are only two conditions. Skipping anota2seq analysis.\n")
    proceed_with_analysis <- FALSE
} else if (num_conditions > 2 && min_replicates < 3) {
    # Run with onlyGroup = TRUE
    cat("WARNING: Less than three replicates detected in at least one condition. Running anota2seq with onlyGroup = TRUE.\n")
    opt$onlyGroup <- TRUE
    proceed_with_analysis <- TRUE
    multiple_contrasts <- TRUE
} else {
    # Unexpected case
    cat("ERROR: Unexpected combination of conditions and replicates. Please check your input data.\n")
    quit(status = 1)
}

# Proceed with analysis only if conditions are met
if (proceed_with_analysis) {
    
    ################################################
    ################################################
    ## Run anota2seqRun()                         ##
    ################################################
    ################################################

    # Separate matrix into riboseq and rnaseq data

    riboseq_samples <- sample.sheet[[opt$sample_id_col]][sample.sheet[['type']] == 'riboseq']
    rnaseq_samples <- sample.sheet[[opt$sample_id_col]][sample.sheet[['type']] == 'rnaseq']

    if (! is.null(opt$samples_pairing_col)){
        riboseq_samples <- riboseq_samples[order(sample.sheet[riboseq_samples, opt$samples_pairing_col])]
        rnaseq_samples <- rnaseq_samples[order(sample.sheet[rnaseq_samples, opt$samples_pairing_col])]
    }

    riboseq_data <- count.table[,riboseq_samples]
    rnaseq_data <- count.table[,rnaseq_samples]

    # Make the anota2seqDataSet

    anota2seqDataSetFromMatrix_args <- list(
        dataP = riboseq_data,
        dataT = rnaseq_data,
        phenoVec = sample.sheet[rnaseq_samples, opt$sample_treatment_col],
        dataType = opt$dataType,
        normalize = opt$normalize,
        transformation = opt$transformation,
        filterZeroGenes = opt$filterZeroGenes,
        varCutOff = opt$varCutOff
    )

    if (! is.null(opt$samples_batch_col)){
        anota2seqDataSetFromMatrix_args$batchVec <- samples[rnaseq_samples, opt$samples_batch_col]
    }

    ads <- do.call(anota2seqDataSetFromMatrix, anota2seqDataSetFromMatrix_args)

    # Run anota2seqRun
    if(multiple_contrasts == TRUE){
        # Create a contrast matrix with multiple contrasts
        contrast_matrix <- matrix(
            nrow=num_conditions,
            ncol=num_conditions - 1,
            dimnames=list(names(condition_counts),c()),
            0
        )
        # set the reference level to -1
        contrast_matrix[which(rownames(contrast_matrix)==opt$reference_level),] <- -1
        # set each target level to 1
        contrast_matrix[which(rownames(contrast_matrix)==opt$target_level),1] <- 1 # ganrantueed target_vs_reference to be the first contrast
        other_conditions <- names(condition_counts)[!names(condition_counts) %in% c(opt$reference_level, opt$target_level)]
        for (i in 1:length(other_conditions)){
            contrast_matrix[which(rownames(contrast_matrix)==other_conditions[i]),i+1] <- 1
        }
    } else {
        # Create a contrast matrix with a single contrast
        contrast_matrix <- matrix(
            nrow=2,
            ncol=1,
            dimnames=list(c(opt$reference_level, opt$target_level),c()),
            c(-1,1)
        )
    }

    ads <- anota2seqRun(
        ads,
        contrasts = contrast_matrix,
        performQC = opt$performQC,
        onlyGroup = opt$onlyGroup,
        performROT = opt$performROT,
        generateSingleGenePlots = opt$generateSingleGenePlots,
        analyzeBuffering = opt$analyzeBuffering,
        analyzemRNA = opt$analyzemRNA,
        useRVM = opt$useRVM,
        correctionMethod = opt$correctionMethod,
        useProgBar = FALSE
    )

    ################################################
    ################################################
    ## Generate outputs                           ##
    ################################################
    ################################################

    contrast.name <-
        paste(opt$target_level, opt$reference_level, sep = "_vs_")
    cat("Saving results for ", contrast.name, " ...\n", sep = "")

    for (analysis in c("translated mRNA", "total mRNA", "translation", "buffering", "mRNA abundance")){

        output <- anota2seqGetOutput(
            ads,
            analysis = analysis,
            output = opt$output,
            getRVM = opt$getRVM,
            selContrast = 1
        )

        write.table(
            output,
            file = paste(opt$output_prefix, sub(' ', '_', analysis), 'anota2seq.results.tsv', sep = '.'),
            col.names = TRUE,
            row.names = FALSE,
            sep = '\t',
            quote = FALSE
        )
    }

    # Fold change plot

    png(
        file = paste(opt$output_prefix, 'fold_change.png', sep = '.'),
        width = 720,
        height = 720
    )
    anota2seqPlotFC(
        ads,
        visualizeRegModes = "all",
        selContrast = 1,
        contrastName = contrast.name,
        plotToFile = FALSE
    )
    dev.off()

    # R object for other processes to use

    saveRDS(ads, file = paste(opt$output_prefix, 'Anota2seqDataSet.rds', sep = '.'))

    # Rename files files output by anota2seqPerformQC()

    sapply(list.files(pattern = "^ANOTA2SEQ_"), function(f) file.rename(f, sub("^ANOTA2SEQ_", paste0(opt$output_prefix, '.'), f)))

}

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(opt$output_prefix, "R_sessionInfo.log", sep = '.'))
# print warning messages dependent on conditions and replicates
if (num_conditions >= 2 && min_replicates >= 3) {
    # Proceed with analysis as normal
    print("Given at least three replicates, full anota2seq analysis was performed.")
} else if (num_conditions == 2 && min_replicates < 3) {
    # Skip analysis with warning
    print("WARNING: anota2seq needs at least three replicates within each condition when there are only two conditions. Skipping anota2seq analysis.")
} else if (num_conditions > 2 && min_replicates < 3) {
    # Run with onlyGroup = TRUE, this is suboptimal but worked
    print("WARNING: Less than three replicates detected in at least one condition. Running anota2seq with onlyGroup = TRUE.")
    print("In general, it is highly recommended by the authors of anota2seq to have at least three replicates in each condition.")
    print("However, the analysis can still be performed if you have less than three replicates but more than two conditions.")
    print("In this case, it is possible to suppress the omnibus interaction analysis and only perform the omnibus treatment analysis by setting onlyGroup = TRUE.")
    print("Please be aware that this analysis can reduce statistical power and the results should be interpreted with caution.")
    print("For any related questions, please refer to the anota2seq paper (https://doi.org/10.1093/nar/gkz223).")
}
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
anota2seq.version <- as.character(packageVersion('anota2seq'))

writeLines(
    c(
        '"${task.process}":',
        paste('    bioconductor-anota2seq:', anota2seq.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################