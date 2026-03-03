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

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

# Set defaults and classes
opt <- list(
    output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    processed_vcf_rds = '$processed_vcf_rds',
    verbose = TRUE,
    add_ratios = TRUE,
    hmm = NULL,
    field_DP = NULL,
    cpus = $task.cpus
)

opt_types <- lapply(opt, class)

# Apply parameter overrides
args_opt <- parse_args('$task.ext.args')
for (ao in names(args_opt)){
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

# Check if required parameters have been provided
required_opts <- c('processed_vcf_rds', 'output_prefix')
missing <- required_opts[!unlist(lapply(opt[required_opts], is_valid_string)) | !required_opts %in% names(opt)]

if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check required file inputs are valid
for (file_input in c('processed_vcf_rds')){
    if (! is_valid_string(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }

    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

# Check optional file inputs and load if provided
if (!is.null(opt\$hmm)) {
    if (!file.exists(opt\$hmm)) {
        stop(paste0('HMM file does not exist: ', opt\$hmm))
    }
    
    hmm <- readRDS(opt\$hmm)
    
    # Validate HMM structure
    required_fields <- c("States", "Symbols", "startProbs", "transProbs", "emissionProbs")
    missing_fields <- setdiff(required_fields, names(hmm))
    
    if (length(missing_fields) > 0) {
        stop(paste0('Invalid HMM structure. Missing fields: ', paste(missing_fields, collapse = ", ")))
    }
} else {
    hmm <- NULL
}

################################################
################################################
## LOAD LIBRARIES                            ##
################################################
################################################

suppressPackageStartupMessages({
    library(UPDhmm)
    library(BiocParallel)
    library(parallel)
})

################################################
################################################
## LOAD PROCESSED VCF                        ##
################################################
################################################

processedVcf <- readRDS(opt\$processed_vcf_rds)

################################################
################################################
## SETUP PARALLEL PROCESSING                 ##
################################################
################################################

# Configure BiocParallel backend
if (opt\$cpus > 1) {
    bp_param <- MulticoreParam(workers = opt\$cpus)
} else {
    bp_param <- SerialParam()
}

################################################
################################################
## CALCULATE UPD EVENTS                      ##
################################################
################################################

updEvents <- tryCatch({
    calculateEvents(
        processedVcf, 
        hmm = hmm,
        field_DP = opt\$field_DP,
        add_ratios = opt\$add_ratios,
        BPPARAM = bp_param, 
        verbose = opt\$verbose
    )
}, error = function(e) {
    cat("[ERROR] ", conditionMessage(e), "\\n")
    quit(status = 1)
})


################################################
################################################
## SAVE OUTPUTS                              ##
################################################
################################################

# Save as text file
results_file <- paste0(opt\$output_prefix, ".upd_events.txt")

write.table(
    updEvents, 
    file = results_file, 
    sep = "\\t", 
    row.names = FALSE, 
    quote = FALSE
)


# Save as RDS file
rds_file <- paste0(opt\$output_prefix, ".upd_events.rds")

saveRDS(updEvents, file = rds_file)


################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

session_info_file <- paste0(opt\$output_prefix, ".R_sessionInfo.log")

sink(session_info_file)
sessionInfo()
sink()

