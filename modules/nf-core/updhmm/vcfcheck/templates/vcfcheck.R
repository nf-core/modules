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
    vcf = '$vcf',
    tbi = '$tbi',
    genome_build = 'hg38',
    proband_id = '$meta.proband_id',    
    mother_id = '$meta.mother_id',      
    father_id = '$meta.father_id'       
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
required_opts <- c('vcf', 'output_prefix', 'genome_build', 'proband_id', 'mother_id', 'father_id')
missing <- required_opts[!unlist(lapply(opt[required_opts], is_valid_string)) | !required_opts %in% names(opt)]

if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid
for (file_input in c('vcf')){
    if (! is_valid_string(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }

    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

################################################
################################################
## LOAD LIBRARIES                            ##
################################################
################################################

suppressPackageStartupMessages({
    library(UPDhmm)
    library(VariantAnnotation)
})

################################################
################################################
## READ VCF FILE                             ##
################################################
################################################

vcf <- readVcf(opt\$vcf, opt\$genome_build)


################################################
################################################
## RUN VCF CHECK                             ##
################################################
################################################

tryCatch({
    processedVcf <- vcfCheck(vcf, 
                           proband = opt\$proband_id, 
                           mother = opt\$mother_id, 
                           father = opt\$father_id)
}, error = function(e) {
    cat("ERROR in vcfCheck:", conditionMessage(e), "\\n")
    quit(status = 1)
})


################################################
################################################
## SAVE OUTPUTS                              ##
################################################
################################################

processed_rds_file <- paste0(opt\$output_prefix, ".processed.rds")
saveRDS(processedVcf, file = processed_rds_file)
  

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

session_info_file <- paste0(opt\$output_prefix, ".R_sessionInfo.log")

sink(session_info_file)
print(sessionInfo())
sink()


################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
updhmm.version <- as.character(packageVersion('UPDhmm'))
variantannotation.version <- as.character(packageVersion('VariantAnnotation'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-updhmm:', updhmm.version),
        paste('    bioconductor-variantannotation:', variantannotation.version)
    ),
    'versions.yml'
)