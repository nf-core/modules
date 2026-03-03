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
    upd_events_rds = '$upd_events_rds',
    min_ME = 2,
    min_size = 500000
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
required_opts <- c('upd_events_rds', 'output_prefix')
missing <- required_opts[!unlist(lapply(opt[required_opts], is_valid_string)) | !required_opts %in% names(opt)]

if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid
for (file_input in c('upd_events_rds')){
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
})

################################################
################################################
## LOAD UPD EVENTS                           ##
################################################
################################################

updEvents <- readRDS(opt\$upd_events_rds)

################################################
################################################
## COLLAPSE UPD EVENTS                       ##
################################################
################################################

collapsed <- tryCatch({
    collapseEvents(
      updEvents, 
      min_ME = opt\$min_ME, 
      min_size = opt\$min_size
    )
}, error = function(e) {
    cat("[ERROR] collapseEvents failed\\n")
    cat("[ERROR] ", conditionMessage(e), "\\n")
    quit(status = 1)
})

################################################
################################################
## SAVE OUTPUTS                              ##
################################################
################################################

# Save as text file
collapsed_file <- paste0(opt\$output_prefix, ".upd_collapsed.txt")

write.table(
      collapsed, 
      file = collapsed_file, 
      sep = "\\t", 
      row.names = FALSE, 
      quote = FALSE
  )


# Save as RDS file
collapsed_rds_file <- paste0(opt\$output_prefix, ".upd_collapsed.rds")

saveRDS(collapsed, file = collapsed_rds_file)


################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

session_info_file <- paste0(opt\$output_prefix, ".R_sessionInfo.log")

sink(session_info_file)
sessionInfo()
sink()
