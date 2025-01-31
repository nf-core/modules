#!/usr/bin/env Rscript


# Functions ---------------------------------------------------------------

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



# PARSE PARAMETERS FROM NEXTFLOW ------------------------------------------

# I've defined these in a single array like this so that we could go back to an
# optparse-driven method in future with module bin/ directories, rather than
# the template

# Set defaults and classes

opt <- list(
  output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
  diff_results = '$differential_result',
  pkn_network = "$pkn_network",
  feature_id_col = "gene_id",
  reference_level = "$meta.reference",
  treatment_level = "$meta.target",
  fold_change_col = "logFC",
  p_value_column = 'P.Value',
  diff_feature_id_col = "ID",
  fold_change_threshold = 2,
  p_value_threshold = 0.05,
  unlog_foldchanges = TRUE,
  palette_name = 'Set1'
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
  if (! ao %in% names(opt)){
    stop(paste("Invalid option:", ao))
  }
  else{

    # Preserve classes from defaults where possible
    if (! is.null(opt[[ao]])){
      args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
    }
    opt[[ao]] <- args_opt[[ao]]
  }
}

# Check if required parameters have been provided

required_opts <- c("diff_results", "pkn_network")
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]

if (length(missing) > 0){
  stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid

for (file_input in c("diff_results", "pkn_network")){
  if (is.null(opt[[file_input]])) {
    stop(paste("Please provide", file_input), call. = FALSE)
  }

  if (! file.exists(opt[[file_input]])){
    stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
  }
}


# Finish loading libraries ------------------------------------------------


library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)


# READ IN Limma differential results --------------------------------------
diff.table <-
  read_delim_flexible(
    file = opt$diff_results,
    header = TRUE,
    row.names = opt$diff_feature_id_col,
    check.names = FALSE
  ) %>%
  rownames_to_column(var = opt$diff_feature_id_col)

# Load Regulons -----------------------------------------------------------
# net <- get_collectri(organism='human', split_complexes=FALSE)
net = read.csv(opt$pkn_network, na.strings = c("", "NA"))

# Prepare DEGs for ULM ----------------------------------------------------
deg <- diff.table %>%
  select(ID, logFC, t, P.Value) %>%
  filter(!is.na(t)) %>%
  column_to_rownames(var = "ID") %>%
  as.matrix()
# Run ULM -----------------------------------------------------------------
contrast_acts <- run_ulm(mat=deg[, 't', drop=FALSE], net=net, .source='source', .target='target',
                         .mor='mor', minsize = 5)
# Generate outputs --------------------------------------------------------
contrast.name <- paste(opt$target_level, opt$reference_level, sep = "_vs_")
cat("Saving results for ", contrast.name, " ...\n", sep = "")

# Differential expression table- note very limited rounding for consistency of
# results

out_df <- contrast_acts
write.table(
  out_df,
  file = paste(opt$output_prefix, 'dorothea.results.tsv', sep = '.'),
  col.names = TRUE,
  row.names = FALSE,
  sep = '\t',
  quote = FALSE
)

# R SESSION INFO  ---------------------------------------------------------

sink(paste(opt$output_prefix, "R_sessionInfo.log", sep = '.'))
print(sessionInfo())
sink()


# VERSIONS FILE -----------------------------------------------------------
r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
decoupleR.version <- as.character(packageVersion('decoupleR'))

writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r.version),
    paste('    bioconductor-decoupler:', decoupleR.version)
  ),
  'versions.yml')


