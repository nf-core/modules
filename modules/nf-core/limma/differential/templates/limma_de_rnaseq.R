#!/usr/bin/env Rscript
library(limma)
library(edgeR)
library(readr)
library(dplyr)
library(tibble)

library(IHW)

########## READ DATA
# setwd("/home/michal_zenczak/Projects/differentialabundance-pipeline/other/FINAL")

## Inputs for pairwise analysis
# counts <- "raw_counts_pairwise.tsv"
# samplesheet <- "samplesheet_pairwise.tsv"
# contrasts <- "contrasts_pairwise.tsv"

## Inputs for mixedmodel analysis
# counts <- "raw_counts_mixedmodel.tsv"
# samplesheet <- "samplesheet_mixedmodel.tsv"
# contrasts <- "contrasts_mixedmodel.tsv"

# contrast_table <- read_tsv(contrasts)
# contrast_names <- colnames(contrast_table)

# i=1
# contrast_variable <- contrast_table\$variable[i]
# target <- contrast_table\$target[i]
# reference <- contrast_table\$reference[i]
# block_var <- contrast_table\$block[i]

# prefix <- contrast_table\$id[i]
# sample_id_col <- "sampleId"
# limma_type <- "mixedmodel"
# limma_type <- "pairwise"

## Functions
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

# Validates the input parameters against the provided sample sheet.
# Ensures the contrast variable exists in the sample sheet,
# verifies that the reference and target are present within that variable,
# and checks if blocking variables, if provided, are valid columns in the sample sheet.
validate_sample_sheet <- function(sample_sheet, contrast_variable, reference_level, target_level, blocking_variables) {
  # Check if the contrast variable is in the sample sheet
  if (!contrast_variable %in% colnames(sample_sheet)) {
    stop(
      paste0(
        'Chosen contrast variable \"',
        contrast_variable,
        '\" not in sample sheet'
      )
    )
  }
  
  # Check if the reference and target levels are in the contrast variable column
  if (any(!c(reference_level, target_level) %in% sample_sheet[[contrast_variable]])) {
    stop(
      paste(
        'Please choose reference and target levels that are present in the',
        contrast_variable,
        'column of the sample sheet'
      )
    )
  }
  
  # If blocking variables are provided, check if they are in the sample sheet
  if (!is.null(blocking_variables)) {
    blocking_vars <- make.names(unlist(strsplit(blocking_variables, split = ';')))
    if (!all(blocking_vars %in% colnames(sample_sheet))) {
      missing_block <- paste(blocking_vars[!blocking_vars %in% colnames(sample_sheet)], collapse = ',')
      stop(
        paste(
          'Blocking variables', missing_block,
          'do not correspond to sample sheet columns.'
        )
      )
    }
  }
}

# Set defaults and classes
opt <- list(
  output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
  count_file = '$intensities',
  pair_id_col = 'pair_id',
  sample_file = '$samplesheet',
  contrast_variable = '$contrast_variable',
  reference_level = '$reference',
  target_level = '$target',
  blocking_variables = '$meta.blocking',
  probe_id_col = "probe_id",
  sample_id_col = "sample",
  analysis_type = NULL,
  subset_to_contrast_samples = FALSE,
  exclude_samples_col = NULL,
  exclude_samples_values = NULL,
  geneFilter = 10,
  IHW_correction = TRUE,
  do_contrast = TRUE,
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
  confint = FALSE,             # topTable
  number = Inf                 # topTable
)
opt_types <- lapply(opt, class)

args_opt <- parse_args('$task.ext.args')

# Update opt based on parsed_args
opt <- modifyList(opt, args_opt)

print(opt)

# Set correct types for options
for (name in names(opt)) {
  if (!is.null(opt[[name]])) {
    class(opt[[name]]) <- opt_types[[name]]
  }
}

## CHECKS
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

## READ IN COUNTS FILE AND SAMPLE METADATA 
count.table <-
  read_delim_flexible(
    file = opt\$count_file,
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
sample.sheet <- read_delim_flexible(file = opt\$sample_file)

# Correct sample.sheet colnames
colnames(sample.sheet) <- make.names(colnames(sample.sheet))

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
# counts table. Assuming they are, subset and sort the count table to
# match the sample sheet
missing_samples <-
  sample.sheet[!rownames(sample.sheet) %in% colnames(count.table), opt\$sample_id_col]

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

########## START ANALYSIS
DGE <- DGEList(count.table)
dim(DGE)

## Filter genes
keep <- which(apply(DGE\$counts, 1, max) >= opt\$geneFilter)
filtered_counts <- DGE\$counts[keep, ]
DGE\$counts <- filtered_counts
dim(DGE)

# Calculate normalization factors
DGE <- calcNormFactors(DGE)

if (opt\$analysis_type == "pairwise") {
  # Pairwise analysis setup
  contrast_variable <- make.names(opt\$contrast_variable)
  reference_level <- make.names(opt\$reference_level)
  target_level <- make.names(opt\$target_level)

  validate_sample_sheet(sample.sheet, contrast_variable, reference_level, target_level, opt\$blocking_variables)

  sample.sheet\$pair_IDs <- as.integer(as.factor(sample.sheet[[opt\$pair_id_col]]))
  model <- paste('~ 0 + pair_IDs +', contrast_variable)

  design <- model.matrix(as.formula(model), data=sample.sheet)
  colnames(design) <- sub(contrast_variable, paste0(contrast_variable, '.'), colnames(design))

  voom_args = list(counts = DGE, design = design, plot = FALSE)
  voom_1 <- do.call(voom, voom_args)
  # Write the normalized counts matrix to a TSV file
  normalized_counts <- voom_1\$E
  normalized_counts_with_genes <- data.frame(Gene = rownames(normalized_counts), normalized_counts, row.names = NULL)
  colnames(normalized_counts_with_genes)[1] <- opt\$probe_id_col
  write.table(normalized_counts_with_genes,
              file = paste(opt\$output_prefix, "normalised_counts.tsv", sep = '.'),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)

  lmfit_args = list(object = voom_1, design = design)
  fit <- do.call(lmFit, lmfit_args)

  if (opt\$do_contrast) {
    contrast <- paste(paste(contrast_variable, c(opt\$reference_level, opt\$target_level), sep='.'), collapse='-')
    contrast_matrix <- makeContrasts(contrasts=contrast, levels=design)
    fit <- contrasts.fit(fit, contrast_matrix)
  }

  ebayes_args = list(fit = fit)
  fit <- do.call(eBayes, ebayes_args)

} else if (opt\$analysis_type == "mixedmodel") {

  ## CHECK CONTRAST SPECIFICATION
  contrast_variable <- make.names(opt\$contrast_variable)
  reference_level <- make.names(opt\$reference_level)
  target_level <- make.names(opt\$target_level)
  sample.sheet[[contrast_variable]] = paste(sample.sheet[[strsplit2(contrast_variable, split = "\\\\.")[1]]],
                                            sample.sheet[[strsplit2(contrast_variable, split = "\\\\.")[2]]], sep = ".")

  validate_sample_sheet(sample.sheet, contrast_variable, reference_level, target_level, opt\$blocking_variables)

  blocking.vars <- sample.sheet[[opt\$blocking_variables]]

  # Mixed model analysis setup
  model <- paste('~ 0 +', contrast_variable)
  design <- model.matrix(as.formula(model), data=sample.sheet)
  colnames(design) <- sub(contrast_variable, paste0(contrast_variable, '.'), colnames(design))

  voom_args = list(counts = DGE, design = design, plot = FALSE)
  voom_1 <- do.call(voom, voom_args)
    # Write the normalized counts matrix to a TSV file
  normalized_counts <- voom_1\$E
  normalized_counts_with_genes <- data.frame(Gene = rownames(normalized_counts), normalized_counts, row.names = NULL)
  colnames(normalized_counts_with_genes)[1] <- opt\$probe_id_col
  write.table(normalized_counts_with_genes,
              file = paste(opt\$output_prefix, "normalised_counts.tsv", sep = '.'),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)

  corfit_args <- list(object = voom_1, design = design, block = blocking.vars)
  corfit = do.call(duplicateCorrelation, corfit_args)

  voom_args <- list(counts = DGE, design = design, plot = FALSE, correlation = corfit\$consensus.correlation)
  voom_1 <- do.call(voom, voom_args)

  corfit_args <- list(object = voom_1, design = design, block = blocking.vars)
  corfit = do.call(duplicateCorrelation, corfit_args)

  lmfit_args = list(object = voom_1, design = design, block = blocking.vars, correlation = corfit\$consensus.correlation)
  fit <- do.call(lmFit, lmfit_args)

  if (opt\$do_contrast) {
    contrast <- paste(paste(contrast_variable, c(reference_level, target_level), sep='.'), collapse='-')
    contrast_matrix <- makeContrasts(contrasts=contrast, levels=design)
    fit <- contrasts.fit(fit, contrast_matrix)
  }

  ebayes_args = list(fit = fit)
  fit <- do.call(eBayes, ebayes_args)

} else {
  stop(
    paste("Bad '--analysis_type' parameter. It can be either 'paired' or 'mixedmodel'. ")
  )
}

## Generate outputs
# Differential expression table

if (opt\$IHW_correction) {

  # Define the column name for the group comparison
  name_con_ref <- glue::glue("{contrast_variable}.{reference_level}-{contrast_variable}.{target_level}")

  # Calculate the median gene expression
  median_expression <- apply(DGE\$counts, 1, median)

  # Apply IHW for multiple testing correction using median expression as covariate
  ihw_result <- ihw(fit\$p.value[, name_con_ref] ~ median_expression, alpha = 0.05)
  ihw_result@df <- as_tibble(ihw_result@df, rownames = "Gene")

  # Convert fit\$p.value to a tibble, including rownames as a new column 'Gene'
  fit_df <- as_tibble(fit\$p.value, rownames = "Gene")

  fit_df <- fit_df %>%
    left_join(ihw_result@df %>% select(Gene,adj_pvalue),by="Gene") %>%
    rename(ihw.adj.p.value = adj_pvalue)

  # Convert back to a matrix if necessary (not needed for merging, but for consistency later)
  fit\$p.value <- as.matrix(fit_df %>% select(-Gene)) # Keep the original structure without 'Gene' column

  # Generate the topTable including the new adjusted p-values
  toptable_args <- list(fit = fit, coef = name_con_ref, number = opt\$number, adjust.method = opt\$adjust.method,
                        p.value = opt\$p.value, lfc = opt\$lfc, confint = opt\$confint)

  comp.results <- do.call(topTable,toptable_args)

  # Add the IHW-adjusted p-values to the topTable results by merging based on 'Gene'
  comp.results <- comp.results %>%
    rownames_to_column("Gene") %>%  # Ensure 'Gene' column exists in comp.results
    left_join(fit_df, by = "Gene") %>%  # Merge based on 'Gene' column
    select(-sym(name_con_ref))

  # Optional: Rearrange columns if necessary
  comp.results_test <- comp.results %>%
    select(Gene, everything())  # Move 'Gene' column to the front

} else {
   toptable_args <- list(fit = fit, number = opt\$number, adjust.method = opt\$adjust.method,
                        p.value = opt\$p.value, lfc = opt\$lfc, confint = opt\$confint)

   comp.results <- do.call(topTable,toptable_args)
}

names(comp.results)[names(comp.results) == 'Gene'] <- opt\$probe_id_col

write.table(
  comp.results,
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
plotMD(fit)
dev.off()

# R object for other processes to use
saveRDS(fit, file = paste(opt\$output_prefix, 'MArrayLM.limma.rds', sep = '.'))

# Save model to file
write(model, file = paste(opt\$output_prefix, 'limma.model.txt', sep = '.'))


## R SESSION INFO
sink(paste(opt\$output_prefix, "R_sessionInfo.log", sep = '.'))
print(sessionInfo())
sink()


## VERSIONS FILE
r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
limma.version <- as.character(packageVersion('limma'))

writeLines(
  c(
    '"task.process":',
    paste('    r-base:', r.version),
    paste('    bioconductor-limma:', limma.version)
  ),
  'versions.yml'
)
