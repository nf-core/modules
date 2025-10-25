#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)
library(readr)

#Load prefix
prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')

# Load the TSV file and keep only $gene_symbol_col column + counts
gene_expression_matrix <- readr::read_tsv('$input_file') %>%
    as.data.frame() %>%
    dplyr::select(
        dplyr::all_of('$gene_symbol_col'), # Keep the '$gene_symbol_col' column
        where(~ !is.character(.))    # Include all non-string columns
    ) %>%
    tibble::column_to_rownames('$gene_symbol_col')

# Check if the data is log-transformed or TPM-transformed
# Check range of values
min_value <- min(gene_expression_matrix, na.rm = TRUE)
max_value <- max(gene_expression_matrix, na.rm = TRUE)

# Detect log-transformed data (values typically within a compressed range)
if (max_value < 100 && min_value >= 0) {
    warning("The data appears to be log-transformed. Please provide TPM-transformed data.")
}

# Detect TPM-transformed data (column sums should be close to 1,000,000)
column_sums <- colSums(gene_expression_matrix, na.rm = TRUE)
if (!all(abs(column_sums - 1e6) < 1e4)) {
    warning("The data does not appear to be properly TPM-transformed. Ensure the data is normalized.")
}

# Generate results
result <- immunedeconv::${function}(gene_expression_matrix, method = '$method')

# Save the result to a CSV file
readr::write_tsv(result, paste0(prefix,'.deconvolution_results.tsv'))

# Plot and save results
# Plot 1: Stacked bar chart
plot1 <- result %>%
    gather(sample, fraction, -cell_type) %>%
    ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
    geom_bar(stat = 'identity') +
    coord_flip() +
    scale_fill_brewer(palette = 'Paired') +
    scale_x_discrete(limits = rev(levels(result)))

# Save Plot 1
ggsave(paste0(prefix,'.plot1_stacked_bar_chart.png'), plot = plot1, dpi = 300, width = 10, height = 8)

# Plot 2: Points with facets
plot2 <- result %>%
    gather(sample, score, -cell_type) %>%
    ggplot(aes(x = sample, y = score, color = cell_type)) +
    geom_point(size = 4) +
    facet_wrap(~cell_type, scales = 'free_x', ncol = 3) +
    scale_color_brewer(palette = 'Paired', guide = FALSE) +
    coord_flip() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Save Plot 2
ggsave(paste0(prefix,'.plot2_points_with_facets.png'), plot = plot2, dpi = 300, width = 12, height = 10)

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

immunedeconv_version <- as.character(packageVersion('immunedeconv'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-immunedeconv:', immunedeconv_version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################