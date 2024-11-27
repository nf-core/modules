#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(fs)
library(ComplexHeatmap)
library(circlize)
library(yaml)
library(ragg)

###Command line argument parsing###
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
    stop("Usage: compute_gene_heatmap.R <annotated_endo_data.tsv> <normalized_counts.tsv> <sample_id_col> or compute_gene_heatmap.R <annotated_endo_data.tsv> <normalized_counts.tsv> <heatmap_genes_to_filter.yaml> <sample_id_col>", call.=FALSE)
}
input_counts_annotated <- args[1]
input_counts <- args[2]
id_col <- tail(args, 1)

#Read annotated counts
# HEADER is always RCC_FILE + GENES + SAMPLE_ID and additional metadata such as GROUP TREATMENT OTHER_METADATA
counts <- read.table(input_counts_annotated, sep="\t", check.names = FALSE, header=TRUE, stringsAsFactors = FALSE)

if (length(args) == 4) {
    input_genes <- args[3]
    genes <- read_yaml(input_genes)
} else {
    genes <- read.table(input_counts, sep="\t", check.names = FALSE, header=TRUE, stringsAsFactors = FALSE) %>% dplyr::filter(CodeClass == 'Endogenous') %>% .$Name
}

#Select counts of interest
counts_selected <- counts %>% dplyr::select(all_of(genes))

#log2+1
counts_selected <- log2(counts_selected + 1)

#Find max
colMax <- function(data) sapply(data, max, na.rm = TRUE)
#Find min
colMin <- function(data) sapply(data, min, na.rm = TRUE)

max_value <- max(colMax(counts_selected))
min_value <- min(colMin(counts_selected))

#Save as PDF
prefix <- ""
if (grepl("wo_HKnorm", input_counts_annotated)) {
    prefix <- "wo_HKnorm_"
}

agg_png(file = paste0(prefix, "gene_heatmap_mqc.png"), width = 1200, height = 2000, unit = "px")

#Add proper row names
counts_matrix <- as.matrix(counts_selected)
row.names(counts_matrix) <- counts[[id_col]]

Heatmap(counts_matrix,
        name = "Gene-Count Heatmap",
        column_title = "Gene log2(count + 1)",
        row_order = order(row.names(counts_matrix)),
        row_title_rot = 90,
        row_title = "SampleID",
        row_dend_reorder = FALSE,
        show_row_dend = FALSE,
        row_names_side = "left",
        show_column_dend = FALSE,
        col = colorRamp2(c(min_value, max_value), c("#f7f7f7", "#67a9cf"))
    )

dev.off()
