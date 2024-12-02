#!/usr/bin/env Rscript
library(optparse)
library(dplyr)
library(ggplot2)
library(fs)
library(NACHO)
library(readr)
library(tidyr)

# Parse Arguments
norm_methods <- c("GLM", "GEO")
option_list <- list(
    make_option(
        c("--input_rcc_path"),
        type = "character",
        default = "./" ,
        help = "Path to the folder that contains the RCC input file(s)",
        metavar = "character"),
    make_option(
        c("--input_samplesheet"),
        type = "character",
        default = NULL ,
        help = "Path to the sample sheet file",
        metavar = "character"),
    make_option(
        c("--norm_method"),
        type = "character",
        default = "GLM",
        help = paste0("Normalization method. One of ", paste(norm_methods, collapse = " "), paste = " "),
        metavar = "character")
)

# Parse the command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Validate mandatory arguments
if (is.null(opt$input_rcc_path)) {
    stop("Error: The --input_rcc_path parameter is mandatory and must be specified.")
}

if (is.null(opt$input_samplesheet)) {
    stop("Error: The --input_samplesheet parameter is mandatory and must be specified.")
}

# Validate that --norm_method is one of the allowed values
if (!(opt$norm_method %in% norm_methods)) {
    stop(paste("Error: The --norm_method parameter must be one of:", paste(norm_methods, collapse = " ")))
}

input_rcc_path    <- opt$input_rcc_path
input_samplesheet <- opt$input_samplesheet
norm_method       <- opt$norm_method

# Create filelist for NachoQC

list_of_rccs <- dir_ls(path = input_rcc_path, glob = "*.RCC")
print(list_of_rccs)

# Core Code
## Read data
nacho_data <- load_rcc(data_directory = input_rcc_path,
                    ssheet_csv = input_samplesheet,
                    id_colname = "RCC_FILE_NAME",
                    normalisation_method = norm_method)

output_base <- "./"

get_counts <- function(
    nacho,
    codeclass = "Endogenous",
    rownames = "RCC_FILE_NAME",
    colnames = c("Name", "Accession")
) {
    nacho[["nacho"]] %>%
    dplyr::select(c("RCC_FILE_NAME", "Name", "Count_Norm", "CodeClass")) %>%
    tidyr::pivot_wider(names_from = "RCC_FILE_NAME", values_from = "Count_Norm")
}

## Write out normalized counts
norm_counts <- as.data.frame(get_counts(nacho_data))
write_tsv(norm_counts, file = "normalized_counts.tsv")

## Create non-hk normalized counts too
nacho_data_no_hk <- load_rcc(data_directory = input_rcc_path,
    ssheet_csv = input_samplesheet,
    id_colname = "RCC_FILE_NAME",
    normalisation_method = norm_method,
    housekeeping_norm = FALSE)

## Export non-hk tables
norm_counts_without_hks <- as.data.frame(get_counts(nacho_data_no_hk))
write_tsv(norm_counts_without_hks, file = "normalized_counts_wo_HKnorm.tsv")
