#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(fs)
library(NACHO)
library(readr)
library(tidyr)
library(nfcore.utils)

################################################
## Functions                                  ##
################################################

get_counts <- function(
    nacho,
    codeclass = "Endogenous",
    rownames = "RCC_FILE_NAME",
    colnames = c("Name", "Accession")
) {
    nacho[["nacho"]] |>
        dplyr::select(c("RCC_FILE_NAME", "Name", "Count_Norm", "CodeClass")) |>
        tidyr::pivot_wider(names_from = "RCC_FILE_NAME", values_from = "Count_Norm")
}

################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################

opt <- list(
    output_prefix = "${prefix}",
    input_rcc_path = "input",
    input_samplesheet = "${sample_sheet}",
    norm_method = "GLM"
)

opt_valid <- process_inputs(
    opt,
    args = '${args}',
    keys_to_nullify = c("norm_method", "output_prefix"),
    expected_folder = c("input_rcc_path"),
    expected_files = c("input_samplesheet"),
    required_opts = c("input_rcc_path", "input_samplesheet", "output_prefix")
)

# Validate that --norm_method is one of the allowed values
norm_methods <- c("GLM", "GEO")
if (!(opt[["norm_method"]] %in% norm_methods)) {
    stop(paste("Error: The --norm_method parameter must be one of:", paste(norm_methods, collapse = " ")))
}

input_rcc_path    <- opt[["input_rcc_path"]]
input_samplesheet <- opt[["input_samplesheet"]]
norm_method       <- opt[["norm_method"]]
prefix            <- opt[["output_prefix"]]

# Create filelist for NachoQC

list_of_rccs <- dir_ls(path = input_rcc_path, glob = "*.RCC")
print(list_of_rccs)

# Core Code
## Read data
nacho_data <- load_rcc(
    data_directory = input_rcc_path,
    ssheet_csv = input_samplesheet,
    id_colname = "RCC_FILE_NAME",
    normalisation_method = norm_method
)

output_base <- "./"

## Write out normalized counts
norm_counts <- as.data.frame(get_counts(nacho_data))
write_tsv(norm_counts, file = paste0(prefix, ".tsv"))

## Create non-hk normalized counts too
nacho_data_no_hk <- load_rcc(
    data_directory = input_rcc_path,
    ssheet_csv = input_samplesheet,
    id_colname = "RCC_FILE_NAME",
    normalisation_method = norm_method,
    housekeeping_norm = FALSE
)

## Export non-hk tables
norm_counts_without_hks <- as.data.frame(get_counts(nacho_data_no_hk))
write_tsv(norm_counts_without_hks, file = paste0(prefix, "_wo_HKnorm.tsv"))

process_end(
    packages = list(
        "r-nacho" = "NACHO",
        "r-dplyr" = "dplyr",
        "r-ggplot2" = "ggplot2",
        "r-tidyr" = "tidyr",
        "r-readr" = "readr",
        "r-fs" = "fs"
    ),
    task_name = "${task.process}",
    versions_path = "versions.yml",
    log_path = "${prefix}.R_sessionInfo.log"
)
