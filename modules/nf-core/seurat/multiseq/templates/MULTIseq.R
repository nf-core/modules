#!/usr/bin/env Rscript

# Tool from Seurat released under the MIT license.
# Ported to nf-core/modules with template by Mariana Gonzales Andre

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
    hto_matrix = '$hto_matrix',
    rna_matrix = '$rna_matrix',
    produce_plots = '$produce_plots'
)

opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')

# All values from ext.args are stored as strings
# Function to transform strings to the correct class
convert_element <- function(x) {
  if (is.character(x)) {
    # Try to convert to numeric
    num_value <- suppressWarnings(as.numeric(x))
    if (!is.na(num_value)) {
      return(num_value)
    }
    
    # Try to convert to boolean
    bool_value <- tolower(x)
    if (bool_value %in% c("true", "t", "yes", "y", "1")) {
      return(TRUE)
    } else if (bool_value %in% c("false", "f", "no", "n", "0")) {
      return(FALSE)
    }
  }
  
  # If no conversion was possible, return the original value
  return(x)
}
process_list <- function(input_list) {
  return(lapply(input_list, convert_element))
}
opt_args_transformed <- process_list(args_opt)

################################################
################################################
## Create Seurat object                       ##
################################################
################################################
#### Seurat object creation following tutorial from: https://satijalab.org/seurat/articles/hashing_vignette
library(Seurat)

# Read data
rna_mtx <- Read10X(opt\$rna_matrix)

hto_mtx <- Read10X(opt\$hto_matrix)

# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
joint_matrices <- intersect(colnames(rna_mtx), colnames(hto_mtx))

# Subset RNA and HTO counts by joint cell barcodes
rna_mtx <- rna_mtx[, joint_matrices]
hto_mtx <- as.matrix(hto_mtx[, joint_matrices])

# Confirm that the HTO have the correct names
rownames(hto_mtx)

# Setup Seurat object
seurat_object <- CreateSeuratObject(counts = Matrix::Matrix(as.matrix(rna_mtx), sparse = T))
# Normalize RNA data with log normalization
seurat_object <- NormalizeData(seurat_object)
# Find and scale variable features
seurat_object <- FindVariableFeatures(seurat_object, selection.method = opt_args_transformed\$selection_method, nfeatures = opt_args_transformed\$nfeatures)
seurat_object <- ScaleData(seurat_object, features = VariableFeatures(seurat_object))
# Add HTO data as a new assay independent from RNA
seurat_object[[opt_args_transformed\$assay]] <- CreateAssayObject(counts = hto_mtx)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
seurat_object <- NormalizeData(seurat_object, assay = opt_args_transformed\$assay, normalization.method = opt_args_transformed\$normalization_method)
################################################
## MULTIseq                       ##
################################################
demultiplex <- MULTIseqDemux(seurat_object ,
  assay = opt_args_transformed\$assay,
  quantile = opt_args_transformed\$quantile,
  autoThresh = opt_args_transformed\$autoThresh,
  maxiter = opt_args_transformed\$maxiter,
  qrange = seq(from = opt_args_transformed\$qrange_from, to = opt_args_transformed\$qrange_to, by = opt_args_transformed\$qrange_by),
  verbose = TRUE
)

################################################
## Results                                    ##
################################################
write.csv(demultiplex\$MULTI_ID, paste0(opt\$output_prefix, "_assignment.csv"))