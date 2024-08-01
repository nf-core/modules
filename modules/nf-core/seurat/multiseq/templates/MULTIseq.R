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
    assay = "HTO",
    quantile = 0.7,
    autoThresh = FALSE,
    maxiter = 5,
    qrange_from = 0.1,
    qrange_to = 0.9,
    qrange_by = 0.05,
    verbose = TRUE,
    selection_method = "mean.var.plot",
    normalization_method = "CLR",
    group_cells_feature_scatter = "HTO_maxID",
    produce_feature_scatter = TRUE,
    feature_scatter_feature_1 = "",
    feature_scatter_feature_2 = "",
    produce_ridge_plot = TRUE,
    number_of_features_ridge_plot = 2,
    produce_violin_plot = TRUE,
    group_cells_violin_plot = "HTO_classification.global",
    features_violin_plot = "nCount_RNA",
    pt_size = 0.1,
    log = TRUE,
    produce_tsne_plot = TRUE,
    subset_idents = "Negative",
    subset_invert = TRUE,
    tsne_scale_data_verbose = FALSE,
    run_pca_approx = FALSE,
    run_tsne_dim_max = 2,
    run_tsne_perplexity = 100
)

opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }else{

        # Preserve classes from defaults where possible
        if (! is.null(opt[[ao]])){
            args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        }
        opt[[ao]] <- args_opt[[ao]]
    }
}


################################################
################################################
## Create Seurat object                       ##
################################################
################################################
#### Seurat object creation following tutorial from: https://satijalab.org/seurat/articles/hashing_vignette

# Read data
rna_mtx <- Read10X("/Users/mylenemarianagonzalesandre/Development/data/rna/filtered_feature_bc_matrix")

hto_mtx <- Read10X("/Users/mylenemarianagonzalesandre/Development/data/hto/filtered_feature_bc_matrix")

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
seurat_object <- FindVariableFeatures(seurat_object, selection.method = opt$selection_method, nfeatures = 2000)
seurat_object <- ScaleData(seurat_object, features = VariableFeatures(seurat_object))


# Add HTO data as a new assay independent from RNA
seurat_object[[opt$assay]] <- CreateAssayObject(counts = hto_mtx)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
seurat_object <- NormalizeData(seurat_object, assay = opt$assay, normalization.method = opt$normalization_method)


################################################
## MULTIseq                       ##
################################################
demuliplex <- MULTIseqDemux(seurat_object ,
  assay = opt$assay,
  quantile = opt$quantile,
  autoThresh = opt$autoThresh,
  maxiter = opt$maxiter,
  qrange = seq(from = opt$qrange_from, to = opt$qrange_to, by = opt$qrange_by),
  verbose = TRUE
)
################################################
## Results                                    ##
################################################
write.csv(seurat_object$MULTI_ID, paste0(opt$output_prefix, "_assignment.csv"))