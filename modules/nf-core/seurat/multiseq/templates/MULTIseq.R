#!/usr/bin/env Rscript

# Tool from Seurat released under the MIT license.
# Ported to nf-core/modules with template by Mariana Gonzales Andre

#' Check for Non-Empty, Non-Whitespace String
#'
#' This function checks if the input is non-NULL and contains more than just whitespace.
#' It returns TRUE if the input is a non-empty, non-whitespace string, and FALSE otherwise.
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

#' All values from ext.args are stored as strings
#' Function to transform strings to the correct class
convert_element <- function(x) {
    if (is.character(x)) {
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
    return(x)
}
#' apply the conversion to all elements of a list
process_list <- function(input_list) {
    return(lapply(input_list, convert_element))
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

#' Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
opt_args_transformed <- process_list(args_opt)
opt <- process_list(opt)

################################################
################################################
## Create Seurat object                       ##
################################################
################################################
#### Seurat object creation following tutorial from: https://satijalab.org/seurat/articles/hashing_vignette
library(Seurat)
library(ggplot2)

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
demultiplex <- MULTIseqDemux(seurat_object,assay = opt_args_transformed\$assay, quantile = opt_args_transformed\$quantile, autoThresh = opt_args_transformed\$autoThresh, maxiter = opt_args_transformed\$maxiter, qrange = seq(from = opt_args_transformed\$qrange_from, to = opt_args_transformed\$qrange_to, by = opt_args_transformed\$qrange_by), verbose = opt_args_transformed\$verbose)
################################################
################################################
## Generate plots                             ##
################################################
################################################

if(opt\$produce_plots){
    ### Ridge plot
    Idents(demultiplex) <- opt_args_transformed\$group_cells_feature_scatter
    RidgePlot(demultiplex, assay = opt_args_transformed\$assay , features = rownames(demultiplex[[opt_args_transformed\$assay ]])[1:opt_args_transformed\$number_of_features_ridge_plot], ncol = opt_args_transformed\$number_of_cols_ridge_plot)
    ggsave(paste0(opt\$output_prefix, '_ridge_plot.jpeg'), device = "jpeg", dpi = 500)

    ### Feature scatter plot
    FeatureScatter(demultiplex, feature1 = opt_args_transformed\$feature_scatter_feature_1, feature2 = opt_args_transformed\$feature_scatter_feature_2)
    ggsave(paste0(opt\$output_prefix, '_feature_scatter_plot.jpeg'), device = "jpeg", dpi = 500)

    ### Violin plot
    Idents(demultiplex) <- opt_args_transformed\$group_cells_violin_plot
    VlnPlot(demultiplex, features = opt_args_transformed\$features_violin_plot, pt.size = opt_args_transformed\$pt_size, log = opt_args_transformed\$log)
    ggsave(paste0(opt\$output_prefix, '_violin_plot.jpeg'), device = "jpeg", dpi = 500)

    ### TSNE plot for HTOs
    subset_demultiplexing_results <- subset(demultiplex, idents = opt_args_transformed\$subset_idents, invert = opt_args_transformed\$subset_invert)
    DefaultAssay(subset_demultiplexing_results) <- opt_args_transformed\$assay
    subset_demultiplexing_results <- ScaleData(subset_demultiplexing_results, features = rownames(subset_demultiplexing_results),
    verbose = opt_args_transformed\$tsne_scale_data_verbose)
    subset_demultiplexing_results <- RunPCA(subset_demultiplexing_results, features =  rownames(subset_demultiplexing_results), approx = opt_args_transformed\$run_pca_approx)
    subset_demultiplexing_results <- RunTSNE(subset_demultiplexing_results, dims = 1:opt_args_transformed\$run_tsne_dim_max, perplexity = opt_args_transformed\$run_tsne_perplexity, check_duplicates = opt_args_transformed\$check_duplicates_tsne)
    DimPlot(subset_demultiplexing_results)
    ggsave(paste0(opt\$output_prefix, '_tsne_plot.jpeg'), device = "jpeg", dpi = 500)

    ### TSNE plot classification
    hto_names <- rownames(demultiplex[[opt_args_transformed\$assay]])
    # Extract the singlets
    subset_demultiplexing_singlet <- subset(demultiplex, idents = hto_names)
    # Select the top 1000 most variable features
    subset_demultiplexing_singlet <- FindVariableFeatures(subset_demultiplexing_singlet, selection.method = opt_args_transformed\$selection_method)
    # Scaling RNA data, we only scale the variable features here for efficiency
    subset_demultiplexing_singlet <- ScaleData(subset_demultiplexing_singlet, features = VariableFeatures(subset_demultiplexing_singlet))
    # Run PCA
    subset_demultiplexing_singlet <- RunPCA(subset_demultiplexing_singlet, features = VariableFeatures(subset_demultiplexing_singlet),check_duplicates = opt_args_transformed\$check_duplicates_tsne)

    # We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
    subset_demultiplexing_singlet <- FindNeighbors(subset_demultiplexing_singlet, reduction = "pca", dims = 1:opt_args_transformed\$run_tsne_dim_max)
    subset_demultiplexing_singlet <- FindClusters(subset_demultiplexing_singlet, resolution = opt_args_transformed\$resolution, verbose = opt_args_transformed\$tsne_scale_data_verbose)
    subset_demultiplexing_singlet <- RunTSNE(subset_demultiplexing_singlet, reduction = "pca", dims = 1:opt_args_transformed\$run_tsne_dim_max,check_duplicates = opt_args_transformed\$check_duplicates_tsne)
    DimPlot(subset_demultiplexing_singlet, group.by = opt_args_transformed\$singlet_identities_tsne)
    ggsave(paste0(opt\$output_prefix, '_tsne_classification.jpeg'), device = "jpeg", dpi = 500)
}


################################################
## Results                                    ##
################################################
write.csv(demultiplex\$MULTI_classification, paste0(opt\$output_prefix, "_classification.csv"))
write.csv(demultiplex\$MULTI_ID, paste0(opt\$output_prefix, "_assignment.csv"))

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(opt\$output_prefix, "R_sessionInfo.log", sep = '.'))
print(sessionInfo())
sink()

################################################
## Version                                    ##
################################################
r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
seurat_object_version <- as.character(packageVersion(('Seurat')))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    Seurat Object:', seurat_object_version)
    ),
'versions.yml')
