#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(anndataR)
    library(Seurat)
    library(Matrix)
    library(SoupX)
})

set.seed(0)

# Template-interpolated by main.nf
h5ad_filtered     <- "${h5ad_filtered}"
h5ad_raw          <- "${h5ad_raw}"
prefix            <- "${prefix}"
npcs              <- ${npcs}
cluster_algorithm <- ${cluster_algorithm}

adata     <- read_h5ad(h5ad_filtered)
adata_raw <- read_h5ad(h5ad_raw)

# anndata stores cells x genes; Seurat/SoupX expect genes x cells
counts_filt <- t(adata\$X)
counts_raw  <- t(adata_raw\$X)

if (sum(counts_filt) == 0) stop("Filtered counts matrix is empty - cannot run SoupX.")

srt <- CreateSeuratObject(counts = counts_filt)
srt <- NormalizeData(srt, normalization.method = "LogNormalize", verbose = FALSE)
srt <- FindVariableFeatures(srt, verbose = FALSE)
srt <- ScaleData(srt, verbose = FALSE)
srt <- RunPCA(srt, npcs = npcs, verbose = FALSE)
srt <- FindNeighbors(srt, dims = 1:npcs, verbose = FALSE)
srt <- FindClusters(srt, algorithm = cluster_algorithm, verbose = FALSE)

soupx_groups <- setNames(as.character(Idents(srt)), colnames(srt))
rm(srt); gc()

sc <- SoupChannel(counts_raw, counts_filt, calcSoupProfile = FALSE)
soup_prof <- data.frame(
    row.names = rownames(counts_filt),
    est       = rowSums(counts_filt) / sum(counts_filt),
    counts    = rowSums(counts_filt)
)
sc  <- setSoupProfile(sc, soup_prof)
sc  <- setClusters(sc, soupx_groups)
sc  <- autoEstCont(sc, doPlot = FALSE)
out <- adjustCounts(sc, roundToInt = FALSE)

adata\$layers[["ambient"]] <- t(out)
write_h5ad(adata, paste0(prefix, ".h5ad"))

# versions.yml
versions <- c(
    '"${task.process}":',
    paste0("    r-base: ",     as.character(getRversion())),
    paste0("    soupx: ",      as.character(packageVersion("SoupX"))),
    paste0("    anndatar: ",   as.character(packageVersion("anndataR"))),
    paste0("    seurat: ",     as.character(packageVersion("Seurat"))),
    paste0("    leidenbase: ", as.character(packageVersion("leidenbase")))
)
writeLines(versions, "versions.yml")
