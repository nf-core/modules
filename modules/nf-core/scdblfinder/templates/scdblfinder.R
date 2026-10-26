#!/usr/bin/env Rscript

library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
library(anndataR)

set.seed(123)

sce <- read_h5ad("${h5ad}")\$as_SingleCellExperiment()

# Expected doublet rate. If none is given, scDblFinder estimates it.
dbr <- suppressWarnings(as.numeric(trimws("${dbr}")))
if (is.na(dbr)) {
    dbr <- NULL
}

# Optional column that splits the cells by sample or capture
batch_col <- trimws("${batch_col ?: ''}")
samples <- NULL
if (nzchar(batch_col)) {
    if (!(batch_col %in% colnames(colData(sce)))) {
        stop(
            "Column '", batch_col, "' is not in the cell metadata. Available columns: ",
            paste(colnames(colData(sce)), collapse = ", ")
        )
    }
    samples <- as.vector(colData(sce)[[batch_col]])
    if (anyNA(samples)) {
        stop("Column '", batch_col, "' contains NA values.")
    }
}

workers <- max(1L, as.integer("${task.cpus}"))
args <- list(
    assay(sce, 1),
    BPPARAM = MulticoreParam(workers = workers, RNGseed = 123)
)
if (!is.null(dbr)) {
    args\$dbr <- dbr
}
if (!is.null(samples)) {
    args\$samples <- samples
}

res <- do.call(scDblFinder, args)

# Add the scDblFinder results to the original object so the other
# cell metadata is kept
cols <- grep("^scDblFinder", colnames(colData(res)), value = TRUE)
for (nm in cols) {
    sce[[nm]] <- res[[nm]]
}

message("scDblFinder results:")
print(table(sce\$scDblFinder.class))

write_h5ad(sce, "${prefix}.h5ad")
write.csv(as.data.frame(colData(sce))[, cols, drop = FALSE], "${prefix}.csv")

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
scdblfinder.version <- as.character(packageVersion('scDblFinder'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    scdblfinder:', scdblfinder.version)
    ),
    'versions.yml'
)
