#!/usr/bin/env Rscript

library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
library(anndataR)

set.seed(123)

sce <- read_h5ad("${h5ad}")\$as_SingleCellExperiment()

# Extra arguments for scDblFinder from task.ext.args, written as R named
# arguments, e.g. dbr = 0.08, samples = "batch", clusters = TRUE
extra <- list(${args})

# A single string given as 'samples' is read as the name of a column in the
# cell metadata (obs), so the cells can be split by sample or capture
if (is.character(extra\$samples) && length(extra\$samples) == 1) {
    batch_col <- extra\$samples
    if (!(batch_col %in% colnames(colData(sce)))) {
        stop(
            "Column '", batch_col, "' is not in the cell metadata. Available columns: ",
            paste(colnames(colData(sce)), collapse = ", ")
        )
    }
    extra\$samples <- as.vector(colData(sce)[[batch_col]])
    if (anyNA(extra\$samples)) {
        stop("Column '", batch_col, "' contains NA values.")
    }
}

workers <- max(1L, as.integer("${task.cpus}"))
res <- do.call(
    scDblFinder,
    c(
        list(assay(sce, 1), BPPARAM = MulticoreParam(workers = workers, RNGseed = 123)),
        extra
    )
)

# Add the scDblFinder results to the original object so the other
# cell metadata is kept
cols <- grep("^scDblFinder", colnames(colData(res)), value = TRUE)
for (nm in cols) {
    sce[[nm]] <- res[[nm]]
}

message("scDblFinder results:")
print(table(sce\$scDblFinder.class))

write_h5ad(sce, "${prefix}.h5ad", compression = "gzip")

write.csv(as.data.frame(colData(sce))[, cols, drop = FALSE], "${prefix}.scores.csv")

# One boolean column (TRUE = doublet), named after the prefix and indexed by cell
# barcode. This is the same layout the other doublet detection modules return, so
# the calls from several methods can be combined without any conversion.
calls <- data.frame(sce\$scDblFinder.class == "doublet", row.names = colnames(sce))
colnames(calls) <- "${prefix}"
write.csv(calls, "${prefix}.predictions.csv")

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    scdblfinder:', as.character(packageVersion('scDblFinder'))),
        paste('    anndatar:', as.character(packageVersion('anndataR'))),
        paste('    singlecellexperiment:', as.character(packageVersion('SingleCellExperiment'))),
        paste('    biocparallel:', as.character(packageVersion('BiocParallel'))),
        paste('    rhdf5:', as.character(packageVersion('rhdf5')))
    ),
    'versions.yml'
)
