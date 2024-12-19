#!/usr/bin/env Rscript

library(scds)
library(SingleCellExperiment)

set.seed(0)

sce <- readRDS("${rds}")

# Annotate doublet using binary classification based doublet scoring:
sce <- bcds(sce, retRes = TRUE, estNdbl=TRUE)

# Annotate doublet using co-expression based doublet scoring:
try({
    sce <- cxds(sce, retRes = TRUE, estNdbl=TRUE)
})

# If cxds worked, run hybrid, otherwise use bcds annotations
if ("cxds_score" %in% colnames(colData(sce))) {
    # Combine both annotations into a hybrid annotation
    sce <- cxds_bcds_hybrid(sce, estNdbl=TRUE)

    predictions <- colData(sce)[, 'hybrid_call', drop=FALSE]
} else {
    predictions <- colData(sce)[, 'bcds_call', drop=FALSE]
}

saveRDS(sce, "${prefix}.rds")

colnames(predictions) <- "${prefix}"
write.csv(predictions, "${prefix}.csv")


# versions file
r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
scds.version <- as.character(packageVersion('scds'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r:', paste(version\$major, version\$minor, sep = ".")),
        paste('    scds:', as.character(packageVersion('scds'))),
        paste('    singlecellexperiment:', as.character(packageVersion('SingleCellExperiment')))
    ),
'versions.yml')
