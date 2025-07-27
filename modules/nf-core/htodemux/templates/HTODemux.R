#!/usr/bin/env Rscript

################################################
################################################
## USE PARAMETERS FROM NEXTFLOW               ##
################################################
################################################

# cast parameters from nextflow
seuratObj = '$seurat_object'
assay = '$assay'
options(digits=5)
quantile = as.double('$quantile')
init = NULL
if ('$init' != "NULL") {
    init = as.integer('$init')
}
nstarts = as.integer('$nstarts')
kfunc = '$kfunc'
nsamples = as.integer('$nsamples')
seed = as.integer('$seed')
verbose = as.logical('$verbose')
prefix = '$prefix'

# check if the file exists
if (! file.exists(seuratObj)){
    stop(paste0(seuratObj, ' is not a valid file'))
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(Seurat)

################################################
################################################
## Main Process                               ##
################################################
################################################

# Loading Seurat object
hashtag <- readRDS(seuratObj)

# Demultiplex cells based on HTO enrichment
hashtag <- HTODemux(hashtag, assay = assay, positive.quantile = quantile, init = init, nstarts = nstarts, kfunc = kfunc, seed = seed, verbose = verbose)


################################################
################################################
## SAVING RESULTS                             ##
################################################
################################################

# create a data frame to save the used parameters in a csv file
if (is.null(init)) {
  init <- "NULL"
}

Argument <- c("seuratObject", "quantile", "kfunc", "nstarts", "nsamples", "seed", "init", "assay", "verbose")
Value <- c(seuratObj, quantile, kfunc, nstarts, nsamples, seed, init, assay, verbose)
params <- data.frame(Argument, Value)

write.csv(params, paste0(prefix ,"_params_htodemux.csv"))

# create csv files to save the results from HTODemux()
donors <- rownames(hashtag[[assay]])
assignment <- hashtag[[paste0(assay, "_classification")]]
assignment[[paste0(assay, "_classification")]][!assignment[[paste0(assay, "_classification")]] %in% c(donors, "Negative")] <- "Doublet"
write.csv(assignment, paste0(prefix ,"_assignment_htodemux.csv"))
write.csv(hashtag[[paste0(assay, "_classification.global")]], paste0(prefix ,"_classification_htodemux.csv"))
saveRDS(hashtag, file = paste0(prefix ,"_htodemux.rds"))

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
seurat.version <- as.character(packageVersion('Seurat'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    seurat:', seurat.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
