#!/usr/bin/env Rscript

################################################
################################################
## Fucntions                                  ##
################################################
################################################

# Helper function for NULL condition
null_if <- function(x, val = "NULL") if (x == val) NULL else x

################################################
################################################
## USE PARAMETERS FROM NEXTFLOW               ##
################################################
################################################

# cast parameters from nextflow

# input parameters
hto_matrix <- '$hto_matrix'
methods <- '$methods' # COMBINED, RAW, CLUSTER
preprocessing <- '$preprocessing' #    preprocess_bff = task.ext.preprocess_bff ?: "FALSE"

# ext.args
methodsForConsensus <- null_if('$methodsForConsensus')  # NULL, RAW, CLUSTER
cellbarcodeWhitelist <- null_if('$cellbarcodeWhitelist')
prefix <- 'prefix'
metricsFile <- paste0(prefix, "_metrics_bff.csv")
doTSNE <- as.logical('$doTSNE')
barcodeWhitelist <- null_if('$barcodeWhitelist')
doHeatmap <- as.logical('$doHeatmap')
perCellSaturation <- null_if('$perCellSaturation')
majorityConsensusThreshold <- null_if('$majorityConsensusThreshold')
chemistry <- '$chemistry'
callerDisagreementThreshold <- null_if('$callerDisagreementThreshold')

# check if the file exists
if (! file.exists(hto_matrix)){
    stop(paste0(hto_matrix, ' is not a valid file'))
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(DropletUtils)
library(Seurat)
library(cellhashR)
library(here)
library(dplyr)
library(ggplot2)
library(argparse)

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
