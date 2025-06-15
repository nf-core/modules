#!/usr/bin/env Rscript

################################################
################################################
## USE PARAMETERS FROM NEXTFLOW               ##
################################################
################################################

# check fileHto

parser <- ArgumentParser("Parameters for BFF")
parser$add_argument("--fileHto", help="Path to file HTO count matrix.")
parser$add_argument("--methods", help='A vector of one or more calling methods to use.', default="combined_bff")
parser$add_argument("--methodsForConsensus", help='By default, a consensus call will be generated using all methods', default=NULL)
parser$add_argument("--cellbarcodeWhitelist", help='A vector of expected cell barcodes. This allows reporting on the total set of expected barcodes, not just those in the filtered count matrix',default=NULL)
parser$add_argument("--metricsFile", help='If provided, summary metrics will be written to this file.', default="metrics_cell_hash_r.csv")
parser$add_argument("--doTSNE", help='If true, tSNE will be run on the resulting hashing calls after each caller.', default=TRUE)
parser$add_argument("--preprocess_bff", help='If given, the PreProcess function by CellHashR is executed', action="store_true")
parser$add_argument("--barcodeWhitelist", help='A vector of barcode names to retain, used for preprocessing step', default=NULL)
parser$add_argument("--doHeatmap", help='f true, Seurat::HTOHeatmap will be run on the results of each caller.', default=TRUE)
parser$add_argument("--perCellSaturation", help='An optional dataframe with the columns cellbarcode and saturation',default=NULL)
parser$add_argument("--majorityConsensusThreshold", help='This applies to calculating a consensus call when multiple algorithms are used',default=NULL)
parser$add_argument("--chemistry", help='This string is passed to EstimateMultipletRate. Should be either 10xV2 or 10xV3.', default="10xV3")
parser$add_argument("--callerDisagreementThreshold", help='If provided, the agreement rate will be calculated between each caller and the simple majority call, ignoring discordant and no-call cells.',default=NULL)



# cast parameters from nextflow
hto_matrix = '$hto_matrix'
methods = '$methods' #COMBINED, RAW, CLUSTER
methodsForConsensus = 'NULL' #RAW, CLUSTER


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
