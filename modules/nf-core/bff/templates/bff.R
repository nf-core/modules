#!/usr/bin/env Rscript

################################################
################################################
## Fucntions                                  ##
################################################
################################################

# Helper function for NULL condition
string_to_null <- function(x, val = "NULL") if (x == val) NULL else x
null_to_string <- function(x, val = "NULL") if (is.null(x)) val else x

################################################
################################################
## USE PARAMETERS FROM NEXTFLOW               ##
################################################
################################################

# cast parameters from nextflow

# input parameters
hto_matrix <- '$hto_matrix'
methods <- '$methods' # COMBINED, RAW, CLUSTER
preprocessing <- '$preprocessing'

# ext.args
methodsForConsensus <- string_to_null('$methodsForConsensus')  # NULL, RAW, CLUSTER
cellbarcodeWhitelist <- string_to_null('$cellbarcodeWhitelist')
prefix <- '$prefix'
metricsFile <- paste0(prefix, "_metrics_bff.csv")
doTSNE <- as.logical('$doTSNE')
barcodeWhitelist <- string_to_null('$barcodeWhitelist')
doHeatmap <- as.logical('$doHeatmap')
perCellSaturation <- string_to_null('$perCellSaturation')
majorityConsensusThreshold <- string_to_null('$majorityConsensusThreshold')
chemistry <- '$chemistry'
callerDisagreementThreshold <- string_to_null('$callerDisagreementThreshold')

# check if the file exists
if (! file.exists(hto_matrix)){
    stop(paste0(hto_matrix, ' is not a valid file'))
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(Seurat) # for Read10X()
library(cellhashR) # for ProcessCountMatrix() and GenerateCellHashingCalls()

################################################
################################################
## Main Process                               ##
################################################
################################################

# perform a basic preprocessing or direcly read the data
if(preprocessing){

  if(is.null(barcodeWhitelist)){
    barcodes_list <- NULL
  }else{
    # separate the barcodes by comma and remove leading/trailing whitespace from each word
    barcodes_list <- trimws(strsplit(barcodeWhitelist, ",")[[1]])
  }

  # perform preprocessing
  counts <- ProcessCountMatrix(rawCountData = hto_matrix, barcodeWhitelist = barcodes_list)
}else{
    # perform preprocessing
  counts <- Read10X(hto_matrix)
}

# determine the methods that should be used
methods_input <- switch(
  methods,
  "COMBINED" = c("bff_raw", "bff_cluster"),
  "RAW" = c("bff_raw"),
  "CLUSTER" = c("bff_cluster"),
  stop("'methods' must be either 'RAW', 'CLUSTER', or 'COMBINED'")
)

# determine the final consensus call based on RAW, CLUSTER or BOTH (NULL)
methodsForConsensus_input = NULL
if (!is.null(methodsForConsensus)) {
    methodsForConsensus_input <- switch(methodsForConsensus,
        "RAW" = "bff_raw",
        "CLUSTER" = "bff_cluster",
        stop("'methodsForConsensus' must be either 'RAW', 'CLUSTER', or NULL")
    )
}

results <- GenerateCellHashingCalls(barcodeMatrix = counts, methods = methods_input, doTSNE = doTSNE, doHeatmap = doHeatmap, methodsForConsensus = methodsForConsensus_input,cellbarcodeWhitelist = cellbarcodeWhitelist,metricsFile = metricsFile, perCellSaturation = perCellSaturation, majorityConsensusThreshold = majorityConsensusThreshold, chemistry = chemistry, callerDisagreementThreshold = callerDisagreementThreshold)

################################################
################################################
## SAVING RESULTS                             ##
################################################
################################################

#saving parameters in a dataframe
Argument <- c("hto_matrix", "methods", "methodsForConsensus", "cellbarcodeWhitelist", "metricsFile", "perCellSaturation","majorityConsensusThreshold","callerDisagreementThreshold", "doTSNE","doHeatmap","chemistry")
Value <- c(hto_matrix, null_to_string(methods), null_to_string(methodsForConsensus), null_to_string(cellbarcodeWhitelist), metricsFile, null_to_string(perCellSaturation), null_to_string(majorityConsensusThreshold), null_to_string(callerDisagreementThreshold), doTSNE, doHeatmap, chemistry)
params <- data.frame(Argument, Value)
write.csv(params, paste0(prefix ,"_params_bff.csv"))

# save the results of GenerateCellHashingCalls or an empty dataframe if NULL
if(is.null(results)){
    results <- data.frame()
}
write.csv(results, paste0(prefix, "_assignment_bff.csv"), row.names=FALSE)

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- paste(R.version[['major']],R.version[['minor']], sep = ".")
seurat.version <- as.character(packageVersion('Seurat'))
cellhashR.version <- as.character(packageVersion('cellhashR'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    r-seurat:', seurat.version),
        paste('    cellhashR:', cellhashR.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
