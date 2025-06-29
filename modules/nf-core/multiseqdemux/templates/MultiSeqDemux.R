#!/usr/bin/env Rscript

################################################
################################################
## USE PARAMETERS FROM NEXTFLOW               ##
################################################
################################################

# cast parameters from nextflow
seuratObj = '$seurat_object'
quantile = as.numeric('$quantile')
autoThresh = as.logical('$autoThresh')
maxiter = as.numeric('$maxiter')
qrangeFrom = as.numeric('$qrangeFrom')
qrangeTo = as.numeric('$qrangeTo')
qrangeBy = as.numeric('$qrangeBy')
verbose = as.logical('$verbose')
assay ='$assay'
prefix = '$prefix'

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

# Demultiplex cells
if (autoThresh == TRUE) {
  hashtag <- MULTIseqDemux(hashtag, assay = assay, quantile = quantile, autoThresh = TRUE, maxiter = maxiter, qrange = seq(from = qrangeFrom, to = qrangeTo, by = qrangeBy), verbose = verbose)
} else {
  hashtag <- MULTIseqDemux(hashtag, assay = assay, quantile = quantile, verbose = verbose)
}

################################################
################################################
## SAVING RESULTS                             ##
################################################
################################################

# create a data frame to save the used parameters in a csv file
Argument <- c("seuratObjectPath", "quantile", "autoThresh", "maxiter", "qrangeFrom", "qrangeTo", "qrangeBy", "verbose", "assay")
Value <- c(seuratObj, quantile, autoThresh, maxiter, qrangeFrom, qrangeTo, qrangeBy, verbose, assay)
params <- data.frame(Argument, Value)
write.csv(params, paste0(prefix ,"_params_multiseqdemux.csv"))

# save the results from MULTIseqDemux()
write.csv(hashtag\$MULTI_ID, paste0(prefix , "_res_multiseqdemux.csv"))
saveRDS(hashtag, file = paste0(prefix ,"_multiseqdemux.rds"))

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
        paste('    r-seurat:', seurat.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
