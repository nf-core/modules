#!/usr/bin/env Rscript

################################################
################################################
<<<<<<< HEAD
<<<<<<< HEAD
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

=======
## Functions                                  ##
=======
## USE PARAMETERS FROM NEXTFLOW               ##
>>>>>>> da0d66277 (adopted the feedback from the review)
################################################
################################################

# cast parameters from nextflow
seuratObj = '$seurat_object'
assay = '$assay'
options(digits=5)
quantile = as.double('$quantile')
init = '$init'
nstarts = as.integer('$nstarts')
kfunc = 'clara'
nsamples = as.integer('$nsamples')
seed = as.integer('$seed')
verbose = as.logical('$verbose')
prefix = '$prefix'

# check if the file exists
if (! file.exists(seuratObj)){
    stop(paste0(seuratObj, ' is not a valid file'))
}

<<<<<<< HEAD

>>>>>>> df971f6e6 (add template)
=======
>>>>>>> da0d66277 (adopted the feedback from the review)
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
<<<<<<< HEAD
<<<<<<< HEAD
hashtag <- readRDS(seuratObj)

# Demultiplex cells based on HTO enrichment
hashtag <- HTODemux(hashtag, assay = assay, positive.quantile = quantile, init = init, nstarts = nstarts, kfunc = kfunc, seed = seed, verbose = verbose)

=======
hashtag <- readRDS(opt\$seuratObj)
=======
hashtag <- readRDS(seuratObj)
>>>>>>> da0d66277 (adopted the feedback from the review)

# Demultiplex cells based on HTO enrichment
# zu einem int, boolean casten
if (kfunc == "clara") {
    hashtag <- HTODemux(hashtag, assay = assay, positive.quantile = quantile, init = init, nstarts = nstarts, kfunc = "clara", seed = seed, verbose = verbose)
} else {
    hashtag <- HTODemux(hashtag, assay = assay, positive.quantile = quantile, init = init, nstarts = nstarts, kfunc = "kmeans", seed = seed, verbose = verbose)
}
>>>>>>> df971f6e6 (add template)

################################################
################################################
## SAVING RESULTS                             ##
################################################
################################################

# create a data frame to save the used parameters in a csv file
<<<<<<< HEAD
<<<<<<< HEAD
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
=======
if (is.null(opt\$init)) {
=======
if (is.null(init)) {
>>>>>>> da0d66277 (adopted the feedback from the review)
  init <- "NULL"
}

Argument <- c("seuratObject", "quantile", "kfunc", "nstarts", "nsamples", "seed", "init", "assay", "verbose")
Value <- c(seuratObj, quantile, kfunc, nstarts, nsamples, seed, init, assay, verbose)
params <- data.frame(Argument, Value)

write.csv(params, paste0(prefix ,"_params_htodemux.csv"))

# create csv files to save the results from HTODemux()
<<<<<<< HEAD
donors <- rownames(hashtag[[opt\$assay]])
assignment <- hashtag[[paste0(opt\$assay, "_classification")]]
assignment[[paste0(opt\$assay, "_classification")]][!assignment[[paste0(opt\$assay, "_classification")]] %in% c(donors, "Negative")] <- "Doublet"
write.csv(assignment, paste0(opt\$prefix ,"_assignment_htodemux.csv"))
write.csv(hashtag[[paste0(opt\$assay, "_classification.global")]], paste0(opt\$prefix ,"_classification_htodemux.csv"))
<<<<<<< HEAD
saveRDS(hashtag, file = paste0(opt\$prefix ,"htodemux.rds"))
>>>>>>> df971f6e6 (add template)
=======
saveRDS(hashtag, file = paste0(opt\$prefix ,"_htodemux.rds"))
>>>>>>> 37633daa7 (add stub)
=======
donors <- rownames(hashtag[[assay]])
assignment <- hashtag[[paste0(assay, "_classification")]]
assignment[[paste0(assay, "_classification")]][!assignment[[paste0(assay, "_classification")]] %in% c(donors, "Negative")] <- "Doublet"
write.csv(assignment, paste0(prefix ,"_assignment_htodemux.csv"))
write.csv(hashtag[[paste0(assay, "_classification.global")]], paste0(prefix ,"_classification_htodemux.csv"))
saveRDS(hashtag, file = paste0(prefix ,"_htodemux.rds"))
>>>>>>> da0d66277 (adopted the feedback from the review)

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
<<<<<<< HEAD
<<<<<<< HEAD
        paste('    seurat:', seurat.version)
=======
        paste('    bioconductor-deseq2:', seurat.version)
>>>>>>> df971f6e6 (add template)
=======
        paste('    seurat:', seurat.version)
>>>>>>> 65d5301bf (add verbose and params to stub)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
