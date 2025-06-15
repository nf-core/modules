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
kfunc = 'clara'
nsamples = as.integer('$nsamples')
seed = as.integer('$seed')
verbose = as.logical('$verbose')
prefix = '$prefix'

# check if the file exists
if (! file.exists(seuratObj)){
    stop(paste0(seuratObj, ' is not a valid file'))
}

## check here

# Set defaults and classes

opt <- list(
    seuratObj = '$seurat_object',
    quantile = 0.7,
    autoThresh = TRUE,
    maxiter = 5,
    qrangeFrom = 0.1,
    qrangeTo = 0.9,
    qrangeBy = 0.05,
    verbose = TRUE,
    assay = 'HTO',
    prefix = '$prefix'
)
opt_types <- lapply(opt, class)

# apply parameter overrides
args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }else{

        # Preserve classes from defaults where possible
        if (! is.null(opt[[ao]])){
            args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        }
        opt[[ao]] <- args_opt[[ao]]
    }
}

# Check if required parameters have been provided
required_opts <- c('seuratObj')
missing <- required_opts[!unlist(lapply(opt[required_opts], is_valid_string)) | !required_opts %in% names(opt)]
if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# check if the input path is a valid string
if (! is_valid_string(opt[['seuratObj']])) {
    stop(paste("Please provide", file_input), call. = FALSE)
}

# check if the file exists
if (! file.exists(opt[['seuratObj']])){
    stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
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
hashtag <- readRDS(opt\$seuratObj)

# Demultiplex cells
if (opt\$autoThresh == TRUE) {
  hashtag <- MULTIseqDemux(hashtag, assay = opt\$assay, quantile = opt\$quantile, autoThresh = TRUE, maxiter = opt\$maxiter, qrange = seq(from = opt\$qrangeFrom, to = opt\$qrangeTo, by = opt\$qrangeBy), verbose = opt\$verbose)
} else {
  hashtag <- MULTIseqDemux(hashtag, assay = opt\$assay, quantile = opt\$quantile, verbose = opt\$verbose)
}

################################################
################################################
## SAVING RESULTS                             ##
################################################
################################################

# create a data frame to save the used parameters in a csv file
Argument <- c("seuratObjectPath", "quantile", "autoThresh", "maxiter", "qrangeFrom", "qrangeTo", "qrangeBy", "verbose", "assay")
Value <- c(opt\$seuratObj, opt\$quantile, opt\$autoThresh, opt\$maxiter, opt\$qrangeFrom, opt\$qrangeTo, opt\$qrangeBy, opt\$verbose, opt\$assay)
params <- data.frame(Argument, Value)
write.csv(params, paste0(opt\$prefix ,"_params_multiseqdemux.csv"))

# save the results from MULTIseqDemux()
write.csv(hashtag\$MULTI_ID, paste0(opt\$prefix , "_res_multiseqdemux.csv"))
saveRDS(hashtag, file = paste0(opt\$prefix ,"_multiseqdemux.rds"))

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
