#!/usr/bin/env Rscript

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}


################################################
################################################
## Pull in module inputs                      ##
################################################
################################################

inputFile <- '$fragments'
# is this a vector of input files?
sampleNames <- gsub(".fragments.tsv.gz","", gsub(".*/","", inputFile))
sampleName <- "$meta.id"
# get sample names
genome <- "$meta.genome"
# "hg19","hg38","mm9", and "mm10" currently supported
threads <- $task.cpus


library(ArchR)

addArchRGenome(genome)
addArchRThreads(threads = threads)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFile,
  sampleNames = sampleNames,
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# need to make minTSS and minFrags variables?


################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(output_prefix, "R_sessionInfo.log", sep = '.'))
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
archr.version <- as.character(packageVersion('ArchR'))

writeLines(
    c(
        '"${task.process}":',
        paste('    archR:', archr.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
