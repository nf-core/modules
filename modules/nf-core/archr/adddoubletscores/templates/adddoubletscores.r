################################################
################################################
## Pull in module inputs                      ##
################################################
################################################

inputFile <- "$fragments"
sampleName <- "$meta.id"
# get sample names
genome <- "$genome"
# "hg19","hg38","mm9", and "mm10" currently supported
threads <- "$task.cpus"


library(ArchR)

addArchRGenome(genome)
addArchRThreads(threads = threads)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# need to make minTSS and minFrags variables?

saveRDS(ArrowFiles, file = "arrowfiles.rds")

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(output_prefix, "R_sessionInfo.log", sep = "."))
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[["version.string"]], " ")[[1]][3]
archr.version <- as.character(packageVersion("ArchR"))

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