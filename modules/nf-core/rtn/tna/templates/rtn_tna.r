#!/usr/bin/env Rscript

# Ported to nf-core/modules with template by Jonathan Manning

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

# Load
library("RTN")
library("snow")

################################################
################################################
## Pull in module inputs                      ##
################################################
################################################

tni_object <- readRDS('tni.rds')
# degs_log2 <- read_tsv
# degs <- read_tsv
# degs_annotation <-

output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
threads <- $task.cpus
args_opt <- parse_args('$task.ext.args')

n_perm <- ifelse('n_permutations' %in% names(args_opt), strtoi(args_opt[['n_permutations']]), 10)

# Debug messages (stdout)
sink(stdout(), type = "message") # sink messages to stdout
message("Expression matrix file   : ", input_expr_matrix)
message("Nb permutations          : ", n_perm)
message("Nb threads               : ", threads)
message("Output basename          : ", output_prefix)
if ('tfs' %in% names(args_opt)) {
    message("TFs                      : ", args_opt[['tfs']])
    tfs <- strsplit(args_opt[['tfs']], ',')
} else {
    # Load data
    data(tfsData)
    tfs <- tfsData\$Lambert2018\$SYMBOL
}
sink(NULL, type="message") # close the sink

# Input 1: 'object', a TNI object with regulons
# Input 2: 'phenotype', a named numeric vector, usually log2 differential expression levels
# Input 3: 'hits', a character vector, usually a set of differentially expressed genes
# Input 4: 'phenoIDs', an optional data frame with gene anottation mapped to the phenotype
rtna <- tni2tna.preprocess(object = tni_object, phenotype = degs_log2, hits = degs, phenoIDs = degs_annotation)

# TNA analysis here

saveRDS(rtni, file = "tni.rds")
saveRDS(rtni_permutation, file = "tni_permutated.rds")
saveRDS(rtni_bootstrapped, file = "tni_bootstrapped.rds")
saveRDS(rtni_filtered, file = "tni_filtered.rds")

# Plot
#pdf(paste0(output_prefix, "_RTN.pdf"))
#tni.graph(rtni_filtered, regulatoryElements = c("FOXM1", "E2F2"))
#title("Regulatory Transcriptional Network")
#mtext(output_prefix, side=3)
#dev.off()
#cat(
#    paste("- Threads::", threads),
#    fill=TRUE, labels=output_prefix,
#    file=paste0(output_prefix, "_intercept_slope.txt"), append=FALSE
#)

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
rtn.version <- as.character(packageVersion('RTN'))

writeLines(
    c(
        '"${task.process}":',
        paste('    bioconductor-rtn:', rtn.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
