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

input_expr_matrix <- '$expression_matrix'
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

# Preprocess
# Input 1: 'expData', a named gene expression matrix (genes on rows, samples on cols);
# Input 2: 'regulatoryElements', a vector listing genes regarded as TFs
# Input 3: 'rowAnnotation', an optional data frame with gene annotation
# Input 4: 'colAnnotation', an optional data frame with sample annotation

exp_data <- read.csv(input_expr_matrix, sep='\t')
rownames(exp_data) <- exp_data[,1]
rowAnnotation <- exp_data[,1:2]
colnames(rowAnnotation) <- c('PROBEID', 'SYMBOL')
rowAnnotation\$SYMBOL <- toupper(rowAnnotation\$SYMBOL)
exp_data[,1:2] <- NULL

# Regulatory Transcriptional Network Inference
tfs <- c('ENSG00000125798', 'ENSG00000125816')
rtni <- tni.constructor(expData = as.matrix(exp_data),
                        regulatoryElements = tfs,
                        rowAnnotation = rowAnnotation)

options(cluster=snow::makeCluster(spec=threads, "SOCK"))

# Please set nPermutations >= 1000
rtni_permutation <- tni.permutation(rtni, nPermutations = n_perm, pValueCutoff = 1e-4)

# Unstable interactions are subsequently removed by bootstrap analysis using the
# tni.bootstrap() function, which creates a consensus bootstrap network, referred
# here as refnet (reference network).
rtni_bootstrapped <- tni.bootstrap(rtni_permutation)

stopCluster(getOption("cluster"))

# remove the weakest interaction in any triplet formed by two TFs and a common
# target gene, preserving the dominant TF-target pair (ARACNe)
rtni_filtered <- tni.dpi.filter(rtni_bootstrapped)

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
