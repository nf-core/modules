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

################################################
################################################
## Pull in module inputs                      ##
################################################
################################################

input_expr_matrix <- '$expression_matrix'
output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
threads <- $task.cpus
args_opt <- parse_args('$task.ext.args')

n_perm <- ifelse('n_permutations' %in% args_opt), args_opt[['n_permutations']], 10)

# Debug messages (stderr)
message("Expression matrix file   : ", input_expr_matrix)
message("Nb permutations          : ", n_perm)
message("Nb threads               : ", threads)
message("Output basename          : ", output_prefix)

# Load / install packages
library("RTN")
library("snow")
# Load data
data(tfsData)

# Preprocess
# Input 1: 'expData', a named gene expression matrix (genes on rows, samples on cols);
# Input 2: 'regulatoryElements', a vector listing genes regarded as TFs
# Input 3: 'rowAnnotation', an optional data frame with gene annotation
# Input 4: 'colAnnotation', an optional data frame with sample annotation
tfs <- tfsData$Lambert2018$SYMBOL
# tfs <- c("ENSMUSG00000001517", "ENSMUSG00000018983", "ENSMUSG00000016477",
#         "ENSMUSG00000039153", "ENSMUSG00000020415")

exp_data <- read.csv(expression_matrix,
                     sep='\t')
rownames(exp_data) <- exp_data$Geneid
rowAnnotation <- exp_data[,1:2]
colnames(rowAnnotation) <- c('PROBEID', 'SYMBOL')
rowAnnotation$SYMBOL <- toupper(rowAnnotation$SYMBOL)
exp_data[,1:2] <- NULL

# Regulatory Transcriptional Network Inference
rtni <- tni.constructor(expData = as.matrix(exp_data),
                        regulatoryElements = tfs,
                        rowAnnotation = rowAnnotation)

#options(cluster=snow::makeCluster(spec=threads, "SOCK"))

# Please set nPermutations >= 1000
#rtni_permutation <- tni.permutation(rtni, nPermutations = n_perm, pValueCutoff = 1e-7)

# Unstable interactions are subsequently removed by bootstrap analysis using the
# tni.bootstrap() function, which creates a consensus bootstrap network, referred
# here as refnet (reference network).
#rtni_bootstrapped <- tni.bootstrap(rtni_permutation)

#stopCluster(getOption("cluster"))

# remove the weakest interaction in any triplet formed by two TFs and a common
# target gene, preserving the dominant TF-target pair (ARACNe)
rtni_filtered <- tni.dpi.filter(rtni_bootstrapped)

# Summary of the resulting regulatory network
tni.regulon.summary(rtni, regulatoryElements = c("RUNX2"))
# tni.regulon.summary(rtni, regulatoryElements = "ENSMUSG00000001517")

# Plot
pdf(paste0(output_prefix, "_RTN.pdf"))
tni.graph(rtni, regulatoryElements = c("AEBP2", "AKAP8"))
title("Regulatory Transcriptional Network")
mtext(output_prefix, side=3)
dev.off()
cat(
    paste("- Threads::", threads),
    fill=TRUE, labels=output_prefix,
    file=paste0(output_prefix, "_intercept_slope.txt"), append=FALSE
)

write(line,file=paste0(output_prefix, "_duprateExpDensCurve_mqc.txt"),append=TRUE)
write.table(
    cbind(curve_x, curve_y),
    file=paste0(output_prefix, "_duprateExpDensCurve_mqc.txt"),
    quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE,
)

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
