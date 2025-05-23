#!/usr/bin/env Rscript

# Ported to nf-core/modules with template by Jonathan Manning

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list,
                        function(x) scan(text=x,
                                         what='character',
                                         quiet = TRUE
                                         )
                       )

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals,
                             function(x) x[2]),
                             names = lapply(args_vals,
                                            function(x) x[1]
                                            )
                             )
    parsed_args[! is.na(parsed_args)]
}

# Load
library("RTN")
library("snow")

output_prefix = ifelse('$task.ext.prefix' == 'null',
                       '$meta.id',
                       '$task.ext.prefix')
threads <- $task.cpus

################################################
################################################
## Handling parameters                        ##
################################################
################################################

args_opt <- parse_args('$task.ext.args')

n_perm <- ifelse('n_permutations' %in% names(args_opt),
                 strtoi(args_opt[['n_permutations']]),
                 10)
p_value_cutoff <- ifelse('p_value_cutoff' %in% names(args_opt),
                         as.numeric(args_opt[['p_value_cutoff']]),
                         as.numeric(1e-4))

plot <- ifelse('plot' %in% names(args_opt),
               TRUE,
               FALSE)

################################################
################################################
## Pull in module inputs                      ##
################################################
################################################

if ('$geneAnnotation' != '') {
  rowAnnotation <- read.csv('$geneAnnotation', sep='\t', check.names = FALSE)
  cat("Using optional gene annotation information\n")
} else {
  rowAnnotation <- NULL
}

if ('$sampleAnnotation' != '') {
  colAnnotation <- read.csv('$sampleAnnotation', sep="\t", check.names = FALSE)
  cat("Using optional column annotation information\n")
} else {
  colAnnotation <- NULL
}

exp_data <- read.csv('$expression_matrix', sep='\t', check.names = FALSE)
tfs <- read.delim('$tfs', header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)[[1]]

# Debug messages (stdout)
sink(stdout(), type = "message") # sink messages to stdout
message("Expression matrix file   : ", '$expression_matrix')
message("TF file                  : ", '$tfs')
message("Nb permutations          : ", n_perm)
message("P-value cut off          : ", p_value_cutoff)
message("Nb threads               : ", threads)
message("Output basename          : ", output_prefix)
sink(NULL, type="message") # close the sink

# Regulatory Transcriptional Network Inference

# Preprocess
# Input 1: 'expData', a named gene expression matrix (genes as rownames, samples as colnames);
# Input 2: 'regulatoryElements', a vector listing genes regarded as TFs. Gene IDs must be in the same format as in expData
# Input 3: 'rowAnnotation', an optional data frame with gene annotation
# Input 4: 'colAnnotation', an optional data frame with sample annotation

rtni <- tni.constructor(expData = as.matrix(exp_data),
                        regulatoryElements = tfs,
                        rowAnnotation = rowAnnotation,
                        colAnnotation = colAnnotation
                        )

options(cluster=snow::makeCluster(spec=threads, "SOCK"))

# Please set nPermutations >= 1000
rtni_permutation <- tni.permutation(rtni,
                                    nPermutations = n_perm,
                                    pValueCutoff = p_value_cutoff
                                    )

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
if (plot) {
  pdf(paste0(output_prefix, "_RTN.pdf"))
  tni.graph(rtni_filtered, regulatoryElements = tfs)
  title("Regulatory Transcriptional Network")
  mtext(output_prefix, side=3)
  dev.off()
  cat(
      paste("- Threads::", threads),
      fill=TRUE, labels=output_prefix,
      file=paste0(output_prefix, "_intercept_slope.txt"), append=FALSE
  )
}

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
